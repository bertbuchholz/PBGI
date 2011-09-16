#ifndef CUBERASTERBUFFER_H
#define CUBERASTERBUFFER_H

#include <algorithm>
#include <vector>
#include <map>

#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <utilities/Debug_info.h>
#include <utilities/Frame_buffer.h>


__BEGIN_YAFRAY

template <class T>
inline T into_range(T const& min, T const& max, T const& val)
{
    if (val < min) return min;
    else if (val > max) return max;
    return val;
}

template <class T>
inline bool is_in_range(T const& min, T const& max, T const& val)
{
    if (val < min) return false;
    else if (val > max) return false;
    return true;
}

struct Gi_point_info
{
    enum Type { Node = 0, Far_surfel, Near_surfel };

    Type type;
    color_t color;
    vector3d_t position;
    vector3d_t receiver_position;
    float depth;
    vector3d_t direction; // from the receiving point to the GI point
    float solid_angle;
    vector3d_t disc_normal;
    float radius;
};

struct Cube_cell
{
    int plane;
    int pos[2];

    friend bool operator==(Cube_cell const& lhs, Cube_cell const& rhs)
    {
        return lhs.plane == rhs.plane && lhs.pos[0] == rhs.pos[0] && lhs.pos[1] == rhs.pos[1];
    }
};


/*

framebuffer types:
cube, hemisphere

- jittered, jitter each cube side by a random 2d vector
- stochastic splatting, splat with russian roulette depending on the alpha value


 splatting types:
 - the single pixel the point is projected onto (most simple)
 - trace ray through pixels intersecting with the disc (surfel only)
 - using the solid angle, approximating the area with axis-aligned rectangles

*/


class Cube_raster_buffer
{
public:
    enum Type { Simple = 0, Accumulating, Reweighting, Distance_weighted, Front_stacked };
    enum Splat_type { Single_pixel = 0, Disc_tracing, AA_square, SA_tracing };

    static std::map<std::string, Splat_type> enum_splat_type_map;
    static std::map<std::string, Type> enum_type_map;

    static const vector3d_t cube_plane_normals[];

    typedef vector3d_t Point;
    typedef color_t Color;

    Cube_raster_buffer()
    {
        add_point_map.resize(3);
    }

    Cube_raster_buffer & operator= (Cube_raster_buffer const& cbr)
    {
        if (this != &cbr)
        {
            for (int i = 0; i < 6; ++i)
            {
                buffers[i] = cbr.buffers[i]->clone();
            }

            total_energy  = cbr.total_energy;
            _resolution   = cbr._resolution;
            _resolution_2 = cbr._resolution_2;
            _cell_centers = cbr._cell_centers;
            add_point_map = cbr.add_point_map;
        }

        return *this;
    }

    ~Cube_raster_buffer()
    {
        for (int i = 0; i < 6; ++i)
        {
            delete buffers[i];
        }
    }

    void setup_simple(int resolution)
    {
        setup(Simple, resolution, Single_pixel, Single_pixel, Single_pixel);
    }

    void setup(Type fb_type,
               int resolution,
               Splat_type const node_splat_type,
               Splat_type const surfel_far_splat_type,
               Splat_type const surfel_near_splat_type)
    {
        _resolution = resolution;
        _resolution_2 = _resolution / 2;

        std::map<Splat_type, Add_point_function_ptr> splat_type_to_function_map;
        splat_type_to_function_map[Single_pixel] = &Cube_raster_buffer::add_point_single_pixel;
        splat_type_to_function_map[Disc_tracing] = &Cube_raster_buffer::add_point_disc_tracing;
        splat_type_to_function_map[AA_square]    = &Cube_raster_buffer::add_point_aa_square;
        splat_type_to_function_map[SA_tracing]   = &Cube_raster_buffer::add_point_solid_angle_rays;

        add_point_map[Gi_point_info::Node]        = splat_type_to_function_map[node_splat_type];
        add_point_map[Gi_point_info::Far_surfel]  = splat_type_to_function_map[surfel_far_splat_type];
        add_point_map[Gi_point_info::Near_surfel] = splat_type_to_function_map[surfel_near_splat_type];

        for (int i = 0; i < 6; ++i)
        {
            if (fb_type == Simple)
            {
                buffers[i] = new Simple_frame_buffer(_resolution, i);
            }
            else if (fb_type == Accumulating)
            {
                buffers[i] = new Accumulating_frame_buffer(_resolution, i);
            }
            else if (fb_type == Reweighting)
            {
                buffers[i] = new Reweighting_frame_buffer(_resolution, i);
            }
            else if (fb_type == Distance_weighted)
            {
                buffers[i] = new Distance_weighted_frame_buffer(_resolution, i);
            }
            else if (fb_type == Front_stacked)
            {
                buffers[i] = new Front_stacked_frame_buffer(_resolution, i);
            }

            else
            {
                std::cout << "Unknown FB type ..." << std::endl;
                assert(false);
            }
        }

        Cube_cell c;
        _cell_centers.resize(_resolution * _resolution * 6);

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            c.plane = plane_index;

            for (int u = -_resolution_2; u < _resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -_resolution_2; v < _resolution_2; ++v)
                {
                    c.pos[1] = v;

                    _cell_centers[get_serial_index(c)] = calc_cell_center(c);
                    _cube_cells.push_back(c);
                }
            }
        }
    }

    int get_resolution() const
    {
        return _resolution;
    }

    void clear()
    {
        for (int i = 0; i < 6; ++i)
        {
            buffers[i]->clear();
        }
    }


    inline int get_longest_axis(Point const& dir) const
    {
        int longest_axis = 0;
        if (std::abs(dir[1]) > std::abs(dir[0]) && std::abs(dir[1]) > std::abs(dir[2])) longest_axis = 1;
        else if (std::abs(dir[2]) > std::abs(dir[0])) longest_axis = 2;
        return longest_axis;
    }

    Cube_cell get_cell(Point const& dir) const
    {
        Cube_cell result;

        int longest_axis = get_longest_axis(dir);

        int plane_index = longest_axis * 2;
        int axis_0 = (longest_axis + 1) % 3;
        int axis_1 = (longest_axis + 2) % 3;
        float longest_extent = dir[longest_axis];

        if (dir[longest_axis] < 0)
        {
            ++plane_index;
            longest_extent = -longest_extent;
        }

        result.plane = plane_index;

        assert(longest_extent != 0.0f);

        float inv_longest_extent = 1.0f / longest_extent;

        float axis_0_pos = dir[axis_0] * inv_longest_extent;
        float axis_1_pos = dir[axis_1] * inv_longest_extent;

        result.pos[0] = std::floor(axis_0_pos * _resolution_2);
        result.pos[1] = std::floor(axis_1_pos * _resolution_2);

        result.pos[0] = std::max(-_resolution_2, std::min(result.pos[0], _resolution_2 - 1));
        result.pos[1] = std::max(-_resolution_2, std::min(result.pos[1], _resolution_2 - 1));

        return result;
    }

    bool get_cell(Point const& dir, Cube_cell & cell) const
    {
        int longest_axis = get_longest_axis(dir);

        int plane_index = longest_axis * 2;
        int axis_0 = (longest_axis + 1) % 3;
        int axis_1 = (longest_axis + 2) % 3;
        float longest_extent = dir[longest_axis];

        if (dir[longest_axis] < 0)
        {
            ++plane_index;
            longest_extent = -longest_extent;
        }

        cell.plane = plane_index;

        assert(longest_extent != 0.0f);

        float inv_longest_extent = 1.0f / longest_extent;

        float axis_0_pos = dir[axis_0] * inv_longest_extent;
        float axis_1_pos = dir[axis_1] * inv_longest_extent;

        cell.pos[0] = std::floor(axis_0_pos * _resolution_2);
        cell.pos[1] = std::floor(axis_1_pos * _resolution_2);

        if (cell.pos[0] < -_resolution_2 || cell.pos[0] >= _resolution_2) return false;
        if (cell.pos[1] < -_resolution_2 || cell.pos[1] >= _resolution_2) return false;

        return true;
    }


    int get_serial_index(Cube_cell const& c) const
    {
        int serial =
                (c.pos[1] + _resolution_2) +
                (c.pos[0] + _resolution_2) * _resolution +
                 c.plane * _resolution * _resolution;

        return serial;
    }

    Point const& get_cell_center(Cube_cell const& c) const
    {
        int serial = get_serial_index(c);

        return _cell_centers[serial];
    }

    std::vector<Cube_cell> const& get_cube_cells()
    {
        return _cube_cells;
    }

    std::vector<Point> get_cell_corners(Cube_cell const& c) const
    {
        std::vector<Point> result;

        Point cell_center = get_cell_center(c);

        float dirs[] = { -1, -1, 1, -1, 1, 1, -1, 1 };

        for (int i = 0; i < 4; ++i)
        {
            Cube_cell neighbor(c);

            neighbor.pos[0] += dirs[i * 2 + 0];
            neighbor.pos[1] += dirs[i * 2 + 1];

            Point neighbor_center = calc_cell_center(neighbor);

            result.push_back((cell_center + neighbor_center) / 2.0f);
        }

        return result;
    }


    // p0, p1: edge to test
    // v0: point on the plane, n: plane normal
    // alpha: return value, intersection point: p0 + alpha * (p1 - p0)
    bool linePlaneIntersection(Point const& p0, Point const& p1, Point const& v0, Point const& n, float & alpha) const
    {
        float d1 = dot(n, v0 - p0);
        float d2 = dot(n, p1 - p0);

        if (std::abs(d2) < 1e-5) return false;

        alpha = d1 / d2;

        return (alpha >= 0.0f && alpha <= 1.0f);
    }


    // dir: direction of line, always starting in the origin
    // n: plane normal
    // alpha: return value,
    // intersection point:  alpha * dir
    bool linePlaneIntersection(Point const& dir, Point const& n, float & alpha) const
    // bool linePlaneIntersection(Point const& dir, Point const& v0, Point const& n, float & alpha)
    {
        // float d1 = dot(n, v0);
        float const d1 = 1.0f;
        float const d2 = n * dir;

        if (std::abs(d2) < 1e-10) return false;

        alpha = d1 / d2;

        return (alpha >= 0.0f && alpha <= 1.0f);
    }

    // only check in the positive halfspace of the line, line starts at origin
    // source: http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    bool sphere_line_intersection(Point const& sphere_center, float const radius, Point const& line_dir)
    {
        float const c_l = sphere_center * line_dir;

        if (c_l <= 0.0f) return false;

        return (c_l * c_l) - sphere_center.lengthSqr() + (radius * radius) >= 0.0f;
    }


    float round(float const d) const
    {
        return std::floor(d + 0.5f);
    }

    float next_raster_in_positive_direction(float u) const
    {
        return std::floor(u + 1.0f);
    }

    void add_point_single_pixel(Gi_point_info const& point_info, GiPoint const* node = NULL)
    {
        Point dir = point_info.direction;
        Color const& color = point_info.color;
        float const depth = point_info.depth;
        float const radius = point_info.radius;
        dir.normalize();

        Cube_cell c;
        if (!get_cell(dir, c)) return;

        buffers[c.plane]->add_point(c.pos[0], c.pos[1], color, 1.0f, depth, radius, node);
    }

    // raytrace a disc through every pixel, assumes receiving point is in the origin and the disc's center relative to it
    void add_point_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL)
    {
        Color const& color = point_info.color;
        Point const& disc_normal = point_info.disc_normal;
        Point disc_center = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
        float const disc_radius = point_info.radius;
        // float const depth = point_info.depth;

        Cube_cell c;

        float disc_radius_sqr = disc_radius * disc_radius;

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            c.plane = plane_index;

            for (int u = -_resolution_2; u < _resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -_resolution_2; v < _resolution_2; ++v)
                {
                    c.pos[1] = v;

                    Point cell_center = get_cell_center(c);
                    cell_center.normalize();

                    float alpha = -1.0f;

                    linePlaneIntersection(Point(0.0f), cell_center, disc_center, disc_normal, alpha);

                    if (alpha > 0.0f)
                    {
                        Point hit_point = alpha * cell_center;

                        float dist_sqr = (hit_point - disc_center).lengthSqr();

                        if (dist_sqr < disc_radius_sqr)
                        {
                            buffers[plane_index]->add_point(u, v, color, 1.0f, alpha, disc_radius, node);
                        }
                    }
                }
            }
        }
    }



    void add_point_aa_square(Gi_point_info const& point_info, GiPoint const* node = NULL)
    {
        float const depth = point_info.depth;
        color_t const& color = point_info.color;
        float const solid_angle = point_info.solid_angle;
        vector3d_t const& direction = point_info.direction;
        float const radius = point_info.radius;

        if (std::abs(solid_angle) < 1e-10f)
        {
//            std::cout << "solid angle 0" << std::endl;
            return;
        }

        if (std::isnan(depth))
        {
            std::cout << "depth nan" << std::endl;
            return;
        }

        if (std::isnan(color.R) || std::isnan(color.G) || std::isnan(color.B))
        {
            std::cout << "color nan" << std::endl;
            return;
        }

        if (std::isnan(solid_angle))
        {
            std::cout << "solid angle nan" << std::endl;
            return;
        }

        if (std::isnan(direction[0]) ||
            std::isnan(direction[1]) ||
            std::isnan(direction[2])
           )
        {
            std::cout << "dir nan" << std::endl;
            return;
        }


        Point dir = direction;
        dir.normalize();

        // --- new square projection ---

        float const pixel_width = 1.0f / float(_resolution_2);

        // find intersection point on plane
        for (int i = 0; i < 3; ++i)
        {
            int normal_axis = i;

            if (std::abs(dir[normal_axis]) < 1e-10f) continue;

            int plane_index = normal_axis * 2;

            if (dir[normal_axis] < 0.0f) ++plane_index;

            assert(dir * cube_plane_normals[plane_index] > 0.0f);

            float alpha = 0.0f;
            linePlaneIntersection(dir, cube_plane_normals[plane_index], alpha);
            if (alpha < 1.0f)
            {
                std::cout << alpha << " " << dir << " " <<  cube_plane_normals[plane_index] << std::endl;
            }
            assert(alpha >= 1.0f);

            Point hit_point = alpha * dir;
            Point hit_dir = hit_point;
            hit_dir.normalize();


            float const square_area = solid_angle; // * (dir * normals[plane_index]);
            float const square_width_2 = std::sqrt(square_area) / 2.0f;

            int const axis_0 = (normal_axis + 1) % 3;
            int const axis_1 = (normal_axis + 2) % 3;

            float u_minus = (hit_point[axis_0] - square_width_2) * _resolution_2;
            float u_plus  = (hit_point[axis_0] + square_width_2) * _resolution_2;

            float v_minus = (hit_point[axis_1] - square_width_2) * _resolution_2;
            float v_plus  = (hit_point[axis_1] + square_width_2) * _resolution_2;

            u_minus = into_range(-float(_resolution_2), float(_resolution_2), u_minus);
            v_minus = into_range(-float(_resolution_2), float(_resolution_2), v_minus);

            u_plus = into_range(-float(_resolution_2), float(_resolution_2), u_plus);
            v_plus = into_range(-float(_resolution_2), float(_resolution_2), v_plus);

            float const u_center = hit_point[axis_0] * _resolution_2;
            float const v_center = hit_point[axis_1] * _resolution_2;

            float u = u_minus;

            float grid_area = 1.0f;
            // float grid_area = (1.0f / float(resolution)) * (1.0f / float(resolution)) * (resolution_2 * resolution_2);


            float debug_acc_area = 0.0f;

            while (u < u_plus)
            {
                float next_u = next_raster_in_positive_direction(u);
                if (next_u > u_plus) next_u = u_plus;

                float v = v_minus;

                if (std::abs(u - std::floor(u + 0.5f)) < 1e-5f) u = std::floor(u + 0.5f);

                while (v < v_plus)
                {
                    if (std::abs(v - std::floor(v + 0.5f)) < 1e-5f) v = std::floor(v + 0.5f);

                    float next_v = next_raster_in_positive_direction(v);
                    if (next_v > v_plus) next_v = v_plus;

                    //                std::cout << "u: " << u << " " << next_u << " v: " << v << " " << next_v << std::endl;

                    float const area = (next_u - u) * (next_v - v);

                    debug_acc_area += area;

                    float fill_ratio = area / grid_area;


                    int const cell_u = into_range(-_resolution_2, _resolution_2 - 1, int(std::floor(u)));
                    int const cell_v = into_range(-_resolution_2, _resolution_2 - 1, int(std::floor(v)));

                    /*
                    if (square_width_2 * 2.0f > pixel_width * 3.0f)
                    {
                        float const sigma = (square_width_2 * 2.0f) / 2.0f;
                        float const dist_from_center = std::sqrt((cell_u - u_center) * (cell_u - u_center) + (cell_v - v_center) * (cell_v - v_center));
                        float const x = dist_from_center;
                        float const gauss = std::min(1.0, 1.0f / std::sqrt(M_2PI * sigma * sigma) * std::exp(-(x*x) / (2.0f * sigma * sigma)));
                        fill_ratio *= gauss;
                    }
                    */

                    /*
                    Cube_cell c;
                    c.pos[0] = cell_u;
                    c.pos[1] = cell_v;
                    c.plane = plane_index;

                    vector3d_t cell_dir = get_cell_center(c);
                    cell_dir.normalize();

                    float const alpha = std::max(0.0f, 2.0f * ((dir * cell_dir) - 0.5f));
                    fill_ratio *= alpha;
                    */

                    // buffers[plane_index]->add_point(cell_u, cell_v, color, ratio, depth, radius, node);
                    buffers[plane_index]->add_point(cell_u, cell_v, color, fill_ratio, depth, radius, node);
                    // buffers[plane_index].add_point(cell_u, cell_v, color * inv_cell_area, ratio, depth);

                    v = next_v;
                }

                u = next_u;
            }
        }
    }


    void add_point_solid_angle_rays(Gi_point_info const& point_info, GiPoint const* node = NULL)
    {
        Color const& color = point_info.color;
        Point const& disc_normal = (point_info.receiver_position - point_info.position).normalize();
        Point disc_center = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
        float const disc_radius = std::sqrt(point_info.solid_angle * point_info.depth * point_info.depth / M_PI);
        float const radius = point_info.radius;
        // float const depth = point_info.depth;

        Cube_cell c;

        float disc_radius_sqr = disc_radius * disc_radius;

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            c.plane = plane_index;

            for (int u = -_resolution_2; u < _resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -_resolution_2; v < _resolution_2; ++v)
                {
                    c.pos[1] = v;

                    Point cell_center = get_cell_center(c);
                    cell_center.normalize();

                    float alpha;

                    linePlaneIntersection(Point(0.0f), cell_center, disc_center, disc_normal, alpha);

                    if (alpha > 0.0f)
                    {
                        Point hit_point = alpha * cell_center;

                        float dist_sqr = (hit_point - disc_center).lengthSqr();

                        if (dist_sqr < disc_radius_sqr)
                        {
                            buffers[plane_index]->add_point(u, v, color, 1.0f, alpha, radius, node);
                        }
                    }
                }
            }
        }


        /*
        float const solid_angle = point_info.solid_angle;
        Color const& color = point_info.color;
        Point sphere_center = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the sphere's center
        float const sphere_radius = std::sqrt(solid_angle * point_info.depth * point_info.depth / M_PI);
        float const depth = point_info.depth;

        Cube_cell c;

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            c.plane = plane_index;

            for (int u = -_resolution_2; u < _resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -_resolution_2; v < _resolution_2; ++v)
                {
                    c.pos[1] = v;

                    Point cell_center = get_cell_center(c);
                    cell_center.normalize();

                    if (sphere_line_intersection(sphere_center, sphere_radius, cell_center))
                    {
                        buffers[plane_index]->add_point(u, v, color, 1.0f, depth, node);
                    }
                }
            }
        }
        */
    }



    void add_point_solid_angle_rays_sub_pixel(Gi_point_info const& point_info, GiPoint const* node = NULL)
    {
        Color const& color = point_info.color;
        Point const& disc_normal = (point_info.receiver_position - point_info.position).normalize();
        Point disc_center = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
        float const disc_radius = std::sqrt(point_info.solid_angle * point_info.depth * point_info.depth / M_PI);
        float const radius = point_info.radius;
        // float const depth = point_info.depth;

        Cube_cell c;

        float disc_radius_sqr = disc_radius * disc_radius;

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            c.plane = plane_index;

            for (int u = -_resolution_2; u < _resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -_resolution_2; v < _resolution_2; ++v)
                {
                    c.pos[1] = v;

                    Point cell_center = get_cell_center(c);
                    cell_center.normalize();

                    float alpha;

                    linePlaneIntersection(Point(0.0f), cell_center, disc_center, disc_normal, alpha);

                    if (alpha > 0.0f)
                    {
                        Point hit_point = alpha * cell_center;

                        float dist_sqr = (hit_point - disc_center).lengthSqr();

                        if (dist_sqr < disc_radius_sqr)
                        {
                            buffers[plane_index]->add_point(u, v, color, 1.0f, alpha, radius, node);
                        }
                    }
                }
            }
        }
    }


    void add_point(Gi_point_info const& point_info, GiPoint const* gi_point = NULL)
    {
        (this->*add_point_map[point_info.type])(point_info, gi_point);
    }

    void set_color(Cube_cell const& c, Color const& color)
    {
        buffers[c.plane]->set_color(c.pos[0], c.pos[1], color);
    }

    void add_color(Cube_cell const& c, Color const& color)
    {
        Color current_color = get_color(c);
        set_color(c, current_color + color);
    }

    Color get_color(Cube_cell const& c, Debug_info * debug_info = NULL) const
    {
        return buffers[c.plane]->get_color(c.pos[0], c.pos[1], debug_info);
    }

    /*
    friend std::ostream & operator<<(std::ostream & s, Cube_raster_buffer const& c)
    {
        for (int i = 0; i < 6; ++i)
        {
            s << c.buffers[i] << " ";
        }

        return s;
    }

    friend std::istream & operator>>(std::istream & s, Cube_raster_buffer & c)
    {
        for (int i = 0; i < 6; ++i)
        {
            s >> c.buffers[i];
        }

        return s;
    }
    */

    float get_solid_angle(Cube_cell const& c) const
    {
        // TODO: add solid angle calculation for different cells depending on distance and angle
        float radius = (1.0f / float(_resolution));
        Point p = get_cell_center(c);
        p.normalize();

        // float corrected_area = 2.0f * M_PI * radius * radius * M_PI * 0.5f;

        float solid_angle = (cube_plane_normals[c.plane] * p) * radius * radius;

        return solid_angle;
    }

    Color get_diffuse(Point const& normal) const
    {
        Color diffuse(0.0f);

        Cube_cell c;

        float inv_cell_area = 1.0f / float(_resolution_2 * _resolution_2);

        for (int i = 0; i < 6; ++i)
        {
            c.plane = i;

            for (int x = -_resolution_2; x < _resolution_2; ++x)
            {
                c.pos[0] = x;

                for (int y = -_resolution_2; y < _resolution_2; ++y)
                {
                    c.pos[1] = y;

                    Point p = get_cell_center(c);
                    p.normalize();

                    float cos_sp_normal_cell_dir = std::max(0.0f, p * normal); // lambertian

                    if (cos_sp_normal_cell_dir < 0.001f) continue;

                    float cos_plane_normal_cell_dir = std::max(0.0f, cube_plane_normals[c.plane] * p); // solid angle of pixel

                    if (cos_plane_normal_cell_dir < 0.001f) continue;

                    Color const& color = get_color(c);

                    diffuse += color * cos_sp_normal_cell_dir * cos_plane_normal_cell_dir;
                    // diffuse += color;
                }
            }
        }

        return diffuse * inv_cell_area;
    }

    float get_total_energy()
    {
        return total_energy;
    }

    void get_debug_info(Debug_info & debug_info)
    {
        Cube_cell c;

        for (int i = 0; i < 6; ++i)
        {
            c.plane = i;

            for (int x = -_resolution_2; x < _resolution_2; ++x)
            {
                c.pos[0] = x;

                for (int y = -_resolution_2; y < _resolution_2; ++y)
                {
                    c.pos[1] = y;

                    get_color(c, &debug_info);
                }
            }
        }
    }

private:
    Point calc_cell_center(Cube_cell const& c) const
    {
        Point result;

        bool negative = (c.plane % 2 != 0);
        int plane = c.plane / 2;

        if (negative)
        {
            result[plane] = -1.0f;
        }
        else
        {
            result[plane] = 1.0f;
        }

        int axis_0 = (plane + 1) % 3;
        int axis_1 = (plane + 2) % 3;

        result[axis_0] = (c.pos[0] + 0.5f) / float(_resolution_2);
        result[axis_1] = (c.pos[1] + 0.5f) / float(_resolution_2);

        return result;
    }

    typedef void (Cube_raster_buffer::*Add_point_function_ptr)(Gi_point_info const& point_info, GiPoint const* gi_point);

    Abstract_frame_buffer* buffers[6];
    float total_energy;
    int _resolution;
    int _resolution_2;

    std::vector<Point> _cell_centers;
    std::vector<Cube_cell> _cube_cells;

    std::vector<Add_point_function_ptr> add_point_map;
};

__END_YAFRAY

#endif // CUBERASTERBUFFER_H

