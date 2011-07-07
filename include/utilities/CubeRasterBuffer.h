#ifndef CUBERASTERBUFFER_H
#define CUBERASTERBUFFER_H

#include <algorithm>
#include <vector>

#include <core_api/color.h>


__BEGIN_YAFRAY

class GiPoint;

template <class T>
inline T into_range(T const& min, T const& max, T const& val)
{
    if (val < min) return min;
    else if (val > max) return max;
    return val;
}

struct Cube_cell
{
    int plane;
    int pos[2];
};

struct Color_depth_pixel
{
    Color_depth_pixel() :
        depth(1e10f),
        filling_degree(0.0f),
        node(NULL)
    {

    }

    color_t color;
    float depth;
    float filling_degree;
    GiPoint const* node;

    friend bool operator<(Color_depth_pixel const& lhs, Color_depth_pixel const& rhs)
    {
        return (lhs.depth < rhs.depth);
    }

};

class Abstract_frame_buffer
{
public:
    Abstract_frame_buffer() {}

    Abstract_frame_buffer(int size) : _size(size)
    { }

    virtual ~Abstract_frame_buffer() {}

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, GiPoint const* node = NULL) = 0;

    virtual float accumulate(std::vector<GiPoint const*> * gi_points = NULL) = 0;

    virtual color_t const& get_color(int const x, int const y) const = 0;

    virtual void set_size(int size) = 0;

    virtual void clear() = 0;

    virtual Abstract_frame_buffer* clone() = 0;

protected:
    int _size;
};

class Simple_frame_buffer : public Abstract_frame_buffer
{
public:
    Simple_frame_buffer() {}

    Simple_frame_buffer(int size) : Abstract_frame_buffer(size)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Simple_frame_buffer* sfb = new Simple_frame_buffer(*this);

        return sfb;
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, GiPoint const* node = NULL)
    {
        int pixel = (x + _size / 2) + _size * (y + _size / 2);

        assert(pixel < _size * _size);

        if (depth < _data[pixel].depth)
        {
            _data[pixel].color = color;
            _data[pixel].depth = depth;
            _data[pixel].filling_degree = filling_degree;
            _data[pixel].node = node;
        }
    }

    virtual void set_size(int size)
    {
        _size = size;
        _data.clear();
        _data.resize(size * size);
    }

    virtual void clear()
    {
        _data.clear();
    }


    virtual float accumulate(std::vector<GiPoint const*> * gi_points = NULL)
    {
//        std::cout << "accumulate()" << std::endl;

        float total_energy = 0.0f;

        for (int pixel = 0; pixel < _size * _size; ++pixel)
        {
            Color_depth_pixel const& c = _data[pixel];
            total_energy += c.color.energy();

            if (gi_points && c.node)
            {
                gi_points->push_back(c.node);
            }
        }

        return total_energy;
    }

    virtual color_t const& get_color(int const x, int const y) const
    {
        int const pixel = (x + _size / 2) + _size * (y + _size / 2);
        return _data[pixel].color;
    }

protected:
    std::vector<Color_depth_pixel> _data;
};

class Accumulating_frame_buffer : public Abstract_frame_buffer
{
public:
    Accumulating_frame_buffer(int size) : Abstract_frame_buffer(size)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Accumulating_frame_buffer* afb = new Accumulating_frame_buffer(*this);

        return afb;
    }

    virtual void set_size(int size)
    {
        _size = size;

        _data.clear();
        _data.resize(size * size);

        _data_accumulated.clear();
        _data_accumulated.resize(size * size);
    }

    virtual void clear()
    {
        _data.clear();
        _data_accumulated.clear();
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, GiPoint const* node = NULL)
    {
        int pos = (x + _size / 2) + _size * (y + _size / 2);

        assert(pos < _size * _size);

        Color_depth_pixel c;
        c.color = color;
        c.depth = depth;
        c.filling_degree = filling_degree;
        c.node = node;

        _data[pos].push_back(c);
    }

    virtual float accumulate(std::vector<GiPoint const*> * gi_points = NULL)
    {
//        std::cout << "accumulate()" << std::endl;

        int samples_sum = 0;
        float total_energy = 0.0f;

        for (int pixel = 0; pixel < _size * _size; ++pixel)
        {
            std::sort(_data[pixel].begin(), _data[pixel].end());

            samples_sum += _data[pixel].size();

            float acc_filling_degree = 0.0f;

            unsigned int i = 0;

            _data_accumulated[pixel] = color_t(0.0f);

            for (unsigned int i = 0; i < _data[pixel].size(); ++i)
            {
                if (gi_points)
                {
                    Color_depth_pixel const& c = _data[pixel][i];
                    gi_points->push_back(c.node);
                }
            }

            bool pixel_full = false;

            while (!pixel_full && i < _data[pixel].size())
            {
                Color_depth_pixel const& c = _data[pixel][i];

                if (acc_filling_degree + c.filling_degree > 1.0f)
                {
                    _data_accumulated[pixel] += c.color * (1.0f - acc_filling_degree);
                    total_energy += (c.color * (1.0f - acc_filling_degree)).energy();
                    pixel_full = true;

//                    if (gi_points)
//                    {
//                        gi_points->push_back(c.node);
//                    }


                }
                else
                {
                    acc_filling_degree += c.filling_degree;
                    _data_accumulated[pixel] += c.color * c.filling_degree;
                    // data_accumulated[pixel] += c.color * 1.0f / float(data[pixel].size());
                    // data_accumulated[pixel] += c.color;
                    total_energy += (c.color * c.filling_degree).energy();

//                    if (gi_points)
//                    {
//                        gi_points->push_back(c.node);
//                    }

//                    std::cout << "pixel: " << pixel << " " << data_accumulated[pixel] << std::endl;
                }


                ++i;
            }
        }

        // std::cout << "accumulate(), samples: " << samples_sum << std::endl;

        return total_energy;
    }

    virtual color_t const& get_color(int const x, int const y) const
    {
        int const pixel = (x + _size / 2) + _size * (y + _size / 2);
        return _data_accumulated[pixel];
    }


    /*
    friend std::ostream & operator<<(std::ostream & s, Frame_buffer const& f)
    {
        for (int i = 0; i < size * size; ++i)
        {
            s << f.data_accumulated[i].R << " " << f.data_accumulated[i].G << " " << f.data_accumulated[i].B << " ";
        }

        return s;
    }

    friend std::istream & operator>>(std::istream & s, Frame_buffer & f)
    {
        for (int i = 0; i < size * size; ++i)
        {
            s >> f.data_accumulated[i].R >> f.data_accumulated[i].G >> f.data_accumulated[i].B;
        }

        return s;
    }
*/

protected:
    std::vector< std::vector<Color_depth_pixel> > _data;
    std::vector<color_t> _data_accumulated;
};


class Cube_raster_buffer
{
public:
    enum Type { Simple, Accumulating };

    typedef vector3d_t Point;
    typedef color_t Color;

    Cube_raster_buffer & operator= (Cube_raster_buffer const& cbr)
    {
        if (this != &cbr)
        {
            for (int i = 0; i < 6; ++i)
            {
                buffers[i] = cbr.buffers[i]->clone();
            }

            total_energy = cbr.total_energy;
            _resolution = cbr._resolution;
            _resolution_2 = cbr._resolution_2;
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

    void setup(Type type, int resolution)
    {
        _resolution = resolution;
        _resolution_2 = _resolution / 2;

        for (int i = 0; i < 6; ++i)
        {
            if (type == Simple)
            {
                buffers[i] = new Simple_frame_buffer(_resolution);
            }
            else if (type == Accumulating)
            {
                buffers[i] = new Accumulating_frame_buffer(_resolution);
            }
            else
            {
                std::cout << "Unknown FB type ..." << std::endl;
                assert(false);
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

    Point get_cell_center(Cube_cell const& c) const
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

    // p0, p1: edge to test
    // v0: point on the plane, n: plane normal
    // alpha: return value, intersection point: p0 + alpha * (p1 - p0)
    bool linePlaneIntersection(Point const& p0, Point const& p1, Point const& v0, Point const& n, float & alpha)
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
    bool linePlaneIntersection(Point const& dir, Point const& n, float & alpha)
    // bool linePlaneIntersection(Point const& dir, Point const& v0, Point const& n, float & alpha)
    {
        // float d1 = dot(n, v0);
        float const d1 = 1.0f;
        float const d2 = n * dir;

        if (std::abs(d2) < 1e-10) return false;

        alpha = d1 / d2;

        return (alpha >= 0.0f && alpha <= 1.0f);
    }

    float round(float const d)
    {
        return std::floor(d + 0.5f);
    }

    float next_raster_in_positive_direction(float u)
    {
        return std::floor(u + 1.0f);
    }

    inline int add_point(Point const& direction, Color const& color, float const solid_angle, float const depth, bool use_rays = false)
    {
        if (use_rays)
        {
            return add_point_rays(direction, color, solid_angle, depth);
        }
        else
        {
            add_point_square_rasterization(direction, color, solid_angle, depth);
        }
        return 0;
    }


    // raytrace a disc through every pixel
    void add_point_single_cell(Point const& direction, Color const& color, float const depth, GiPoint const* node = NULL)
    {
        Point normals[6] =
        {
            Point( 1.0f,  0.0f,  0.0f),
            Point(-1.0f,  0.0f,  0.0f),
            Point( 0.0f,  1.0f,  0.0f),
            Point( 0.0f, -1.0f,  0.0f),
            Point( 0.0f,  0.0f,  1.0f),
            Point( 0.0f,  0.0f, -1.0f),
        };

        Point dir = direction;
        dir.normalize();

        int longest_axis = get_longest_axis(dir);

        if (std::abs(dir[longest_axis]) < 1e-10f) return;

        int plane_index = longest_axis * 2;

        if (dir[longest_axis] < 0.0f) ++plane_index;

        assert(dir * normals[plane_index] > 0.0f);

        float alpha = 0.0f;
        linePlaneIntersection(dir, normals[plane_index], alpha);
        if (alpha < 1.0f)
        {
            std::cout << alpha << " " << dir << " " <<  normals[plane_index] << std::endl;
        }
        assert(alpha >= 1.0f);

        Point hit_point = alpha * dir;

        int const axis_0 = (longest_axis + 1) % 3;
        int const axis_1 = (longest_axis + 2) % 3;

        int const u = std::floor(hit_point[axis_0] * _resolution_2);
        int const v = std::floor(hit_point[axis_1] * _resolution_2);

        if (u < -_resolution_2 || u > _resolution_2 - 1) return;
        if (v < -_resolution_2 || v > _resolution_2 - 1) return;

        // std::cout << "hitpoint: " << hit_point << "u: " << u << " v: " << v << std::endl;

        buffers[plane_index]->add_point(u, v, color, 1.0f, depth, node);
    }


    // raytrace a disc through every pixel, assumes receiving point is in the origin and the disc's center relative to it
    void add_point_exact(Color const& color,
                         Point const& disc_normal, Point const& disc_center, float const disc_radius,
                         float const depth,
                         GiPoint const* node = NULL)
    {
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
                            buffers[plane_index]->add_point(u, v, color, 1.0f, alpha, node);
                        }
                    }
                }
            }
        }
    }


    // only check in the positive halfspace of the line, line starts at origin
    bool sphere_line_intersection(Point const& center, float const radius, Point const& line_dir)
    {
        float const c_l = center * line_dir;

        if (c_l <= 0.0f) return false;

        return (c_l * c_l) - center.lengthSqr() + (radius * radius) >= 0.0f;
    }

    // direction of the sphere, its color, radius and distance as seen from the receiving point
    int add_point_from_sphere_with_rays(Point const& direction, Color const& color, float const radius, float const distance, GiPoint const* node = NULL)
    {
        Cube_cell c;

        int hit_cells = 0;

        Point dir = direction;
        dir.normalize();

        Point sphere_center = dir * distance;

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

                    if (sphere_line_intersection(sphere_center, radius, cell_center))
                    {
                        buffers[plane_index]->add_point(u, v, color, 1.0f, distance, node);
                        ++hit_cells;
                    }
                }
            }
        }

        return hit_cells;
    }



    // raytrace against the solid angle
    int add_point_rays(Point const& direction, Color const& color, float const solid_angle, float const depth)
    {
        Cube_cell c;

        int hit_cells = 0;

        Point dir = direction;
        dir.normalize();

        float const disc_angle = std::acos(1.0f - solid_angle / (2.0f * M_PI));
        // float cos_disc_angle = 1.0f - solid_angle / (2.0f * M_PI);

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

                    float const angle_dir_cell_center = std::acos(dir * cell_center);
                    // float cos_dir_cell_center = dir * cell_center;

                    if (angle_dir_cell_center < disc_angle)
                    {
                        buffers[plane_index]->add_point(u, v, color, 1.0f, depth);
                        ++hit_cells;
                    }
                }
            }
        }

        return hit_cells;
    }

    void add_point_square_rasterization(Point const& direction, Color const& color, float const solid_angle, float const depth, GiPoint const* node = NULL)
    {
        /*
        if (std::abs(solid_angle) > 1.0f)
        {
//            std::cout << "solid angle > 1" << std::endl;
            return;
        }
        */


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

/*
        Cube_cell c = get_cell(dir);

        float cell_solid_angle = get_solid_angle(c);

         // if fill_degree > 1.0f, this contribution should spill into neighbouring cells
        float fill_degree = 0.0f;
//        float fill_degree = solid_angle / cell_solid_angle;

        buffers[c.plane].add_point(c.pos[0], c.pos[1], color, fill_degree, depth);
        return;
*/

        // --- new square projection ---

        Point normals[6] =
        {
            Point( 1.0f,  0.0f,  0.0f),
            Point(-1.0f,  0.0f,  0.0f),
            Point( 0.0f,  1.0f,  0.0f),
            Point( 0.0f, -1.0f,  0.0f),
            Point( 0.0f,  0.0f,  1.0f),
            Point( 0.0f,  0.0f, -1.0f),
        };

        // float inv_cell_area = 1.0f / float(resolution_2 * resolution_2);

        // find intersection point on plane

        float acc_area = 0.0f;

        for (int i = 0; i < 3; ++i)
        {
            // int longest_axis = get_longest_axis(dir);
            int longest_axis = i;

            if (std::abs(dir[longest_axis]) < 1e-10f) continue;

            int plane_index = longest_axis * 2;

            if (dir[longest_axis] < 0.0f) ++plane_index;

            assert(dir * normals[plane_index] > 0.0f);

            float alpha = 0.0f;
            linePlaneIntersection(dir, normals[plane_index], alpha);
            if (alpha < 1.0f)
            {
                std::cout << alpha << " " << dir << " " <<  normals[plane_index] << std::endl;
            }
            assert(alpha >= 1.0f);

            Point hit_point = alpha * dir;

            // float const square_area = solid_angle * alpha * alpha / (dir * normals[plane_index]);
            // float const square_area = solid_angle * alpha * alpha * (dir * normals[plane_index]);
            float const square_area = solid_angle; // * (dir * normals[plane_index]);
            float const square_width_2 = std::sqrt(square_area / 2.0f);

            int const axis_0 = (longest_axis + 1) % 3;
            int const axis_1 = (longest_axis + 2) % 3;

            float u_minus = (hit_point[axis_0] - square_width_2) * _resolution_2;
            float u_plus  = (hit_point[axis_0] + square_width_2) * _resolution_2;

            float v_minus = (hit_point[axis_1] - square_width_2) * _resolution_2;
            float v_plus  = (hit_point[axis_1] + square_width_2) * _resolution_2;

            u_minus = into_range(-float(_resolution_2), float(_resolution_2), u_minus);
            v_minus = into_range(-float(_resolution_2), float(_resolution_2), v_minus);

            u_plus = into_range(-float(_resolution_2), float(_resolution_2), u_plus);
            v_plus = into_range(-float(_resolution_2), float(_resolution_2), v_plus);

            float u = u_minus;

            float grid_area = 1.0f;
            // float grid_area = (1.0f / float(resolution)) * (1.0f / float(resolution)) * (resolution_2 * resolution_2);

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

                    acc_area += area;

                    float const ratio = area / grid_area;

                    int const cell_u = into_range(-_resolution_2, _resolution_2 - 1, int(std::floor(u)));
                    int const cell_v = into_range(-_resolution_2, _resolution_2 - 1, int(std::floor(v)));

                    buffers[plane_index]->add_point(cell_u, cell_v, color, ratio, depth, node);
                    // buffers[plane_index].add_point(cell_u, cell_v, color * inv_cell_area, ratio, depth);

                    //                std::cout << "add to cell: " << cell_u << " " << cell_v << " area: " << area << " ratio: " << ratio << std::endl;

                    v = next_v;
                }

                u = next_u;
            }
        }

        if (false && std::abs(solid_angle) > 1e-10f)
        {
            float unit_area = solid_angle;
            float acc_area_ratio = acc_area / unit_area;
            // std::cout << "acc_area_ratio: " << acc_area_ratio << std::endl;
            std::cout << "acc_area: " << acc_area << std::endl;
        }
    }

    float accumulate(std::vector<GiPoint const*> * gi_points = NULL)
    {
        total_energy = 0.0f;

        for (int i = 0; i < 6; ++i)
        {
            total_energy += buffers[i]->accumulate(gi_points);
        }

        return total_energy;

        // std::cout << "accumulate(): total_energy: " << total_energy << std::endl;
    }

    Color const& get_color(Cube_cell const& c) const
    {
        // return buffers[c.plane]._data_accumulated[c.pos[0] + resolution_2 + resolution * (c.pos[1] + resolution_2)];
        return buffers[c.plane]->get_color(c.pos[0], c.pos[1]);
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
        Point const normals[6] =
        {
            Point( 1.0f,  0.0f,  0.0f),
            Point(-1.0f,  0.0f,  0.0f),
            Point( 0.0f,  1.0f,  0.0f),
            Point( 0.0f, -1.0f,  0.0f),
            Point( 0.0f,  0.0f,  1.0f),
            Point( 0.0f,  0.0f, -1.0f),
        };

        // TODO: add solid angle calculation for different cells depending on distance and angle
        float radius = (1.0f / float(_resolution));
        Point p = get_cell_center(c);
        p.normalize();

        // float corrected_area = 2.0f * M_PI * radius * radius * M_PI * 0.5f;

        float solid_angle = (normals[c.plane] * p) * radius * radius;

        return solid_angle;
    }

    Color get_diffuse(Point const& normal) const
    {
        Point const normals[6] =
        {
            Point( 1.0f,  0.0f,  0.0f),
            Point(-1.0f,  0.0f,  0.0f),
            Point( 0.0f,  1.0f,  0.0f),
            Point( 0.0f, -1.0f,  0.0f),
            Point( 0.0f,  0.0f,  1.0f),
            Point( 0.0f,  0.0f, -1.0f),
        };

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

                    Color const& color = get_color(c);

                    float cos_sp_normal_cell_dir = std::max(0.0f, p * normal); // lambertian
                    float cos_plane_normal_cell_dir = std::max(0.0f, normals[c.plane] * p); // solid angle of pixel

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

private:
    Abstract_frame_buffer* buffers[6];
    //Accumulating_frame_buffer<resolution, Color> buffers[6];
    float total_energy;
    int _resolution;
    int _resolution_2;
};

__END_YAFRAY

#endif // CUBERASTERBUFFER_H

