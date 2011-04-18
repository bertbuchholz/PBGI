#ifndef CUBERASTERBUFFER_H
#define CUBERASTERBUFFER_H

#include <algorithm>
#include <vector>

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

template <class Color>
struct Color_depth_pixel
{
    Color color;
    float depth;
    float filling_degree;

    friend bool operator<(Color_depth_pixel const& lhs, Color_depth_pixel const& rhs)
    {
        return (lhs.depth < rhs.depth);
    }

};

template <int size, class Color>
struct Frame_buffer
{
    void add_point(int const x, int const y, Color const& color, float const filling_degree, float const depth)
    {
        int pos = (x + size / 2) + size * (y + size / 2);

        assert(pos < size * size);

        Color_depth_pixel<Color> c;
        c.color = color;
        c.depth = depth;
        c.filling_degree = filling_degree;

        data[pos].push_back(c);
    }

    float accumulate()
    {
//        std::cout << "accumulate()" << std::endl;

        int samples_sum = 0;
        float total_energy = 0.0f;

        for (int pixel = 0; pixel < size * size; ++pixel)
        {
            std::sort(data[pixel].begin(), data[pixel].end());

            samples_sum += data[pixel].size();

            float acc_filling_degree = 0.0f;

            unsigned int i = 0;

            data_accumulated[pixel] = Color(0.0f);

            while (i < data[pixel].size())
            {
                Color_depth_pixel<Color> const& c = data[pixel][i];

                if (acc_filling_degree + c.filling_degree > 1.0f)
                {
                    data_accumulated[pixel] += c.color * (1.0f - acc_filling_degree);
                    total_energy += (c.color * (1.0f - acc_filling_degree)).col2bri();
                    break;
                }
                else
                {
                    acc_filling_degree += c.filling_degree;
                    data_accumulated[pixel] += c.color * c.filling_degree;
                    // data_accumulated[pixel] += c.color;
                    total_energy += c.color.col2bri();

//                    std::cout << "pixel: " << pixel << " " << data_accumulated[pixel] << std::endl;
                }

                ++i;
            }
        }

        // std::cout << "accumulate(), samples: " << samples_sum << std::endl;

        return total_energy;
    }

    friend std::ostream & operator<<(std::ostream & s, Frame_buffer const& f)
    {
        /*
        for (int i = 0; i < size * size; ++i)
        {
            s << f.data[i] << " ";
        }
        */

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

    std::vector< Color_depth_pixel<Color> > data[size * size];

    Color data_accumulated[size * size];
};

template <int resolution, class Point_t, class Color_t>
class Cube_raster_buffer
{
public:
    typedef Point_t Point;
    typedef Color_t Color;
    static const int size = resolution;

    static const int resolution_2 = resolution / 2;

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

        result.pos[0] = std::floor(axis_0_pos * resolution_2);
        result.pos[1] = std::floor(axis_1_pos * resolution_2);

        result.pos[0] = std::max(-resolution_2, std::min(result.pos[0], resolution_2 - 1));
        result.pos[1] = std::max(-resolution_2, std::min(result.pos[1], resolution_2 - 1));

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

        result[axis_0] = (c.pos[0] + 0.5f) / float(resolution_2);
        result[axis_1] = (c.pos[1] + 0.5f) / float(resolution_2);

        return result;
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

    void add_point(Point const& direction, Color const& color, float const solid_angle, float const depth)
    {
        if (std::abs(solid_angle) < 1e-10f)
        {
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

        float inv_cell_area = 1.0f / float(resolution_2 * resolution_2);

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
            float const square_width_2 = std::sqrt(square_area / 4.0f);

            int const axis_0 = (longest_axis + 1) % 3;
            int const axis_1 = (longest_axis + 2) % 3;

            float u_minus = (hit_point[axis_0] - square_width_2) * resolution_2;
            float u_plus  = (hit_point[axis_0] + square_width_2) * resolution_2;

            float v_minus = (hit_point[axis_1] - square_width_2) * resolution_2;
            float v_plus  = (hit_point[axis_1] + square_width_2) * resolution_2;

            u_minus = into_range(-float(resolution_2), float(resolution_2), u_minus);
            v_minus = into_range(-float(resolution_2), float(resolution_2), v_minus);

            u_plus = into_range(-float(resolution_2), float(resolution_2), u_plus);
            v_plus = into_range(-float(resolution_2), float(resolution_2), v_plus);

/*
            u_minus = std::max(-float(resolution_2), u_minus);
            v_minus = std::max(-float(resolution_2), v_minus);

            u_plus = std::min(float(resolution_2), u_plus);
            v_plus = std::min(float(resolution_2), v_plus);
*/

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

                    int const cell_u = into_range(-resolution_2, resolution_2 - 1, int(std::floor(u)));
                    int const cell_v = into_range(-resolution_2, resolution_2 - 1, int(std::floor(v)));

                    buffers[plane_index].add_point(cell_u, cell_v, color * inv_cell_area, ratio, depth);

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
            std::cout << "acc_area_ratio: " << acc_area_ratio << std::endl;
        }
    }

    void accumulate()
    {
        total_energy = 0.0f;

        for (int i = 0; i < 6; ++i)
        {
            total_energy += buffers[i].accumulate();
        }

        // std::cout << "accumulate(): total_energy: " << total_energy << std::endl;
    }

    Color const& get_color(Cube_cell const& c) const
    {
        return buffers[c.plane].data_accumulated[c.pos[0] + resolution_2 + resolution * (c.pos[1] + resolution_2)];
    }

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
        float radius = (1.0f / float(size));
        Point p = get_cell_center(c);
        p.normalize();

        float corrected_area = 2.0f * M_PI * radius * radius * M_PI * 0.5f;

        float solid_angle = (normals[c.plane] * p) * corrected_area;

        return solid_angle;
    }

    Color get_diffuse(Point const& normal)
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

        for (int i = 0; i < 6; ++i)
        {
            c.plane = i;

            for (int x = -size / 2; x < size / 2; ++x)
            {
                c.pos[0] = x;

                for (int y = -size / 2; y < size / 2; ++y)
                {
                    c.pos[1] = y;

                    Point p = get_cell_center(c);
                    Color const& color = get_color(c);

                    float lSqr = p.lengthSqr();

                    p.normalize();

                    float cos_sp_normal_cell_dir = std::max(0.0f, p * normal);
                    float cos_plane_normal_cell_dir = std::max(0.0f, normals[c.plane] * p);

                    diffuse += color * cos_sp_normal_cell_dir * cos_plane_normal_cell_dir; // / lSqr;
                }
            }
        }

        return diffuse; // * inv_cell_area;
    }

    float get_total_energy()
    {
        return total_energy;
    }

private:
    Frame_buffer<resolution, Color> buffers[6];
    float total_energy;
};

#endif // CUBERASTERBUFFER_H
