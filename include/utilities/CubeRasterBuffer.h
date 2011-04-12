#ifndef CUBERASTERBUFFER_H
#define CUBERASTERBUFFER_H

#include <algorithm>
#include <vector>

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

        Color_depth_pixel<Color> c;
        c.color = color;
        c.depth = depth;
        c.filling_degree = filling_degree;

        data[pos].push_back(c);
    }

    float accumulate()
    {
        // std::cout << "accumulate()" << std::endl;

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
                    // std::cout << "pixel: " << pixel << " " << data_accumulated[pixel] << std::endl;
                    acc_filling_degree += c.filling_degree;
                    data_accumulated[pixel] += c.color; // * c.filling_degree;
                    total_energy += c.color.col2bri();
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

    Cube_cell get_cell(Point const& dir) const
    {
        Cube_cell result;

        int longest_axis = 0;
        if (std::abs(dir[1]) > std::abs(dir[0]) && std::abs(dir[1]) > std::abs(dir[2])) longest_axis = 1;
        else if (std::abs(dir[2]) > std::abs(dir[0])) longest_axis = 2;

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

        float inv_longest_extent = 1.0f / longest_extent;

        float axis_0_pos = dir[axis_0] * inv_longest_extent;
        float axis_1_pos = dir[axis_1] * inv_longest_extent;

        result.pos[0] = int(axis_0_pos * resolution_2);
        result.pos[1] = int(axis_1_pos * resolution_2);

        if (axis_0_pos < 0) --result.pos[0];
        if (axis_1_pos < 0) --result.pos[1];

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

    void add_point(Point const& dir, Color const& color, float const solid_angle, float const depth)
    {
        Cube_cell c = get_cell(dir);

        float cell_solid_angle = get_solid_angle(c);

         // if fill_degree > 1.0f, this contribution should spill into neighbouring cells
        float fill_degree = 0.0f;
        // float fill_degree = solid_angle / cell_solid_angle;

        buffers[c.plane].add_point(c.pos[0], c.pos[1], color, fill_degree, depth);
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
        // TODO: add solid angle calculation for different cells depending on distance and angle
        float radius = (0.5f / float(size));
        return 2.0f * M_PI * radius * radius;
    }

    Color get_diffuse(Point const& normal)
    {
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

                    p.normalize();

                    float cos_normal_dir = std::max(0.0f, p * normal);

                    diffuse += color * cos_normal_dir;
                }
            }
        }

        return diffuse;
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
