#ifndef FRAME_BUFFER_H
#define FRAME_BUFFER_H

#include <cassert>

#include <core_api/color.h>
#include <utilities/Debug_info.h>


__BEGIN_YAFRAY


struct Color_depth_pixel
{
    Color_depth_pixel() :
        depth(1e10f),
        filling_degree(0.0f),
        radius(0.0f),
        node(NULL)
    {

    }

    color_t color;
    float depth;
    float filling_degree;
    float radius;
    GiPoint const* node;

    friend bool operator<(Color_depth_pixel const& lhs, Color_depth_pixel const& rhs)
    {
        return (lhs.depth < rhs.depth);
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

class Abstract_frame_buffer
{
public:
    Abstract_frame_buffer() {}

    Abstract_frame_buffer(int resolution, int plane) :
        _resolution(resolution),
        _resolution_2(_resolution / 2),
        _debug_plane(plane)
    { }

    virtual ~Abstract_frame_buffer() {}

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, GiPoint const* node = NULL) = 0;

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL) = 0;

    virtual void set_size(int size) = 0;

    virtual void clear() = 0;

    virtual Abstract_frame_buffer* clone() = 0;

protected:
    int _resolution;
    int _resolution_2;
    int _debug_plane;
};

// store only a single color value and depth per pixel
class Simple_frame_buffer : public Abstract_frame_buffer
{
public:
    Simple_frame_buffer() {}

    Simple_frame_buffer(int size, int plane) : Abstract_frame_buffer(size, plane)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Simple_frame_buffer* sfb = new Simple_frame_buffer(*this);

        return sfb;
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, GiPoint const* node = NULL)
    {
        int pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);

        assert(pixel < _resolution * _resolution);

        if (depth < _data[pixel].depth)
        {
            _data[pixel].color = color;
            _data[pixel].depth = depth;
            _data[pixel].filling_degree = filling_degree;
            _data[pixel].node = node;
            _data[pixel].radius = radius;
        }
    }

    virtual void set_size(int size)
    {
        _resolution = size;
        _resolution_2 = _resolution / 2;

        _data.clear();
        _data.resize(_resolution * _resolution);
    }

    virtual void clear()
    {
        _data.clear();
    }


    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        int const pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);

        Color_depth_pixel const& c = _data[pixel];

        if (debug_info && c.node)
        {
            // debug_info->gi_points.push_back(c.node);
            debug_info->gi_points.insert(c.node);

            int cube_x = (pixel % _resolution) - _resolution_2;
            int cube_y = (pixel / _resolution) - _resolution_2;

            if (debug_info->cube_plane == _debug_plane && debug_info->cube_x == cube_x && debug_info->cube_y == cube_y)
            {
                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, 1.0f));
            }
        }

        return c.color;
    }

protected:
    std::vector<Color_depth_pixel> _data;
};

// store a list of color values and corresponding depths and weights (for example filling) per pixel
class Accumulating_frame_buffer : public Abstract_frame_buffer
{
public:
    Accumulating_frame_buffer(int size, int plane) : Abstract_frame_buffer(size, plane)
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
        _resolution = size;
        _resolution_2 = _resolution / 2;

        _data.clear();
        _data.resize(size * size);
    }

    virtual void clear()
    {
        _data.clear();
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, GiPoint const* node = NULL)
    {
        int pos = (x + _resolution / 2) + _resolution * (y + _resolution / 2);

        assert(pos < _resolution * _resolution);

        Color_depth_pixel c;
        c.color = color;
        c.depth = depth;
        c.filling_degree = filling_degree;
        c.radius = radius;
        c.node = node;

        _data[pos].push_back(c);
    }

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        color_t accumulated_color;

        int const pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);

        if (_data[pixel].size() == 0) return accumulated_color;

        bool debug_pixel = false;

        if (debug_info)
        {
            int cube_x = (pixel % _resolution) - _resolution_2;
            int cube_y = (pixel / _resolution) - _resolution_2;

            debug_pixel = debug_info->cube_plane == _debug_plane && debug_info->cube_x == cube_x && debug_info->cube_y == cube_y;
        }

        std::sort(_data[pixel].begin(), _data[pixel].end());

        float acc_filling_degree = 0.0f;


        if (debug_info)
        {
            for (unsigned int i = 0; i < _data[pixel].size(); ++i)
            {
                Color_depth_pixel const& c = _data[pixel][i];
                // debug_info->gi_points.push_back(c.node);
                debug_info->gi_points.insert(c.node);
            }
        }

        unsigned int node_index = 0;
        bool pixel_full = false;

        while (!pixel_full && node_index < _data[pixel].size())
        {
            Color_depth_pixel const& c = _data[pixel][node_index];

            float weight = c.filling_degree;

            if (acc_filling_degree + weight > 1.0f)
            {
                pixel_full = true;
                weight = (1.0f - acc_filling_degree);
            }

            acc_filling_degree += weight;
            accumulated_color += c.color * weight;

            if (debug_pixel)
            {
                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, weight));
            }

            ++node_index;
        }

        if (debug_pixel)
        {
            for (unsigned int i = node_index; i < _data[pixel].size(); ++i)
            {
                Color_depth_pixel const& c = _data[pixel][i];
                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, 0.0f));
            }
        }

        return accumulated_color;
    }


protected:
    std::vector< std::vector<Color_depth_pixel> > _data;
};


class Reweighting_frame_buffer : public Accumulating_frame_buffer
{
public:
    Reweighting_frame_buffer(int size, int plane) : Accumulating_frame_buffer(size, plane)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Reweighting_frame_buffer* afb = new Reweighting_frame_buffer(*this);

        return afb;
    }

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        color_t accumulated_color;

        int const pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);

        if (_data[pixel].size() == 0) return accumulated_color;

        bool debug_pixel = false;

        if (debug_info)
        {
            int cube_x = (pixel % _resolution) - _resolution_2;
            int cube_y = (pixel / _resolution) - _resolution_2;

            debug_pixel = debug_info->cube_plane == _debug_plane && debug_info->cube_x == cube_x && debug_info->cube_y == cube_y;
        }

        std::sort(_data[pixel].begin(), _data[pixel].end());

        float acc_filling_degree = 0.0f;


        if (debug_info)
        {
            for (unsigned int i = 0; i < _data[pixel].size(); ++i)
            {
                Color_depth_pixel const& c = _data[pixel][i];
                // debug_info->gi_points.push_back(c.node);
                debug_info->gi_points.insert(c.node);
            }
        }

        unsigned int node_index = 0;
        bool pixel_full = false;

        while (!pixel_full && node_index < _data[pixel].size())
        {
            Color_depth_pixel const& c = _data[pixel][node_index];

            float weight = c.filling_degree;

            /*
            if (acc_filling_degree + weight > 1.0f)
            {
                pixel_full = true;
                weight = (1.0f - acc_filling_degree);
            }
            */

            acc_filling_degree += weight;
            accumulated_color += c.color * weight;

            if (debug_pixel)
            {
                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, weight));
            }

            ++node_index;
        }

        acc_filling_degree = std::max(1.0f, acc_filling_degree);
        accumulated_color *= 1.0f / float(acc_filling_degree);

        if (debug_pixel)
        {
            for (unsigned int i = node_index; i < _data[pixel].size(); ++i)
            {
                Color_depth_pixel const& c = _data[pixel][i];
                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, 0.0f));
            }
        }

        return accumulated_color;
    }
};


class Distance_weighted_frame_buffer : public Accumulating_frame_buffer
{
public:
    Distance_weighted_frame_buffer(int size, int plane) : Accumulating_frame_buffer(size, plane)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Distance_weighted_frame_buffer* afb = new Distance_weighted_frame_buffer(*this);

        return afb;
    }

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        int const pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);
        std::vector<Color_depth_pixel> & pixel_colors = _data[pixel];

        if (pixel_colors.size() == 0) return color_t(0.0f);

        color_t accumulated_color(0.0f);

        bool debug_pixel = false;

        if (debug_info)
        {
            int cube_x = (pixel % _resolution) - _resolution_2;
            int cube_y = (pixel / _resolution) - _resolution_2;

            debug_pixel = debug_info->cube_plane == _debug_plane && debug_info->cube_x == cube_x && debug_info->cube_y == cube_y;
        }

        std::sort(pixel_colors.begin(), pixel_colors.end());
        std::reverse(pixel_colors.begin(), pixel_colors.end());

        if (debug_info)
        {
            for (unsigned int i = 0; i < pixel_colors.size(); ++i)
            {
                Color_depth_pixel const& c = pixel_colors[i];
                debug_info->gi_points.insert(c.node);
            }
        }



        // float const z_threshold = 0.1f;

        int group_index = 0;

        for (unsigned int i = 0; i < pixel_colors.size(); )
        {
            Color_depth_pixel const& c = pixel_colors[i];

            color_t locally_accumulated_color(0.0f);
            float locally_accumulated_filling_degree(0.0f);
            float weight_sum = 0.0f;

            int const group_start_index = i;

            float previous_depth = c.depth;
            float previous_radius = c.radius;

            int group_size = 0;

            /*
            while (i < pixel_colors.size())
            {
                Color_depth_pixel const& c2 = pixel_colors[i];

                // if (std::abs(c2.depth - c.depth) > c.radius + c2.radius)
                if (std::abs(c2.depth - previous_depth) > previous_radius + c2.radius)
                {
                    break;
                }

                float filling_degree = c2.filling_degree;
                assert(c2.filling_degree <= 1.0f);

                const float weight = filling_degree;
                // locally_accumulated_color += c2.color * weight;
                locally_accumulated_color += c2.color * weight;
                locally_accumulated_filling_degree += filling_degree * weight;
                weight_sum += weight;

                previous_depth = c2.depth;
                previous_radius = c2.radius;

                ++i;
            }

//            weight_sum = std::min(1.0f, weight_sum);

            locally_accumulated_color *= 1.0f / weight_sum;
            locally_accumulated_filling_degree /= weight_sum;
            */

            while (i < pixel_colors.size())
            {
                Color_depth_pixel const& c2 = pixel_colors[i];

                // if (std::abs(c2.depth - c.depth) > c.radius + c2.radius)
                if (std::abs(c2.depth - previous_depth) > previous_radius + c2.radius)
                {
                    break;
                }

                float filling_degree = c2.filling_degree;
                assert(c2.filling_degree <= 1.0f);

                // locally_accumulated_color += c2.color * weight;
                locally_accumulated_color += c2.color * filling_degree;
                locally_accumulated_filling_degree += filling_degree;

                previous_depth = c2.depth;
                previous_radius = c2.radius;

                ++i;
            }

            if (locally_accumulated_filling_degree > 1.0f)
            {
                locally_accumulated_color *= 1.0f / locally_accumulated_filling_degree;
            }

            locally_accumulated_filling_degree = std::min(1.0f, locally_accumulated_filling_degree);

            // accumulated_color = accumulated_color * (1.0f - weight_sum) + locally_accumulated_color * weight_sum;
            accumulated_color = accumulated_color * (1.0f - locally_accumulated_filling_degree) + locally_accumulated_color * locally_accumulated_filling_degree;


            if (debug_pixel)
            {
                for (unsigned int j = 0; j < debug_info->single_pixel_contributors.size(); ++j)
                {
                    debug_info->single_pixel_contributors[j].weight *= (1.0f - weight_sum);
                }

                for (unsigned int j = group_start_index; j < i; ++j)
                {
                    debug_info->single_pixel_contributors.push_back(Node_weight_pair(pixel_colors[j].node, weight_sum, group_index));
                }
            }

            ++group_index;
        }

        /*
        if (debug_pixel)
        {
            for (unsigned int i = node_index; i < _data[pixel].size(); ++i)
            {
                Color_depth_pixel const& c = _data[pixel][i];
                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, 0.0f));
            }
        }
        */

        return accumulated_color;
    }
};


__END_YAFRAY

#endif // FRAME_BUFFER_H
