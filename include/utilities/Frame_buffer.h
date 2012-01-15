#ifndef FRAME_BUFFER_H
#define FRAME_BUFFER_H

#include <cassert>
#include <vector>

#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <utilities/Debug_info.h>


__BEGIN_YAFRAY


struct Color_depth_pixel
{
    Color_depth_pixel() :
        depth(1e10f),
        filling_degree(0.0f),
        radius(0.0f),
        cone_angle(0.0f),
        node(NULL)
    {

    }

    color_t color;
    float depth;
    float filling_degree;
    float radius;
    float cone_angle;
    vector3d_t direction;
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

    virtual Abstract_frame_buffer* clone() = 0;

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL) = 0;

    virtual void    set_color(int const /* x */, int const /* y */, color_t const& /* color */) {}
    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL) = 0;
    virtual color_t get_color_interpolated(float const /* x */, float const /* y */, Debug_info * debug_info = NULL) { assert(false); return color_t(); }

    virtual void set_size(int size) = 0;

    virtual void clear() = 0;

    inline int calc_serial_index(int const x, int const y) const
    {
        return (x + _resolution_2) + _resolution * (y + _resolution_2);
    }

protected:
    int _resolution;
    int _resolution_2;
    int _debug_plane;
};




// store only a single color value and depth per pixel
class Simple_single_value_frame_buffer : public Abstract_frame_buffer
{
public:
    Simple_single_value_frame_buffer() {}

    Simple_single_value_frame_buffer(int size, int plane) : Abstract_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Simple_single_value_frame_buffer(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Simple_single_value_frame_buffer* sfb = new Simple_single_value_frame_buffer(*this);

        return sfb;
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL)
    {
        assert(false);
    }

    virtual void set_size(int size)
    {
        _data.clear();
        _data.resize(_resolution * _resolution);
    }

    virtual void clear()
    {
        _data.clear();
    }

    void set_color(int const x, int const y, color_t const& color)
    {
        int serial_index = calc_serial_index(x, y);

        assert(serial_index < _resolution * _resolution);

        _data[serial_index] = color;
    }

    virtual color_t get_color(int const x, int const y, Debug_info * /* debug_info */ = NULL)
    {
        int serial_index = calc_serial_index(x, y);
        return _data[serial_index];
    }

protected:
    std::vector<color_t> _data;
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

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Simple_frame_buffer(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Simple_frame_buffer* sfb = new Simple_frame_buffer(*this);

        return sfb;
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL)
    {
        int pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);

        assert(pixel < _resolution * _resolution);

        if (depth < _data[pixel].depth)
        {
            _data[pixel].color          = color;
            _data[pixel].depth          = depth;
            _data[pixel].filling_degree = filling_degree;
            _data[pixel].node           = node;
            _data[pixel].radius         = radius;
            _data[pixel].cone_angle     = std::atan2(radius, depth);
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

    void set_color(int const x, int const y, color_t const& color)
    {
        int pixel = (x + _resolution_2) + _resolution * (y + _resolution_2);

        assert(pixel < _resolution * _resolution);

        _data[pixel].color = color;
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

    color_t interpolate(color_t const& i1 , color_t const& i2 , float const offset)
    {
        assert(offset >= 0.0f && offset <= 1.0f);
        return i1 * (1.0f - offset) + i2 * offset;
    }

    virtual color_t get_color_interpolated(float const x, float const y, Debug_info * debug_info = NULL)
    {
        int x0 = int(std::floor(x - 0.5f));
        int y0 = int(std::floor(y - 0.5f));

        int x1 = x0 + 1;
        int y1 = y0 + 1;

        if (x1 >= _resolution_2)
        {
            x1 = x0;
        }
        else if (x0 < -_resolution_2)
        {
            x0 = x1;
        }

        if (y1 >= _resolution_2)
        {
            y1 = y0;
        }
        else if (y0 < -_resolution_2)
        {
            y0 = y1;
        }


        //x0 = std::max(-_resolution_2 + 1, std::min(int(std::floor(x - 0.5f)), _resolution_2 - 2));
        //y0 = std::max(-_resolution_2 + 1, std::min(int(std::floor(y - 0.5f)), _resolution_2 - 2));

        float offset_x = x - x0 - 0.5f;
        float offset_y = y - y0 - 0.5f;

        offset_x = std::max(0.0f, std::min(offset_x, 1.0f));
        offset_y = std::max(0.0f, std::min(offset_y, 1.0f));

        color_t const c1 = get_color(x0, y0);
        color_t const c2 = get_color(x1, y0);
        color_t const c3 = get_color(x0, y1);
        color_t const c4 = get_color(x1, y1);

        color_t const c1_2 = interpolate(c1, c2, offset_x);
        color_t const c3_4 = interpolate(c3, c4, offset_x);

        if (!(x1 < _resolution_2 && y1 < _resolution_2))
        {
            std::cout << "_resolution_2: " << _resolution_2 << " x: " << x << " y: " << y << " x_0: " << x0 << " y_0: " << y0 << std::endl;
            assert(false);
        }

        return interpolate(c1_2, c3_4, offset_y);
    }

protected:
    std::vector<Color_depth_pixel> _data;
};



// store a list of color values and corresponding depths and weights (for example filling) per pixel
class Simple_frame_buffer_without_queue : public Abstract_frame_buffer
{
    struct Color_filled
    {
        color_t color;
        bool filled;
    };

public:
    virtual ~Simple_frame_buffer_without_queue() {}

    Simple_frame_buffer_without_queue(int size, int plane) : Abstract_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Simple_frame_buffer_without_queue(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Simple_frame_buffer_without_queue* afb = new Simple_frame_buffer_without_queue(*this);

        return afb;
    }

    virtual void set_size(int size)
    {
        _data.clear();
        _data.resize(size * size);
    }

    virtual void clear()
    {
        _data.clear();
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL)
    {
        int const serial_index = calc_serial_index(x, y);

        assert(serial_index < _resolution * _resolution);

        if (!_data[serial_index].filled)
        {
            _data[serial_index].color = color;
            _data[serial_index].filled = true;

            if (node)
            {
                if (_single_pixel_contributors.size() == 0)
                {
                    _single_pixel_contributors.resize(_resolution * _resolution);
                }

                Node_weight_pair nwp(node, filling_degree, 0);
                _single_pixel_contributors[serial_index].push_back(nwp);

                _debug_gi_points.insert(node);
            }
        }
    }

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        int const serial_index = calc_serial_index(x, y);

        if (debug_info)
        {
            bool debug_pixel = debug_info->cube_plane == _debug_plane && debug_info->cube_x == x && debug_info->cube_y == y;

            for (std::tr1::unordered_set<yafaray::GiPoint const*>::iterator iter = _debug_gi_points.begin(); iter != _debug_gi_points.end(); ++iter)
            {
                debug_info->gi_points.insert(*iter);
            }

            if (debug_pixel && _single_pixel_contributors[serial_index].size() > 0)
            {
                debug_info->single_pixel_contributors = _single_pixel_contributors[serial_index];
            }
        }

//        std::cout << "Simple_frame_buffer_without_queue::get_color(): " << _data[serial_index].color << std::endl;

        return _data[serial_index].color;
    }

protected:
    std::vector<Color_filled> _data;

    std::tr1::unordered_set<GiPoint const*> _debug_gi_points;
    std::vector< std::vector<Node_weight_pair> > _single_pixel_contributors;
};


// store a list of color values and corresponding depths and weights (for example filling) per pixel
class Accumulating_frame_buffer_without_queue : public Abstract_frame_buffer
{
public:
    virtual ~Accumulating_frame_buffer_without_queue() {}

    Accumulating_frame_buffer_without_queue(int size, int plane) : Abstract_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Accumulating_frame_buffer_without_queue(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Accumulating_frame_buffer_without_queue* afb = new Accumulating_frame_buffer_without_queue(*this);

        return afb;
    }

    virtual void set_size(int size)
    {
        _data.clear();
        _data.resize(size * size);
    }

    virtual void clear()
    {
        _data.clear();
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL)
    {
        int const serial_index = calc_serial_index(x, y);

        assert(serial_index < _resolution * _resolution);

        Color_depth_pixel & c = _data[serial_index];

        bool use_visibility = false;

        if (use_visibility)
        {
            if (c.filling_degree < 1.0f)
            {
                float weight = filling_degree;

                if (c.filling_degree + weight > 1.0f)
                {
                    weight = (1.0f - c.filling_degree);
                }

                c.filling_degree += weight;
                c.color += color * weight;
            }
        }
        else
        {
            c.filling_degree += filling_degree;
            c.color += color * filling_degree;
        }
    }

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        int const serial_index = calc_serial_index(x, y);

        Color_depth_pixel const& c = _data[serial_index];

        if (c.filling_degree > 1.0f)
        {
            return c.color / c.filling_degree;
        }

        return c.color;
    }

protected:
    std::vector<Color_depth_pixel> _data;
};



// store a list of color values and corresponding depths and weights (for example filling) per pixel
class Grouping_frame_buffer_without_queue : public Abstract_frame_buffer
{
    struct Color_weight_group
    {
        color_t color;
        float weight;
    };

public:
    virtual ~Grouping_frame_buffer_without_queue() {}

    Grouping_frame_buffer_without_queue(int size, int plane) : Abstract_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Grouping_frame_buffer_without_queue(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Grouping_frame_buffer_without_queue* afb = new Grouping_frame_buffer_without_queue(*this);

        return afb;
    }

    virtual void set_size(int size)
    {
        _groups.clear();
        _groups.resize(size * size);

        _last_added_pixel_depths.resize(size * size, 1e10);
        _last_added_pixel_radius.resize(size * size, 0);

        _single_pixel_contributors.resize(size * size);
    }

    virtual void clear()
    {
        //_data.clear();
        _groups.clear();
    }

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL)
    {
        int const serial_index = calc_serial_index(x, y);

        assert(serial_index < _resolution * _resolution);

        if (std::abs(depth - _last_added_pixel_depths[serial_index]) > _last_added_pixel_radius[serial_index] * 2.0f)
        {
            // new group
            Color_weight_group cd;
            cd.color  = color * filling_degree;
            cd.weight = filling_degree;
            _groups[serial_index].push_back(cd);
        }
        else
        {
            assert(_groups[serial_index].size() > 0);
            Color_weight_group & cd = _groups[serial_index].back();
            cd.color  += color * filling_degree;
            cd.weight += filling_degree;
        }

        if (node)
        {
            Node_weight_pair nwp(node, filling_degree, _groups[serial_index].size() - 1);
            _single_pixel_contributors[serial_index].push_back(nwp);

            // std::cout << serial_index << " " << depth << " " << _last_added_pixel_depths[serial_index] << " " << radius << std::endl;
            _debug_gi_points.insert(node);
        }

        _last_added_pixel_depths[serial_index] = depth;
        _last_added_pixel_radius[serial_index] = radius;
    }

    virtual color_t get_color(int const x, int const y, Debug_info * debug_info = NULL)
    {
        int const serial_index = calc_serial_index(x, y);
        assert(serial_index < _resolution * _resolution && serial_index >= 0);

        std::vector<Color_weight_group> const& pixel_groups = _groups[serial_index];

        color_t accumulated_color(0.0f);

        for (int i = pixel_groups.size() - 1; i >= 0; --i)
        {
            float weight = pixel_groups[i].weight;

            color_t color = pixel_groups[i].color;

            if (weight > 1.0f)
            {
                color *= 1.0f / weight;
                weight = 1.0f;
            }

            accumulated_color = accumulated_color * (1.0f - weight) + color;
        }

        if (debug_info)
        {
            for (std::tr1::unordered_set<yafaray::GiPoint const*>::iterator iter = _debug_gi_points.begin(); iter != _debug_gi_points.end(); ++iter)
            {
                debug_info->gi_points.insert(*iter);
            }

            bool const debug_pixel = debug_info->cube_plane == _debug_plane && debug_info->cube_x == x && debug_info->cube_y == y;

            if (debug_pixel)
            {
                for (std::size_t i = 0; i < _single_pixel_contributors[serial_index].size(); ++i)
                {
                    int const group_index =_single_pixel_contributors[serial_index][i].group_index;
                    _single_pixel_contributors[serial_index][i].group_weight = pixel_groups[group_index].weight;
                }

                debug_info->single_pixel_contributors = _single_pixel_contributors[serial_index];
            }


            // std::cout << "debug_gi_points: " << debug_gi_points.size() << std::endl;
        }

        return accumulated_color;
    }

protected:
    std::vector<float> _last_added_pixel_depths;
    std::vector<float> _last_added_pixel_radius;
    std::vector< std::vector<Color_weight_group> > _groups;

    std::tr1::unordered_set<GiPoint const*> _debug_gi_points;
    std::vector< std::vector<Node_weight_pair> > _single_pixel_contributors;
};



// store a list of color values and corresponding depths and weights (for example filling) per pixel
class Accumulating_frame_buffer : public Abstract_frame_buffer
{
public:
    virtual ~Accumulating_frame_buffer() {}

    Accumulating_frame_buffer(int size, int plane) : Abstract_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Accumulating_frame_buffer(resolution, plane);
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

    virtual void add_point(int const x, int const y, color_t const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
                           GiPoint const* node = NULL)
    {
        int pixel = (x + _resolution / 2) + _resolution * (y + _resolution / 2);

        assert(pixel < _resolution * _resolution);

        Color_depth_pixel c;
        c.color          = color;
        c.depth          = depth;
        c.filling_degree = filling_degree;
        c.radius         = radius;
        c.direction      = direction;
        c.cone_angle     = std::atan2(radius, depth);
        c.node           = node;

        _data[pixel].push_back(c);
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
        bool const visibility_disabled = true;

        while (!pixel_full && node_index < _data[pixel].size())
        {
            Color_depth_pixel const& c = _data[pixel][node_index];

            float weight = c.filling_degree;

            if (!visibility_disabled && acc_filling_degree + weight > 1.0f)
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

        if (acc_filling_degree > 1.0f)
        {
            accumulated_color *= 1.0f / acc_filling_degree;
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

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Reweighting_frame_buffer(resolution, plane);
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

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Distance_weighted_frame_buffer(resolution, plane);
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
            // float weight_sum = 0.0f;

            int const group_start_index = i;

            float previous_depth = c.depth;
            float previous_radius = c.radius;

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

                // float const stackedness = std::abs(c2.depth - previous_depth) / (5.0f * (previous_radius + c2.radius) * 0.5f);

                // if (std::abs(c2.depth - c.depth) > (c.radius + c2.radius) * 4.0f)
                if (std::abs(c2.depth - c.depth) > (c.radius) * 2.0f)
                // if (std::abs(c2.depth - c.depth) > (c.radius + c2.radius))
                // if (std::abs(c2.depth - previous_depth) > previous_radius + c2.radius)
                // if (std::abs(c2.depth - previous_depth) > (previous_radius + c2.radius) * 1.0f)
                {
                    break;
                }

                // float stacked_weight = std::max(0.0f, 1.0f - stackedness);
                float const stacked_weight = 1.0f;

                float const filling_degree = c2.filling_degree;
                assert(filling_degree <= 1.0f);

                // locally_accumulated_color += c2.color * weight;
                locally_accumulated_color += c2.color * filling_degree * stacked_weight;
                locally_accumulated_filling_degree += filling_degree * stacked_weight;

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
            // accumulated_color = accumulated_color * (1.0f - locally_accumulated_filling_degree) + locally_accumulated_color * locally_accumulated_filling_degree;
            accumulated_color = accumulated_color * (1.0f - locally_accumulated_filling_degree) + locally_accumulated_color;

            if (debug_pixel)
            {
                for (unsigned int j = 0; j < debug_info->single_pixel_contributors.size(); ++j)
                {
                    debug_info->single_pixel_contributors[j].group_weight *= (1.0f - locally_accumulated_filling_degree);
                    // debug_info->single_pixel_contributors[j].group_weight *= (1.0f - weight_sum);
                }

                std::cout << "Group: " << group_index << " fill: " << locally_accumulated_filling_degree << std::endl;

                for (unsigned int j = group_start_index; j < i; ++j)
                {
                    debug_info->single_pixel_contributors.push_back(
                                Node_weight_pair(pixel_colors[j].node, pixel_colors[j].filling_degree, group_index, locally_accumulated_filling_degree));
                    // debug_info->single_pixel_contributors.push_back(Node_weight_pair(pixel_colors[j].node, pixel_colors[j].filling_degree, group_index, weight_sum));


                    std::cout <<
                                 "c.filling_degree: " << pixel_colors[j].filling_degree << " " <<
                                 "c.depth: " << pixel_colors[j].depth << " " <<
                                 "c.radius: " << pixel_colors[j].radius << " " <<
                                 "c.color: " << pixel_colors[j].color << " " <<

                                 std::endl;
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

        if (debug_pixel)
        {
            std::cout << "Pixel acc color: " << accumulated_color << std::endl;
        }

        return accumulated_color;
    }
};




class Forward_distance_weighted_frame_buffer : public Accumulating_frame_buffer
{
public:
    Forward_distance_weighted_frame_buffer(int size, int plane) : Accumulating_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Forward_distance_weighted_frame_buffer(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Forward_distance_weighted_frame_buffer* afb = new Forward_distance_weighted_frame_buffer(*this);

        return afb;
    }


    struct Color_distance_group
    {
        color_t color;
        float weight;
    };


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

        if (debug_info)
        {
            for (unsigned int i = 0; i < pixel_colors.size(); ++i)
            {
                Color_depth_pixel const& c = pixel_colors[i];
                debug_info->gi_points.insert(c.node);
            }
        }


        std::vector<Color_distance_group> color_distance_groups;

        int group_index = 0;

        bool finished = false;

        for (unsigned int i = 0; i < pixel_colors.size() && !finished; )
        {
            Color_depth_pixel const& c = pixel_colors[i];

            color_t locally_accumulated_color(0.0f);
            float locally_accumulated_filling_degree(0.0f);
            // float weight_sum = 0.0f;

            int const group_start_index = i;

            while (i < pixel_colors.size())
            {
                Color_depth_pixel const& c2 = pixel_colors[i];

                if (std::abs(c.depth - c2.depth) > (c.radius) * 2.0f)
                {
                    break;
                }

                float const filling_degree = c2.filling_degree;
                assert(filling_degree <= 1.0f);

                // locally_accumulated_color += c2.color * weight;
                locally_accumulated_color += c2.color * filling_degree;
                locally_accumulated_filling_degree += filling_degree;

                ++i;
            }

            if (locally_accumulated_filling_degree >= 1.0f)
            {
                locally_accumulated_color *= 1.0f / locally_accumulated_filling_degree;
                locally_accumulated_filling_degree = std::min(1.0f, locally_accumulated_filling_degree);
                finished = true;
            }

            Color_distance_group cd;
            cd.color = locally_accumulated_color;
            cd.weight = locally_accumulated_filling_degree;
            color_distance_groups.push_back(cd);

            ++group_index;

            if (debug_pixel)
            {
                for (unsigned int j = 0; j < debug_info->single_pixel_contributors.size(); ++j)
                {
                    debug_info->single_pixel_contributors[j].group_weight *= (1.0f - locally_accumulated_filling_degree);
                    // debug_info->single_pixel_contributors[j].group_weight *= (1.0f - weight_sum);
                }

                std::cout << "Group: " << group_index << " fill: " << locally_accumulated_filling_degree << std::endl;

                for (unsigned int j = group_start_index; j < i; ++j)
                {
                    debug_info->single_pixel_contributors.push_back(
                                Node_weight_pair(pixel_colors[j].node, pixel_colors[j].filling_degree, group_index, locally_accumulated_filling_degree));
                    // debug_info->single_pixel_contributors.push_back(Node_weight_pair(pixel_colors[j].node, pixel_colors[j].filling_degree, group_index, weight_sum));


                    std::cout <<
                                 "c.filling_degree: " << pixel_colors[j].filling_degree << " " <<
                                 "c.depth: " << pixel_colors[j].depth << " " <<
                                 "c.radius: " << pixel_colors[j].radius << " " <<
                                 "c.color: " << pixel_colors[j].color << " " <<

                                 std::endl;
                }
            }

        }

        std::reverse(color_distance_groups.begin(), color_distance_groups.end());

        for (unsigned int i = 0; i < color_distance_groups.size(); ++i)
        {
            accumulated_color = accumulated_color * (1.0f - color_distance_groups[i].weight) + color_distance_groups[i].color;
        }

        if (debug_pixel)
        {
            std::cout << "Pixel acc color: " << accumulated_color << std::endl;
        }

        return accumulated_color;
    }
};




class Distance_weighted_frame_buffer_XXX : public Accumulating_frame_buffer
{
public:
    Distance_weighted_frame_buffer_XXX(int size, int plane) : Accumulating_frame_buffer(size, plane)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Distance_weighted_frame_buffer_XXX* afb = new Distance_weighted_frame_buffer_XXX(*this);

        return afb;
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Distance_weighted_frame_buffer_XXX(resolution, plane);
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
            // float weight_sum = 0.0f;

            int const group_start_index = i;

            float previous_depth = c.depth;
            float previous_radius = c.radius;

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

                // if (std::abs(c2.depth - c.depth) > (c.radius + c2.radius) * 4.0f)
                // if (std::abs(c2.depth - c.depth) > (c.radius + c2.radius))
                // if (std::abs(c2.depth - previous_depth) > previous_radius + c2.radius)
                if (std::abs(c2.depth - previous_depth) > (previous_radius + c2.radius) * 1.0f)
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
                    debug_info->single_pixel_contributors[j].group_weight *= (1.0f - locally_accumulated_filling_degree);
                    // debug_info->single_pixel_contributors[j].group_weight *= (1.0f - weight_sum);
                }

                for (unsigned int j = group_start_index; j < i; ++j)
                {
                    debug_info->single_pixel_contributors.push_back(
                                Node_weight_pair(pixel_colors[j].node, pixel_colors[j].filling_degree, group_index, locally_accumulated_filling_degree));
                    // debug_info->single_pixel_contributors.push_back(Node_weight_pair(pixel_colors[j].node, pixel_colors[j].filling_degree, group_index, weight_sum));
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



class Directional_stacked_frame_buffer : public Accumulating_frame_buffer
{
public:
    Directional_stacked_frame_buffer(int size, int plane) : Accumulating_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Directional_stacked_frame_buffer(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Directional_stacked_frame_buffer* afb = new Directional_stacked_frame_buffer(*this);

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

        if (debug_info)
        {
            for (unsigned int i = 0; i < pixel_colors.size(); ++i)
            {
                Color_depth_pixel const& c = pixel_colors[i];
                debug_info->gi_points.insert(c.node);
            }
        }

        float accumulated_filling_degree = 0.0f;

        for (unsigned int i = 0; i < pixel_colors.size() && accumulated_filling_degree < 1.0f; ++i)
        {
            Color_depth_pixel const& c = pixel_colors[i];

            float stackedness = 0.0f;

            for (unsigned int j = 0; j < i && stackedness < 1.0f; ++j)
            {
                Color_depth_pixel const& c_in_front = pixel_colors[j];

                assert(c_in_front.depth <= c.depth);

                float const angle_between_cones = std::acos(c.direction * c_in_front.direction);
                float const angle_overlap       = std::max(0.0f, c.cone_angle + c_in_front.cone_angle - angle_between_cones);

                float const stackedness_tmp = angle_overlap / (c.cone_angle + c_in_front.cone_angle);

                if (stackedness_tmp > stackedness)
                {
                    stackedness = stackedness_tmp;
                }

                if (debug_pixel)
                {
                    std::cout << "Stacking, i: " << i << " " <<
                                 "j: " << j << " " <<
                                 "angle_between_cones: " << angle_between_cones << " " <<
                                 "angle_overlap: " << angle_overlap << " " <<
                                 "stackedness_tmp: " << stackedness_tmp << " " <<
                                 std::endl;
                }
            }

            float filling_degree = std::max(0.0f, 1.0f - stackedness) * c.filling_degree;

            if (accumulated_filling_degree + filling_degree > 1.0f)
            {
                filling_degree = 1.0f - accumulated_filling_degree + 0.01f;
            }

            accumulated_color += c.color * filling_degree;
            accumulated_filling_degree += filling_degree;

            if (debug_pixel)
            {
                std::cout << "Directional_stacked_frame_buffer: " <<
                             "stackedness: " << stackedness << " " <<
                             "filling_degree: " << filling_degree << " " <<
                             "c.filling_degree: " << c.filling_degree << " " <<
                             "c.depth: " << c.depth << " " <<
                             "c.radius: " << c.radius << " " <<
                             std::endl;

                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, filling_degree));
            }
        }

        return accumulated_color;
    }
};




class Front_stacked_frame_buffer : public Accumulating_frame_buffer
{
public:
    Front_stacked_frame_buffer(int size, int plane) : Accumulating_frame_buffer(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer* create(int resolution, int plane)
    {
        return new Front_stacked_frame_buffer(resolution, plane);
    }

    virtual Abstract_frame_buffer* clone()
    {
        Front_stacked_frame_buffer* afb = new Front_stacked_frame_buffer(*this);

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

        if (debug_info)
        {
            for (unsigned int i = 0; i < pixel_colors.size(); ++i)
            {
                Color_depth_pixel const& c = pixel_colors[i];
                debug_info->gi_points.insert(c.node);
            }
        }

        float accumulated_filling_degree = 0.0f;

        for (unsigned int i = 0; i < pixel_colors.size() && accumulated_filling_degree < 1.0f; ++i)
        {
            Color_depth_pixel const& c = pixel_colors[i];

            float stackedness = 0.0f;

            for (unsigned int j = 0; j < i; ++j)
            {
                Color_depth_pixel const& c2 = pixel_colors[j];

                assert(c2.depth <= c.depth);

                float const stackedness_tmp = std::abs(c2.depth - c.depth) / (2.0f * (c.radius + c2.radius) * 0.5f) * c2.filling_degree / c.filling_degree;

                if (stackedness_tmp > stackedness)
                {
                    stackedness = stackedness_tmp;
                }
            }

            float filling_degree = std::max(0.0f, 1.0f - stackedness) * c.filling_degree;

            if (accumulated_filling_degree  + filling_degree > 1.0f)
            {
                filling_degree = 1.0f - accumulated_filling_degree + 0.01f;
            }

            accumulated_color += c.color * filling_degree;
            accumulated_filling_degree += filling_degree;

            if (debug_pixel)
            {
                std::cout << "Front_stacked_frame_buffer: " <<
                             "stackedness: " << stackedness << " " <<
                             "filling_degree: " << filling_degree << " " <<
                             "c.filling_degree: " << c.filling_degree << " " <<
                             "c.depth: " << c.depth << " " <<
                             "c.radius: " << c.radius << " " <<
                             std::endl;

                debug_info->single_pixel_contributors.push_back(Node_weight_pair(c.node, filling_degree));
            }
        }

        return accumulated_color;
    }
};




static void test_directional_fb()
{
    yafaray::Directional_stacked_frame_buffer fb(8, 0);

    fb.add_point(0, 0, yafaray::color_t(1.0f), 0.5f, 1.0f, 0.5f, yafaray::vector3d_t(1.0f, 0.0f, 0.0f));
    fb.add_point(0, 0, yafaray::color_t(1.0f), 0.5f, 1.0f, 0.5f, yafaray::vector3d_t(1.0f, 1.0f, 0.0f).normalize());

    fb.get_color(0, 0);
}


__END_YAFRAY

#endif // FRAME_BUFFER_H
