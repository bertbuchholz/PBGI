#ifndef FRAME_BUFFER_H
#define FRAME_BUFFER_H

#include <cassert>
#include <vector>

#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <utilities/math_utils.h>

#include <integrators/Pbgi_integrator/Gi_point_info.h>

__BEGIN_YAFRAY



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

template <class Data>
class Abstract_frame_buffer
{
public:
    Abstract_frame_buffer() {}

    Abstract_frame_buffer(int const resolution) :
        _resolution(resolution),
        _resolution_2(_resolution / 2)
    { }

    virtual ~Abstract_frame_buffer() {}

    virtual Abstract_frame_buffer* clone() = 0;


    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info) = 0;

    virtual void set_color(int const /* x */, int const /* y */, Data const& /* color */) {}

    virtual Data get_color(int const x, int const y) const = 0;
    virtual Data get_color_interpolated(float const /* x */, float const /* y */) { assert(false); return Data(); }

    virtual void set_size(int size) = 0;

    virtual void clear() = 0;

    inline int calc_serial_index(int const x, int const y) const
    {
        int serial_index = (x + _resolution_2) + _resolution * (y + _resolution_2);
        assert(serial_index < Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
        return serial_index;
    }

    int get_resolution() const
    {
        return _resolution;
    }

protected:
    int _resolution;
    int _resolution_2;
};




// store only a single color value and depth per pixel
template <class Data>
class Simple_single_value_frame_buffer : public Abstract_frame_buffer<Data>
{
public:
    Simple_single_value_frame_buffer() {}

    Simple_single_value_frame_buffer(int size) : Abstract_frame_buffer<Data>(size)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer<Data>* clone()
    {
        Simple_single_value_frame_buffer* sfb = new Simple_single_value_frame_buffer(*this);

        return sfb;
    }

    virtual void add_point(int const /* x */, int const /* y */, Data const& /* color */, float const /* filling_degree */, Gi_point_info const& /* point_info */)
    {
        assert(false);
    }

    virtual void set_size(int size)
    {
        _data.clear();
        _data.resize(Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
    }

    virtual void clear()
    {
        _data.clear();
    }

    void set_color(int const x, int const y, Data const& color)
    {
        int serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        assert(serial_index < Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);

        _data[serial_index] = color;
    }

    virtual Data get_color(int const x, int const y) const
    {
        int serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);
        return _data[serial_index];
    }

protected:
    std::vector<Data> _data;
};


// store a list of color values and corresponding depths and weights (for example filling) per pixel
template <class Data>
class Simple_frame_buffer_without_queue : public Abstract_frame_buffer<Data>
{
    struct Color_filled
    {
        Color_filled() : color(0.0f), filled(false)
        {}

        Data color;
        bool filled;
    };

public:
    virtual ~Simple_frame_buffer_without_queue() {}

    Simple_frame_buffer_without_queue(int size) : Abstract_frame_buffer<Data>(size)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer<Data>* clone()
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

    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info)
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        if (!_data[serial_index].filled)
        {
            _data[serial_index].color = color;
            _data[serial_index].filled = true;
        }
    }

    virtual Data get_color(int const x, int const y) const
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        return _data[serial_index].color;
    }

protected:
    std::vector<Color_filled> _data;
};




// store a list of color values and corresponding depths and weights (for example filling) per pixel
template <class Data>
class Accumulating_frame_buffer_without_queue : public Abstract_frame_buffer<Data>
{
    struct Data_depth_pixel
    {
        Data_depth_pixel() :
            color(color_t(0.0f)),
            filling_degree(0.0f)
        { }

        Data color;
        float filling_degree;
    };

public:
    virtual ~Accumulating_frame_buffer_without_queue() {}

    Accumulating_frame_buffer_without_queue()
    { }

    Accumulating_frame_buffer_without_queue(int size) : Abstract_frame_buffer<Data>(size)
    {
        set_size(size);
    }

    virtual Abstract_frame_buffer<Data>* clone()
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

    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info)
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        Data_depth_pixel & c = _data[serial_index];

        bool accepted = false;
        float weight = filling_degree;

        if (c.filling_degree < 1.0f)
        {
            if (c.filling_degree + weight > 1.0f)
            {
                weight = 1.0f - c.filling_degree;
            }

            accepted = true;
        }

        if (accepted && weight > 0.0001f)
        {
            c.filling_degree += weight;
            c.color += color * weight;
        }
    }

    virtual Data get_color(int const x, int const y) const
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        Data_depth_pixel const& c = _data[serial_index];

        return c.color;
    }


protected:
    std::vector<Data_depth_pixel> _data;
};



__END_YAFRAY

#endif // FRAME_BUFFER_H
