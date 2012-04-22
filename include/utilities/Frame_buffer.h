#ifndef FRAME_BUFFER_H
#define FRAME_BUFFER_H

#include <cassert>
#include <vector>

#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <utilities/Debug_info.h>
#include <utilities/math_utils.h>
#include <utilities/spherical_harmonics.h>
#include <utilities/Mises_fisher.h>

#include <integrators/Gi_point_info.h>
#include <bert/shared/Registry_parameters.h>
#include <bert/shared/Parameter.h>

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

    Abstract_frame_buffer(int resolution, int plane) :
        _resolution(resolution),
        _resolution_2(_resolution / 2),
        _debug_plane(plane)
    { }

    virtual ~Abstract_frame_buffer() {}

    virtual Abstract_frame_buffer* clone() = 0;

    virtual void set_parameters(Parameter_list const& parameters)
    {
        _resolution = parameters["resolution"]->get_value<int>();
        _resolution += _resolution % 2;
        _resolution_2 = _resolution / 2;
    }

//    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, float const depth, float const radius, vector3d_t const& direction,
//                           Spherical_node_representation const* sf_representation,
//                           Gi_point_base const* node = NULL) = 0;

    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info) = 0;

    virtual void set_color(int const /* x */, int const /* y */, Data const& /* color */) {}

    virtual Data get_color(int const x, int const y, Debug_info * debug_info = NULL) const = 0;
    virtual Data get_color_interpolated(float const /* x */, float const /* y */, Debug_info * debug_info = NULL) { assert(false); return Data(); }

    virtual void set_size(int size) = 0;

    virtual void clear() = 0;

    inline int calc_serial_index(int const x, int const y) const
    {
        int serial_index = (x + _resolution_2) + _resolution * (y + _resolution_2);
        assert(serial_index < Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
        return serial_index;
    }

    void set_plane(int const plane)
    {
        _debug_plane = plane;
    }

    int get_resolution() const
    {
        return _resolution;
    }

    static Parameter_list get_parameters()
    {
        Parameter_list parameters;

        parameters.add_parameter(new Parameter("resolution", 8, 2, 128));

        return parameters;
    }

protected:
    int _resolution;
    int _resolution_2;
    int _debug_plane;
};




// store only a single color value and depth per pixel
template <class Data>
class Simple_single_value_frame_buffer : public Abstract_frame_buffer<Data>
{
public:
    Simple_single_value_frame_buffer() {}

    Simple_single_value_frame_buffer(int size, int plane) : Abstract_frame_buffer<Data>(size, plane)
    {
        set_size(size);
    }

    void set_parameters(Parameter_list const& parameters) {}

//    static Abstract_frame_buffer<Data>* create(int resolution, int plane)
//    {
//        return new Simple_single_value_frame_buffer(resolution, plane);
//    }

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

    virtual Data get_color(int const x, int const y, Debug_info * /* debug_info */ = NULL) const
    {
        int serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);
        return _data[serial_index];
    }

    static Simple_single_value_frame_buffer * create()
    {
        return new Simple_single_value_frame_buffer(1, 1);
    }

protected:
    std::vector<Data> _data;
};


static Base_registration< Abstract_frame_buffer<color_t> > Simple_single_value_frame_buffer_color_t("Accumulating_frame_buffer_without_queue",
                                                                                                         &Simple_single_value_frame_buffer<color_t>::create);



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

    Simple_frame_buffer_without_queue(int size, int plane) : Abstract_frame_buffer<Data>(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer<Data>* create(int resolution, int plane)
    {
        return new Simple_frame_buffer_without_queue(resolution, plane);
    }

    virtual Abstract_frame_buffer<Data>* clone()
    {
        Simple_frame_buffer_without_queue* afb = new Simple_frame_buffer_without_queue(*this);

        return afb;
    }

    void set_parameters(Parameter_list const& parameters) {}

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

#ifdef DEBUG
            if (_single_pixel_contributors.size() == 0)
            {
                _single_pixel_contributors.resize(Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
            }

            Node_weight_pair nwp(point_info, filling_degree, color);
            _single_pixel_contributors[serial_index].push_back(nwp);

            point_info.splatted = true;
#endif
        }

    }

    virtual Data get_color(int const x, int const y, Debug_info * debug_info = NULL) const
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

#ifdef DEBUG
        bool debug_pixel = debug_info->cube_plane == Abstract_frame_buffer<Data>::_debug_plane && debug_info->cube_x == x && debug_info->cube_y == y;

        if (debug_pixel && _single_pixel_contributors[serial_index].size() > 0)
        {
            debug_info->single_pixel_contributors = _single_pixel_contributors[serial_index];
        }
#endif


//        std::cout << "Simple_frame_buffer_without_queue::get_color(): " << _data[serial_index].color << std::endl;

        return _data[serial_index].color;
    }

protected:
    std::vector<Color_filled> _data;

    mutable std::vector< std::vector<Node_weight_pair> > _single_pixel_contributors;
};




// store a list of color values and corresponding depths and weights (for example filling) per pixel
template <class Data>
class Accumulating_frame_buffer_without_queue : public Abstract_frame_buffer<Data>
{

    struct Data_depth_pixel
    {
        Data_depth_pixel() :
            depth(1e10f),
            filling_degree(0.0f),
            radius(0.0f),
            cone_angle(0.0f),
            is_full(false),
            node(NULL),
            debug_depth(0.0f)
        { }

        Data color;
        float depth;
        float filling_degree;
        float radius;
        float cone_angle;
        vector3d_t direction;
        bool is_full;
        Gi_point_base const* node;
        float debug_depth;

        friend bool operator<(Data_depth_pixel const& lhs, Data_depth_pixel const& rhs)
        {
            return (lhs.depth < rhs.depth);
        }
    };

public:
    virtual ~Accumulating_frame_buffer_without_queue() {}

    Accumulating_frame_buffer_without_queue()
    { }

    Accumulating_frame_buffer_without_queue(int size, int plane) : Abstract_frame_buffer<Data>(size, plane)
    {
        set_size(size);

        _use_visibility = true;
        _use_corrected_epsilon_z_buffer = false;
        _use_epsilon_z_buffer = false;
        _use_depth_dependant = false;
        _use_depth_modulation = false;
        _use_stochastic_visibility = false;

        _depth_modulation_max_depth_factor = 2.0f;
        _normalize = true;
    }

    virtual Abstract_frame_buffer<Data>* clone()
    {
        Accumulating_frame_buffer_without_queue* afb = new Accumulating_frame_buffer_without_queue(*this);

        return afb;
    }

    void set_parameters(Parameter_list const& parameters)
    {
        Abstract_frame_buffer<color_t>::set_parameters(parameters);

        _use_corrected_epsilon_z_buffer = false;
        _use_epsilon_z_buffer = false;
        _use_stochastic_visibility = false;

        _use_visibility = parameters["use_visibility"]->get_value<bool>();
        _use_depth_modulation = parameters["use_depth_modulation"]->get_value<bool>();
        _use_depth_dependant = parameters["use_depth_dependant"]->get_value<bool>();
        _normalize = parameters["normalize"]->get_value<bool>();
        _depth_modulation_max_depth_factor = parameters["depth_modulation_max_depth_factor"]->get_value<float>();


//        _use_visibility = true;
//        _use_depth_modulation = false;
//        _normalize = true;

        set_size(this->_resolution);
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

    float apply_z_correction(float const radius, float const distance_from_center)
    {
        return std::sqrt(radius * radius - distance_from_center * distance_from_center);
    }

    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info)
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        Data_depth_pixel & c = _data[serial_index];

        bool accepted = false;
        float weight = filling_degree;

        if (_use_visibility)
        {
            if (c.filling_degree < 1.0f)
            {
                if (c.filling_degree + weight > 1.0f)
                {
                    weight = 1.0f - c.filling_degree;
                }

                accepted = true;
            }
        }
        else if (_use_epsilon_z_buffer)
        {
            if (c.depth > 1e5f)
            {
                c.depth = point_info.depth + point_info.radius;
                // c.depth = point_info.depth + 0.01f;
            }

            if (point_info.depth < c.depth)
            {
                accepted = true;
            }
        }
//        else if (_use_corrected_epsilon_z_buffer)
//        {
//            if (c.depth > 1e5f)
//            {
//                c.depth = point_info.depth - apply_z_correction(point_info.radius, point_info.distance_from_center) + point_info.radius;
//                // c.depth = point_info.depth + 0.01f;
//            }

//            if (point_info.depth - apply_z_correction(point_info.radius, point_info.distance_from_center) < c.depth)
//            {
//                accepted = true;
//                c.debug_depth = point_info.distance_from_center / point_info.radius;
//            }
//        }
        else if (_use_stochastic_visibility)
        {
            if (!c.is_full && c.filling_degree + weight > 1.0f)
            {
                if (int(weight * 100000.0f) % 10 < 5.0f) // russian roulette
                {
                    accepted = true;
                }
                else
                {
                    c.is_full = true;
                }
            }

            if (!c.is_full)
            {
                accepted = true;
            }
        }
        else if (_use_depth_dependant)
        {
            if (c.filling_degree + weight > 1.0f)
            {
                if (c.depth > 1e5f) // first node being splat onto this pixel that completely fills the pixel
                {
                    c.depth  = point_info.depth;
                    c.radius = point_info.radius;
                }
                else
                {
                    float const d_0 = c.depth;
                    float const max_dist = _depth_modulation_max_depth_factor * c.radius;
                    float const depth_modulation = (point_info.depth - d_0) / max_dist;

                    assert(d_0 <= point_info.depth && depth_modulation >= 0.0f);

                    if (depth_modulation < 1.0f)
                    {
                        weight *= wendland_2_1(depth_modulation);
                    }
                    else
                    {
                        accepted = false;
                    }
                }
            }

            accepted = true;
        }
        else if (_use_depth_modulation)
        {
            accepted = true;

            if (c.depth > 1e5f) // first node being splat onto this pixel
            {
                c.depth  = point_info.depth;
                c.radius = point_info.radius;
            }
            else
            {
                float const d_0 = c.depth;
                float const max_dist = _depth_modulation_max_depth_factor * c.radius;
                //float const depth_modulation = (depth - d_0) / (2.0f * d_0 - d_0);
                // float const depth_modulation = (point_info.depth - d_0) / (1.1f * d_0 - d_0);
                float const depth_modulation = (point_info.depth - d_0) / max_dist;

                assert(d_0 <= point_info.depth && depth_modulation >= 0.0f);

                if (depth_modulation < 1.0f)
                {
                    weight = weight * wendland_2_1(depth_modulation);
                }
                else
                {
                    accepted = false;
                }
            }
        }
        else
        {
            accepted = true;
        }

        if (accepted && weight > 0.0001f)
        {
            c.filling_degree += weight;
            c.color += color * weight;

            // c.depth = (c.depth > 1e5f) ? point_info.depth * weight : c.depth + point_info.depth * weight;

#ifdef DEBUG
            if (_single_pixel_contributors.size() == 0)
            {
                _single_pixel_contributors.resize(Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
            }

            Node_weight_pair nwp(point_info, weight, color);
            _single_pixel_contributors[serial_index].push_back(nwp);

            point_info.splatted = true;
#endif

        }
    }

    virtual Data get_color(int const x, int const y, Debug_info * debug_info = NULL) const
    {
        int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

        Data_depth_pixel const& c = _data[serial_index];

        Color result_color = c.color;

        /*
        if (c.filling_degree > 0.001f)
        {
            result_color *= 1.0f / c.filling_degree;
        }
        else
        {
            result_color = color_t(0.0f);
        }
        */

        if (_normalize)
        {
            if (c.filling_degree > 0.1f)
            {
                result_color *= 1.0f / c.filling_degree;
            }
            else if (c.filling_degree > 0.0001f)
            {
                result_color *= 1.0f / c.filling_degree;
                result_color *= c.filling_degree / 0.1f;
            }
            else
            {
                result_color = color_t(0.0f);
            }
        }



        if (debug_info)
        {
//            for (std::tr1::unordered_set<yafaray::Gi_point_base const*>::const_iterator iter = _debug_gi_points.begin(); iter != _debug_gi_points.end(); ++iter)
//            {
//                debug_info->gi_points.insert(*iter);
//            }

            bool const debug_pixel = debug_info->cube_plane == Abstract_frame_buffer<Data>::_debug_plane && debug_info->cube_x == x && debug_info->cube_y == y;

            if (debug_pixel && _single_pixel_contributors.size() > 0)
            {
                debug_info->single_pixel_contributors = _single_pixel_contributors[serial_index];
            }

            // std::cout << "debug_gi_points: " << debug_gi_points.size() << std::endl;

            if (debug_info->debug_return_type == std::string("depth"))
            {
                result_color = color_t(c.debug_depth);

                // result_color = color_t(c.depth);
                // result_color = color_t(0.5f);
            }
            else if (debug_info->debug_return_type == std::string("filling_degree"))
            {
                result_color = color_t(c.filling_degree);
            }
        }

        return result_color;
    }


    static std::string name()
    {
        return "Accumulating_frame_buffer_without_queue";
    }

    static Accumulating_frame_buffer_without_queue * create()
    {
        return new Accumulating_frame_buffer_without_queue;
    }

    static Parameter_list get_parameters()
    {
        Parameter_list parameters;

        parameters.add_parameter(new Parameter("use_visibility",       true));
        parameters.add_parameter(new Parameter("use_depth_modulation", false));
        parameters.add_parameter(new Parameter("use_depth_dependant",  false));
        parameters.add_parameter(new Parameter("normalize",            true));
        parameters.add_parameter(new Parameter("depth_modulation_max_depth_factor", 2.0f, 0.0f, 4.0f));

        return parameters;
    }

protected:
    std::vector<Data_depth_pixel> _data;

    std::vector< std::vector<Node_weight_pair> > _single_pixel_contributors;

    bool _normalize;
    float _depth_modulation_max_depth_factor;

    bool _use_visibility;
    bool _use_corrected_epsilon_z_buffer;
    bool _use_epsilon_z_buffer;
    bool _use_depth_dependant;
    bool _use_depth_modulation;
    bool _use_stochastic_visibility;
};



static Class_parameter_registration< Abstract_frame_buffer<color_t>, Accumulating_frame_buffer_without_queue<color_t> > Accumulating_frame_buffer_without_queue_float(
        &Accumulating_frame_buffer_without_queue<color_t>::create,
        Accumulating_frame_buffer_without_queue<color_t>::get_parameters());


// REGISTER_CLASS_WITH_PARAMETERS("Accumulating_frame_buffer_without_queue", Abstract_frame_buffer<float>, Accumulating_frame_buffer_without_queue<float>);



// store a list of color values and corresponding depths and weights (for example filling) per pixel
template <class Data>
class Parameter_frame_buffer : public Abstract_frame_buffer<Data>
{
    struct Pixel
    {
        Pixel() :
            lobe(Mises_fisher_lobe<color_t>(vector3d_t(0.0f), 0.0f, Data(0.0f))),
            rp_to_node(vector3d_t(0.0f)),
            depth(1e10f),
            filling_degree(0.0f),
            radius(0.0f),
            node(NULL)
        { }

        Mises_fisher_lobe<color_t> lobe;
        vector3d_t rp_to_node;
        float depth;
        float filling_degree;
        float radius;
        Gi_point_base const* node;

        friend bool operator<(Pixel const& lhs, Pixel const& rhs)
        {
            return (lhs.depth < rhs.depth);
        }
    };

public:
    virtual ~Parameter_frame_buffer()
    {  }

    Parameter_frame_buffer(int size, int plane) : Abstract_frame_buffer<Data>(size, plane)
    {
        set_size(size);
    }

    static Abstract_frame_buffer<Data>* create(int resolution, int plane)
    {
        return new Parameter_frame_buffer(resolution, plane);
    }

    virtual Abstract_frame_buffer<Data>* clone()
    {
        Parameter_frame_buffer* afb = new Parameter_frame_buffer(*this);

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

    virtual void add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info);

    virtual Data get_color(int const x, int const y, Debug_info * debug_info = NULL) const;


protected:
    std::vector<Pixel> _data;

    mutable std::vector< std::vector<Node_weight_pair> > _single_pixel_contributors;
    std::vector< std::vector<Mises_fisher_lobe<color_t> > > _pixel_lobes;
};





/*
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


        if (debug_pixel)
        {
            std::cout << "Pixel acc color: " << accumulated_color << std::endl;
        }

        return accumulated_color;
    }
};

*/


template <class Data>
void Parameter_frame_buffer<Data>::add_point(int const x, int const y, Data const& color, float const filling_degree, Gi_point_info const& point_info)
{
    int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

    Pixel & c = _data[serial_index];

    Mises_Fisher_spherical_function<color_t> const* mf_sf = dynamic_cast< Mises_Fisher_spherical_function<color_t> *>(point_info.spherical_function->color);

    assert(mf_sf);

    if (mf_sf->get_lobes().size() == 0) return;

    assert(mf_sf->get_lobes().size() > 0);

    Mises_fisher_lobe<color_t> const& lobe = mf_sf->get_lobes()[0];

    assert(!std::isnan(lobe.mean_dir.x));


    bool const use_visibility = false;
    bool const use_depth_modulation = true;

    bool accepted = false;
    float weight = filling_degree;

    if (use_visibility)
    {
        if (c.filling_degree < 1.0f)
        {
            if (c.filling_degree + weight > 1.0f)
            {
                weight = 1.0f - c.filling_degree;
            }

            accepted = true;
        }
    }
    else if (use_depth_modulation)
    {
        accepted = true;

        if (c.depth > 1e5f) // first node being splat onto this pixel
        {
            c.depth  = point_info.depth;
            c.radius = point_info.radius;
        }
        else
        {
            float const d_0 = c.depth;
            float const max_dist = 2.0f * c.radius;
            //float const depth_modulation = (depth - d_0) / (2.0f * d_0 - d_0);
            // float const depth_modulation = (point_info.depth - d_0) / (1.1f * d_0 - d_0);
            float const depth_modulation = (point_info.depth - d_0) / max_dist;

            assert(d_0 <= point_info.depth && depth_modulation >= 0.0f);

            if (depth_modulation < 1.0f)
            {
                weight = weight * wendland_2_1(depth_modulation);
            }
            else
            {
                accepted = false;
            }
        }
    }
    else
    {
        accepted = true;
    }

    if (accepted && weight > 0.0001f)
    {
        c.filling_degree += weight;
        c.lobe.concentration += lobe.concentration * weight;
        c.lobe.mean_dir += lobe.mean_dir * weight;
        c.lobe.weight += lobe.weight * weight;
        c.rp_to_node += point_info.direction * weight;

#ifdef DEBUG
        if (_single_pixel_contributors.size() == 0)
        {
            _single_pixel_contributors.resize(Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
        }

        if (_pixel_lobes.size() == 0)
        {
            _pixel_lobes.resize(Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution);
        }

        Node_weight_pair nwp(point_info, weight, color);
        _single_pixel_contributors[serial_index].push_back(nwp);
        _pixel_lobes[serial_index].push_back(lobe);

        point_info.splatted = true;
#endif
    }
}


template <class Data>
Data Parameter_frame_buffer<Data>::get_color(int const x, int const y, Debug_info * debug_info) const
{
    int const serial_index = Abstract_frame_buffer<Data>::calc_serial_index(x, y);

    Pixel const& c = _data[serial_index];

    if (c.filling_degree < 0.001f) return Data(0.0f);

    Mises_fisher_lobe<color_t> lobe = c.lobe;

     // needs to be normalized even when buffer not completely filled
    lobe.weight *= (1.0f / c.filling_degree);
    lobe.concentration /= c.filling_degree;
    lobe.mean_dir.normalize();

    vector3d_t rp_to_node = c.rp_to_node;
    rp_to_node.normalize();


    int pixel_serial_index_cube =
            (y + Abstract_frame_buffer<Data>::_resolution_2) +
            (x + Abstract_frame_buffer<Data>::_resolution_2) * Abstract_frame_buffer<Data>::_resolution +
            Abstract_frame_buffer<Data>::_debug_plane * Abstract_frame_buffer<Data>::_resolution * Abstract_frame_buffer<Data>::_resolution;

    vector3d_t pixel_direction = Cube_static_data::cell_centers_luts[Abstract_frame_buffer<Data>::_resolution][pixel_serial_index_cube];
    pixel_direction.normalize();



    // std::cout << "x, y: " << x << ", " << y << " pixel_dir: " << pixel_direction << " color: " << lobe.evaluate(-pixel_direction) << " rp_to_node " << rp_to_node << " color: " << lobe.evaluate(-rp_to_node) << std::endl;

    // return lobe.evaluate(-rp_to_node) * c.filling_degree;
    // return lobe.evaluate(-pixel_direction) * c.filling_degree;
    color_t result_color = lobe.evaluate(-pixel_direction);

    if (debug_info)
    {
//        for (std::tr1::unordered_set<yafaray::Gi_point_base const*>::const_iterator iter = _debug_gi_points.begin(); iter != _debug_gi_points.end(); ++iter)
//        {
//            debug_info->gi_points.insert(*iter);
//        }

        bool const debug_pixel = debug_info->cube_plane == Abstract_frame_buffer<Data>::_debug_plane && debug_info->cube_x == x && debug_info->cube_y == y;

        if (debug_pixel && _single_pixel_contributors.size() > 0)
        {
            debug_info->single_pixel_contributors = _single_pixel_contributors[serial_index];

            std::cout << "Parameter_frame_buffer::get_color(): lobe " << lobe << " -rp to node dir: " << -rp_to_node <<
                         " -pixel_direction: " << -pixel_direction << " color: " << result_color << std::endl;
        }

        if (debug_pixel && _pixel_lobes.size() > 0)
        {
            for (std::size_t i = 0; i < _pixel_lobes[serial_index].size(); ++i)
            {
                std::cout << "Parameter_frame_buffer::get_color(): lobe " << i << " " << _pixel_lobes[serial_index][i] << std::endl;
            }
        }

        // std::cout << "debug_gi_points: " << debug_gi_points.size() << std::endl;

        // std::cout << "Parameter_frame_buffer::get_color(): debug_type: " << debug_info->debug_return_type << std::endl;

        if (debug_info->debug_return_type == std::string("mean_dir"))
        {
            result_color.R = lobe.mean_dir.x * 0.5f + 0.5f;
            result_color.G = lobe.mean_dir.y * 0.5f + 0.5f;
            result_color.B = lobe.mean_dir.z * 0.5f + 0.5f;
        }
        else if (debug_info->debug_return_type == std::string("filling_degree"))
        {
            result_color = color_t(c.filling_degree);
        }
    }

    return result_color;
}



__END_YAFRAY

#endif // FRAME_BUFFER_H
