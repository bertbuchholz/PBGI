#ifndef CUBERASTERBUFFER_H
#define CUBERASTERBUFFER_H

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

#include <bert/shared/Parameter.h>

#include <core_api/background.h>
#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <core_api/surface.h>
#include <utilities/Debug_info.h>
#include <utilities/sample_utils.h>
#include <utilities/mcqmc.h>
#include <integrators/Gi_point_info.h>

__BEGIN_YAFRAY

template <class Data>
class Abstract_frame_buffer;


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

inline bool is_in_range(float const min, float const max, float const val1, float const val2)
{
    if      (val1 / val2 < min) return false;
    else if (val1 / val2 > max) return false;
    return true;
}


inline bool is_in_range(float const min, float const max, color_t const& val1, color_t const& val2)
{
    if      (val1.R / val2.R < min) return false;
    else if (val1.R / val2.R > max) return false;
    if      (val1.G / val2.G < min) return false;
    else if (val1.G / val2.G > max) return false;
    if      (val1.B / val2.B < min) return false;
    else if (val1.B / val2.B > max) return false;
    return true;
}

template <class T>
T interpolate(T const& i1, T const& i2, float const offset)
{
    assert(offset >= 0.0f && offset <= 1.0f);
    return i1 * (1.0f - offset) + i2 * offset;
}

int get_longest_axis(vector3d_t const& dir);
float pyramid_solid_angle(std::vector<vector3d_t> const& corners);
float next_raster_in_positive_direction(float u);
float round(float const d);
bool sphere_line_intersection(vector3d_t const& sphere_center, float const radius, vector3d_t const& line_dir);
bool linePlaneIntersection(vector3d_t const& dir, vector3d_t const& n, float & alpha);
bool linePlaneIntersection(vector3d_t const& p0, vector3d_t const& p1, vector3d_t const& v0, vector3d_t const& n, float & alpha);


//class Spherical_node_representation;

//struct Gi_point_info
//{
//    enum Type { Node = 0, Far_surfel, Near_surfel };

//    Type type;
//    color_t color;
//    vector3d_t position;
//    vector3d_t receiver_position;
//    float depth;
//    vector3d_t direction; // from the receiving point to the GI point
//    float solid_angle;
//    vector3d_t disc_normal;
//    float radius;
//    float weight; // used for mixing when doing variational splatting
//    Spherical_node_representation const* spherical_function;
//};

struct Cube_cell
{
    int plane;
    int pos[2];
    float offset[2]; // 0..1

    Cube_cell() :
        plane(-1)
    {
        pos[0] = 0.0f;
        pos[1] = 0.0f;
        offset[0] = 0.0f;
        offset[1] = 0.0f;
    }

    Cube_cell & operator=(Cube_cell const& c)
    {
        plane = c.plane;
        pos[0] = c.pos[0];
        pos[1] = c.pos[1];
        offset[0] = c.offset[0];
        offset[1] = c.offset[1];

        return *this;
    }

    friend bool operator==(Cube_cell const& lhs, Cube_cell const& rhs)
    {
        return lhs.plane == rhs.plane && lhs.pos[0] == rhs.pos[0] && lhs.pos[1] == rhs.pos[1];
    }

    friend bool operator<(Cube_cell const& lhs, Cube_cell const& rhs)
    {
        if (lhs.plane < rhs.plane) return true;
        if (lhs.plane == rhs.plane)
        {
            if (lhs.pos[0] < rhs.pos[0]) return true;
            if (lhs.pos[0] == rhs.pos[0])
            {
                if (lhs.pos[1] < rhs.pos[1]) return true;
            }
        }

        return false;
    }

    friend std::ostream & operator<< (std::ostream & s, Cube_cell const& c)
    {
        s << c.plane << " " << c.pos[0] << " " << c.pos[1];
        return s;
    }
};


struct Cube_static_data
{
    static std::vector< std::vector<float> > solid_angle_luts; // [buffer resolution][cell serial index]
    static std::vector< std::vector<Cube_cell> > cube_cells_luts;
    static std::vector< std::vector<vector3d_t> > cell_centers_luts;

    static const vector3d_t cube_plane_normals[];

    static int const corresponding_axis_0[];
    static int const corresponding_axis_1[];

    static int get_corresponding_axis(int const axis, int const plane);

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


template <class Data>
class Cube_raster_buffer
{
public:


    typedef vector3d_t Point;
    typedef color_t Color;

    // static std::map<std::string, Frame_buffer_create_function> frame_buffer_create_map;

    Cube_raster_buffer();

    Cube_raster_buffer & operator= (Cube_raster_buffer const& cbr);

    virtual ~Cube_raster_buffer();

    void clear();

    // void setup_simple(int resolution);
    void setup_surfel_buffer(int resolution);

    /*
    void setup(Type fb_type,
               int resolution,
               Splat_type const node_splat_type,
               Splat_type const surfel_far_splat_type,
               Splat_type const surfel_near_splat_type);
               */

    int get_resolution() const;

    Cube_cell get_cell(Point const& dir, bool * ok = NULL) const;

    inline int get_serial_index(Cube_cell const& c) const;

    Point const& get_cell_center(int const serial_index) const;
    Point const& get_cell_center(Cube_cell const& c) const;
    Point        get_cell_direction(int const serial_index) const;
    Point        get_cell_direction(Cube_cell const& c) const;

    std::vector<Cube_cell> const& get_cube_cells() const;

    std::vector<Point> get_cell_corners(Cube_cell const& c) const;

    /*
    void add_point_single_pixel(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_stochastic_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_stochastic_node_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_aa_square(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_solid_angle_rays(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point(Gi_point_info const& point_info, GiPoint const* gi_point = NULL);
    */

    void set_data(Cube_cell const& c, Data const& data);
    void add_data(Cube_cell const& c, Data const& data);
    Data get_data(Cube_cell const& c, Debug_info * debug_info = NULL) const;
    Data get_data(int const plane, int const x, int const y, Debug_info * debug_info = NULL) const;
    Data get_data_interpolated(Cube_cell const& c, Debug_info * debug_info = NULL) const;

    Abstract_frame_buffer<Data> * get_component(int const index) const;
    std::vector<float> component_to_vector(int const component_index, bool const use_first_component_only = false) const;
    void from_component_vector(int const component_index, std::vector<float> const& v, bool const use_first_component_only = false);

    float calc_solid_angle(Cube_cell const& c) const;
    float get_solid_angle (int const serial_index) const;
    float get_solid_angle (Cube_cell const& c) const;
    float get_solid_angle (Point const& dir)   const;

    Data get_diffuse(Point const& normal) const;

    float get_total_energy();
    float get_non_zero_area();

    void get_debug_info(Debug_info & debug_info);

    Cube_raster_buffer<Data> blur();

    Point calc_cell_center(Cube_cell const& c) const;
    Point calc_cell_position_with_offset(Cube_cell const& c) const;

    bool is_black() const;

    void print() const;
    Data find_maximum() const;
    Data calc_total_energy() const;
    Data integrate() const;
    Data integrate(Cube_raster_buffer<Data> const& weights) const;

protected:
    Abstract_frame_buffer<Data>* buffers[6];
    mutable float _total_energy;
    mutable float _non_zero_area;
    int _resolution;
    int _resolution_2;

    mutable bool _get_color_done;
};


class Splat_cube_raster_buffer : public Cube_raster_buffer<color_t>
{
public:
    enum Splat_type { Single_pixel = 0, Disc_tracing, AA_square, Gaussian_splat, SA_tracing, Stocastic_tracing, Stocastic_node_tracing };
    enum Buffer_type { Simple = 0, Accumulating, Distance_weighted, Parameter };

    static std::map<std::string, Splat_type> enum_splat_type_map;
    static std::map<std::string, Buffer_type> enum_type_map;


    void setup_empty();

    void setup(
            Buffer_type fb_type,
            int resolution,
            Splat_type const node_splat_type,
            Splat_type const surfel_far_splat_type,
            Splat_type const surfel_near_splat_type);

    void setup(Parameter_list const& parameters);

    Splat_cube_raster_buffer & operator= (Splat_cube_raster_buffer const& cbr);

    void add_point_single_pixel(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_stochastic_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_stochastic_node_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_aa_square(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_gaussian_splat(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_solid_angle_rays(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point(Gi_point_info const& point_info, GiPoint const* gi_point = NULL);

    void add_background(background_t * background);

    Color get_brdf_response(renderState_t & state, surfacePoint_t const& receiving_point, vector3d_t const& wo) const;

    Splat_cube_raster_buffer blur() const;


private:
    typedef void (Splat_cube_raster_buffer::*Add_point_function_ptr)(Gi_point_info const& point_info, GiPoint const* gi_point);

    std::vector<Add_point_function_ptr> add_point_function_map;
};



__END_YAFRAY

#include <utilities/Frame_buffer.h>

__BEGIN_YAFRAY







template <class Data>
Cube_raster_buffer<Data>::Cube_raster_buffer() :
    _get_color_done(false)
{
    for (int i = 0; i < 6; ++i)
    {
        buffers[i] = NULL;
    }
}

template <class Data>
Cube_raster_buffer<Data> & Cube_raster_buffer<Data>::operator= (Cube_raster_buffer const& cbr)
{
    if (this != &cbr)
    {
        for (int i = 0; i < 6; ++i)
        {
            buffers[i] = cbr.buffers[i]->clone();
        }

        _total_energy          = cbr._total_energy;
        _resolution            = cbr._resolution;
        _resolution_2          = cbr._resolution_2;
        // _cell_centers          = cbr._cell_centers;
        // _cube_cells            = cbr._cube_cells;
        _get_color_done        = cbr._get_color_done;
    }

    return *this;
}

template <class Data>
Cube_raster_buffer<Data>::~Cube_raster_buffer()
{
    for (int i = 0; i < 6; ++i)
    {
        if (buffers[i])
        {
            buffers[i]->clear();
        }
    }

    for (int i = 0; i < 6; ++i)
    {
        delete buffers[i];
    }
}

template <class Data>
void Cube_raster_buffer<Data>::clear()
{
    for (int i = 0; i < 6; ++i)
    {
        buffers[i]->clear();
    }

    // _cell_centers.clear();
    // _cube_cells.clear();

    // add_point_function_map.clear();
}

template <class Data>
void Cube_raster_buffer<Data>::setup_surfel_buffer(int resolution)
{

    _resolution = resolution;
    _resolution_2 = _resolution / 2;

    /*
    add_point_function_map.resize(3);
    add_point_function_map[Gi_point_info::Node]        = Cube_raster_buffer::add_point_single_pixel;
    add_point_function_map[Gi_point_info::Far_surfel]  = Cube_raster_buffer::add_point_single_pixel;
    add_point_function_map[Gi_point_info::Near_surfel] = Cube_raster_buffer::add_point_single_pixel;
    */

    for (int i = 0; i < 6; ++i)
    {
        buffers[i] = new Simple_single_value_frame_buffer<Data>(_resolution, i);
    }
}

/*
template <class Data>
void Cube_raster_buffer<Data>::setup_simple(int resolution)
{
    setup(Simple, resolution, Single_pixel, Single_pixel, Single_pixel);
}
*/

/*
template <class Data>
void Cube_raster_buffer<Data>::setup(Type fb_type,
                               int resolution,
                               Splat_type const node_splat_type,
                               Splat_type const surfel_far_splat_type,
                               Splat_type const surfel_near_splat_type)
{
    _resolution = resolution;
    _resolution_2 = _resolution / 2;

    for (int i = 0; i < 6; ++i)
    {
        // buffers[i] = frame_buffer_create_map[() // FIXME: missing id string, need some restructuring

        if (fb_type == Simple)
        {
            buffers[i] = new Simple_frame_buffer_without_queue<Data>(_resolution, i);
        }
        else if (fb_type == Accumulating)
        {
            buffers[i] = new Accumulating_frame_buffer_without_queue<Data>(_resolution, i);
        }
        else if (fb_type == Distance_weighted)
        {
            buffers[i] = new Grouping_frame_buffer_without_queue<Data>(_resolution, i);
        }
        else
        {
            std::cout << "Unknown FB type ..." << std::endl;
            assert(false);
        }
    }
}
*/

template <class Data>
int Cube_raster_buffer<Data>::get_resolution() const
{
    return _resolution;
}


template <class Data>
Cube_cell Cube_raster_buffer<Data>::get_cell(Point const& dir, bool * ok) const
{
    Cube_cell result;

    int longest_axis = get_longest_axis(dir);

    int plane_index = longest_axis * 2;

    float longest_extent = dir[longest_axis];

    if (longest_extent < 0)
    {
        ++plane_index;
        longest_extent = -longest_extent;
    }

    int const axis_0 = Cube_static_data::corresponding_axis_0[plane_index] / 2;
    int const axis_1 = Cube_static_data::corresponding_axis_1[plane_index] / 2;

    int const sign_axis_0 = (Cube_static_data::corresponding_axis_0[plane_index] % 2 == 0) ? 1 : -1;
    int const sign_axis_1 = (Cube_static_data::corresponding_axis_1[plane_index] % 2 == 0) ? 1 : -1;

    result.plane = plane_index;

    assert(longest_extent != 0.0f);

    float inv_longest_extent = 1.0f / longest_extent;

    float const axis_0_pos = dir[axis_0] * inv_longest_extent * _resolution_2 * sign_axis_0;
    float const axis_1_pos = dir[axis_1] * inv_longest_extent * _resolution_2 * sign_axis_1;

    result.pos[0] = std::floor(axis_0_pos);
    result.pos[1] = std::floor(axis_1_pos);

    result.pos[0] = into_range(-_resolution_2, _resolution_2 - 1, result.pos[0]);
    result.pos[1] = into_range(-_resolution_2, _resolution_2 - 1, result.pos[1]);

    result.offset[0] = axis_0_pos - result.pos[0] - 0.5f;
    result.offset[1] = axis_1_pos - result.pos[1] - 0.5f;

    if (ok)
    {
        *ok = true;

        if (result.pos[0] < -_resolution_2 || result.pos[0] >= _resolution_2) *ok = false;
        if (result.pos[1] < -_resolution_2 || result.pos[1] >= _resolution_2) *ok = false;
    }

    return result;
}


template <class Data>
inline int Cube_raster_buffer<Data>::get_serial_index(Cube_cell const& c) const
{
    int serial =
            (c.pos[1] + _resolution_2) +
            (c.pos[0] + _resolution_2) * _resolution +
            c.plane * _resolution * _resolution;

    return serial;
}


template <class Data>
typename Cube_raster_buffer<Data>::Point const& Cube_raster_buffer<Data>::get_cell_center(Cube_cell const& c) const
{
    int serial = get_serial_index(c);

    return get_cell_center(serial);
}


template <class Data>
inline typename Cube_raster_buffer<Data>::Point const& Cube_raster_buffer<Data>::get_cell_center(int const serial_index) const
{
    return Cube_static_data::cell_centers_luts[_resolution][serial_index];
}


template <class Data>
typename Cube_raster_buffer<Data>::Point Cube_raster_buffer<Data>::get_cell_direction(int const serial_index) const
{
    Point dir = get_cell_center(serial_index);
    dir.normalize();
    return dir;
}


template <class Data>
typename Cube_raster_buffer<Data>::Point Cube_raster_buffer<Data>::get_cell_direction(Cube_cell const& c) const
{
    Point dir = get_cell_center(c);
    dir.normalize();
    return dir;
}


template <class Data>
std::vector<Cube_cell> const& Cube_raster_buffer<Data>::get_cube_cells() const
{
    //return _cube_cells;
    return Cube_static_data::cube_cells_luts[_resolution];
}


template <class Data>
std::vector<typename Cube_raster_buffer<Data>::Point> Cube_raster_buffer<Data>::get_cell_corners(Cube_cell const& c) const
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




template <class Data>
void Cube_raster_buffer<Data>::set_data(Cube_cell const& c, Data const& data)
{
    buffers[c.plane]->set_color(c.pos[0], c.pos[1], data);
}


template <class Data>
void Cube_raster_buffer<Data>::add_data(Cube_cell const& c, Data const& data)
{
    Data current_color = get_data(c);
    set_data(c, current_color + data);
}


template <class Data>
Data Cube_raster_buffer<Data>::get_data(Cube_cell const& c, Debug_info * debug_info) const
{
    return buffers[c.plane]->get_color(c.pos[0], c.pos[1], debug_info);
}


template <class Data>
Data Cube_raster_buffer<Data>::get_data(int const plane, int const x, int const y, Debug_info * debug_info) const
{
    return buffers[plane]->get_color(x, y, debug_info);
}


/*
void get_wrapped_coordinate(int const plane, int const coordinate[2], int const resolution, int * wrapped_coordinate)
{
    wrapped_coordinate[0] = coordinate[0];
    wrapped_coordinate[1] = coordinate[1];

    if (coordinate[0] >= -resolution && coordinate[0] < resolution &&
            coordinate[1] >= -resolution && coordinate[1] < resolution)
    {
        return;
    }

    bool const is_s = coordinate[0] < -resolution || coordinate[0] >= resolution;
    int const changed_coordinate_index = is_s ? 0 : 1;

    int const sign = coordinate[changed_coordinate_index] > 0 ? 1 : -1;

    if (plane == 4) // +z
    {


        if (is_s)
        {
            wrapped_coordinate[0] = sign * 2 * resolution + coordinate[0];
        }
        // keep t as is
    }
    else if (plane / 2 == 0) // +x, -x
    {
        if (is_s) // wrap to z
        {
            wrapped_coord = sign * 2 * resolution + coordinate;
        }
        else // wrap to y
        {

        }
    }
    else if (plane / 2 == 1) // +y, -y
    {
        if (is_s) // wrap to z
        {
            wrapped_coord = sign * 2 * resolution + coordinate;
        }
        else // wrap to y
        {

        }
    }

}

*/


template <class Data>
Data Cube_raster_buffer<Data>::get_data_interpolated(Cube_cell const& c, Debug_info * debug_info) const
{
    // return buffers[c.plane]->get_color_interpolated(c.pos[0] + c.offset[0], c.pos[1] + c.offset[1], debug_info);

    Cube_cell cell = c;

    Data colors[4];

    int const shift_x = (c.offset[0] < 0.0f) ? -1 : 0;
    int const shift_y = (c.offset[1] < 0.0f) ? -1 : 0;

    for (int u = 0; u < 2; ++u)
    {
        for (int v = 0; v < 2; ++v)
        {
            cell.pos[0] = c.pos[0] + u + shift_x;
            cell.pos[1] = c.pos[1] + v + shift_y;

            //Point new_center = get_cell_center(cell);
            Point new_center = calc_cell_center(cell);
            new_center.normalize();

            Cube_cell new_cell = get_cell(new_center);
            colors[u + v * 2] = get_data(new_cell);

            if (debug_info)
            {
                std::cout << "Cube_raster_buffer::get_color_interpolated(): new_cell: " <<  new_cell << " index: " << (u + v * 2) << std::endl;
            }
        }
    }

    float offset_x = (c.offset[0] < 0.0f) ? 1.0f + c.offset[0] : c.offset[0];
    float offset_y = (c.offset[1] < 0.0f) ? 1.0f + c.offset[1] : c.offset[1];

    if (debug_info)
    {
        // std::cout << "Cube_raster_buffer::get_color_interpolated(): offset: " << offset_x << " " << offset_y << std::endl;
    }

    offset_x = std::max(0.0f, std::min(offset_x, 1.0f));
    offset_y = std::max(0.0f, std::min(offset_y, 1.0f));

    Data const c1_2 = interpolate(colors[0], colors[1], offset_x);
    Data const c3_4 = interpolate(colors[2], colors[3], offset_x);

    Data const result = interpolate(c1_2, c3_4, offset_y);

    if (debug_info)
    {
        std::cout << "Cube_raster_buffer::get_color_interpolated(): offset: " << offset_x << " " << offset_y << " result: " << result << std::endl;
    }


    return result;
}

/*
    Color get_color_interpolated(Cube_cell const& c, Debug_info * debug_info = NULL) const
    {
        // return buffers[c.plane]->get_color_interpolated(c.pos[0] + c.offset[0], c.pos[1] + c.offset[1], debug_info);
        Color result;

        float const x = c.pos[0] + c.offset[0];
        float const y = c.pos[1] + c.offset[1];

        int x0 = int(std::floor(x - 0.5f));
        int y0 = int(std::floor(y - 0.5f));

        int x1 = x0 + 1;
        int y1 = y0 + 1;

        if (x0 >= -_resolution_2 && y0 >= -_resolution_2 &&
            x1 <   _resolution_2 && y1 <   _resolution_2)
        {
            // in the middle, no wrapping necessary

            color_t const c1 = get_color(c.plane, x0, y0);
            color_t const c2 = get_color(c.plane, x1, y0);
            color_t const c3 = get_color(c.plane, x0, y1);
            color_t const c4 = get_color(c.plane, x1, y1);

            float offset_x = x - x0 - 0.5f;
            float offset_y = y - y0 - 0.5f;

            offset_x = std::max(0.0f, std::min(offset_x, 1.0f));
            offset_y = std::max(0.0f, std::min(offset_y, 1.0f));

            color_t const c1_2 = interpolate(c1, c2, offset_x);
            color_t const c3_4 = interpolate(c3, c4, offset_x);

            return interpolate(c1_2, c3_4, offset_y);
        }
        else if (x0 < -_resolution_2 && y0 >= -_resolution_2 && y1 < _resolution_2)
        {
            // y in the middle, x0 needs to be wrapped left
            int wrapped_x0    = -1;
            int wrapped_plane = -1;

            color_t const c1 = get_color(wrapped_plane, wrapped_x0, y0);
            color_t const c2 = get_color(plane,         x1,         y0);
            color_t const c3 = get_color(wrapped_plane, wrapped_x0, y1);
            color_t const c4 = get_color(plane,         x1,         y1);

            float offset_x = x - x0 - 0.5f;
            float offset_y = y - y0 - 0.5f;

            offset_x = std::max(0.0f, std::min(offset_x, 1.0f));
            offset_y = std::max(0.0f, std::min(offset_y, 1.0f));

            color_t const c1_2 = interpolate(c1, c2, offset_x);
            color_t const c3_4 = interpolate(c3, c4, offset_x);

            return interpolate(c1_2, c3_4, offset_y);
        }
        else
        {


        }



        if (!(x1 < _resolution_2 && y1 < _resolution_2))
        {
            std::cout << "_resolution_2: " << _resolution_2 << " x: " << x << " y: " << y << " x_0: " << x0 << " y_0: " << y0 << std::endl;
            assert(false);
        }



        return result;
    }
    */


template <class Data>
Abstract_frame_buffer<Data> * Cube_raster_buffer<Data>::get_component(int const index) const
{
    return buffers[index];
}





template <class Data>
void Cube_raster_buffer<Data>::from_component_vector(int const component_index, std::vector<float> const& v, bool const use_first_component_only)
{
    Cube_cell c;

    unsigned int i = 0;

    c.plane = component_index;

    for (int x = -_resolution_2; x < _resolution_2; ++x)
    {
        c.pos[0] = x;

        for (int y = -_resolution_2; y < _resolution_2; ++y)
        {
            c.pos[1] = y;

            // FIXME: vector of float to different types (color, single float etc.)
            // set_data(c, v[i]);

            ++i;
        }
    }
}


template <class Data>
float Cube_raster_buffer<Data>::calc_solid_angle(Cube_cell const& c) const
{
    std::vector<vector3d_t> corners = get_cell_corners(c);
    float const solid_angle = pyramid_solid_angle(corners);
    return solid_angle;
}


template <class Data>
float Cube_raster_buffer<Data>::get_solid_angle(Cube_cell const& c) const
{
    int const serial_index = get_serial_index(c);
    return Cube_static_data::solid_angle_luts[_resolution][serial_index];
}


template <class Data>
float Cube_raster_buffer<Data>::get_solid_angle(int const serial_index) const
{
    return Cube_static_data::solid_angle_luts[_resolution][serial_index];
}


template <class Data>
float Cube_raster_buffer<Data>::get_solid_angle(Point const& dir) const
{
    Cube_cell c = get_cell(dir);
    return get_solid_angle(c);
}


template <class Data>
Data Cube_raster_buffer<Data>::get_diffuse(Point const& normal) const
{
    assert(!_get_color_done);
    _get_color_done = true;

    Color diffuse(0.0f);
    Color total(0.0f);

    _non_zero_area = 0;

    Cube_cell c;

    float debug_weight_sum = 0.0f;

    for (int i = 0; i < 6; ++i)
    {
        c.plane = i;

        for (int x = -_resolution_2; x < _resolution_2; ++x)
        {
            c.pos[0] = x;

            for (int y = -_resolution_2; y < _resolution_2; ++y)
            {
                c.pos[1] = y;

                int const serial_index = get_serial_index(c);
                Point cell_dir = get_cell_direction(serial_index);

                float const cos_sp_normal_cell_dir = cell_dir * normal; // assumes the receiver is a lambertian reflector
                assert(cos_sp_normal_cell_dir <= 1.0001f);
                if (cos_sp_normal_cell_dir < 0.001f) continue;

                float const solid_angle = get_solid_angle(serial_index);
                float const weight = solid_angle; // normalization not to 1 but the hemisphere area (2.0f * M_PI) which it already is
                assert(weight >= 0.0f);

                Color color = get_data(c);

//                total += color; // * weight;

//                if (color.energy() > 0.0000001f)
//                {
//                    // std::cout << "non_zero:" << solid_angle << std::endl;
//                    _non_zero_area += solid_angle;
//                }


                debug_weight_sum += weight;

                diffuse += color * cos_sp_normal_cell_dir * weight;
            }
        }
    }

    // not sure about the normalization, the solid angle weights sum up to 2 pi (over the hemisphere)
    // so it seems reasonable to take the factor out again
//    diffuse *= 1.0f / (2.0f * M_PI);

    // std::cout << "get_diffuse(): " << debug_weight_sum << std::endl;

    _total_energy = total.energy();

    return diffuse;
}


template <class Data>
float Cube_raster_buffer<Data>::get_total_energy()
{
    return _total_energy;
}


template <class Data>
float Cube_raster_buffer<Data>::get_non_zero_area()
{
    return _non_zero_area;
}


template <class Data>
void Cube_raster_buffer<Data>::get_debug_info(Debug_info & debug_info)
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

                get_data(c, &debug_info);
            }
        }
    }
}


template <class Data>
typename Cube_raster_buffer<Data>::Point Cube_raster_buffer<Data>::calc_cell_center(Cube_cell const& c) const
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

    int const axis_0 = Cube_static_data::corresponding_axis_0[c.plane] / 2;
    int const axis_1 = Cube_static_data::corresponding_axis_1[c.plane] / 2;

    int const sign_axis_0 = (Cube_static_data::corresponding_axis_0[c.plane] % 2 == 0) ? 1 : -1;
    int const sign_axis_1 = (Cube_static_data::corresponding_axis_1[c.plane] % 2 == 0) ? 1 : -1;

    result[axis_0] = ((c.pos[0] + 0.5f) * sign_axis_0) / float(_resolution_2);
    result[axis_1] = ((c.pos[1] + 0.5f) * sign_axis_1) / float(_resolution_2);

    return result;
}


template <class Data>
typename Cube_raster_buffer<Data>::Point Cube_raster_buffer<Data>::calc_cell_position_with_offset(Cube_cell const& c) const
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

    int const axis_0 = Cube_static_data::corresponding_axis_0[c.plane] / 2;
    int const axis_1 = Cube_static_data::corresponding_axis_1[c.plane] / 2;

    int const sign_axis_0 = (Cube_static_data::corresponding_axis_0[c.plane] % 2 == 0) ? 1 : -1;
    int const sign_axis_1 = (Cube_static_data::corresponding_axis_1[c.plane] % 2 == 0) ? 1 : -1;

    result[axis_0] = ((c.pos[0] + 0.5f + c.offset[0]) * sign_axis_0) / float(_resolution_2);
    result[axis_1] = ((c.pos[1] + 0.5f + c.offset[1]) * sign_axis_1) / float(_resolution_2);

    return result;
}



template <class Data>
std::vector<float> Cube_raster_buffer<Data>::component_to_vector(int const component_index, bool const use_first_component_only) const
{
    std::vector<float> result;

    Cube_cell c;

    c.plane = component_index;

    for (int x = -_resolution_2; x < _resolution_2; ++x)
    {
        c.pos[0] = x;

        for (int y = -_resolution_2; y < _resolution_2; ++y)
        {
            c.pos[1] = y;

            Color const& color = get_data(c);

            if (use_first_component_only)
            {
                result.push_back(color[0]);
            }
            else
            {
                result.push_back(color[0]);
                result.push_back(color[1]);
                result.push_back(color[2]);
            }
        }
    }

    return result;
}

template <class Data>
bool Cube_raster_buffer<Data>::is_black() const
{
    std::vector<Cube_cell> const& cells = get_cube_cells();

    Data threshold = Data(0.00001f);

    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        if (get_data(cells[i]) > threshold)
        {
            return false;
        }
    }

    return true;
}

template <class Data>
void Cube_raster_buffer<Data>::print() const
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

                Data const& d = get_data(c);

                std::cout << d << " ";
            }

            std::cout << std::endl;
        }

        std::cout << std::endl;
    }
}

template <class Data>
Data Cube_raster_buffer<Data>::find_maximum() const
{
    std::vector<Cube_cell> const& cells = get_cube_cells();

    Data max(-1e5f);

    for (unsigned int i = 0; i < cells.size(); ++i)
    {
        if (get_data(cells[i]) > max)
        {
            max = get_data(cells[i]);
        }
    }

    return max;
}

template <class Data>
Data Cube_raster_buffer<Data>::calc_total_energy() const
{
    std::vector<Cube_cell> const& cells = get_cube_cells();

    Data total(0.0f);

    for (unsigned int i = 0; i < cells.size(); ++i)
    {
        total += get_data(cells[i]);
    }

    return total;
}


template <class Data>
Data Cube_raster_buffer<Data>::integrate() const
{
    std::vector<Cube_cell> const& cells = get_cube_cells();

    Data total(0.0f);

    for (unsigned int i = 0; i < cells.size(); ++i)
    {
        total += get_data(cells[i]) * get_solid_angle(cells[i]);
    }

    return total;
}


template <class Data>
Data Cube_raster_buffer<Data>::integrate(Cube_raster_buffer<Data> const& weights) const
{
    assert(_resolution == weights.get_resolution());

    std::vector<Cube_cell> const& cells = get_cube_cells();

    Data total(0.0f);

    for (unsigned int i = 0; i < cells.size(); ++i)
    {
        total += get_data(cells[i]) * get_solid_angle(cells[i]) * weights.get_data(cells[i]);
    }

    return total;
}


template <class Data>
Cube_raster_buffer<Data> Cube_raster_buffer<Data>::blur()
{
    Cube_raster_buffer<Data> tmp_buffer;
    tmp_buffer.setup_surfel_buffer(_resolution);

    Cube_cell c;

    int const blur_size = 2;

    for (int i = 0; i < 6; ++i)
    {
        c.plane = i;

        for (int x = -_resolution_2; x < _resolution_2; ++x)
        {
            c.pos[0] = x;

            for (int y = -_resolution_2; y < _resolution_2; ++y)
            {
                c.pos[1] = y;

                Data d(0.0f);

                for (int u = -blur_size; u <= blur_size; ++u)
                {
                    for (int v = -blur_size; v <= blur_size; ++v)
                    {
                        Cube_cell c_tmp;
                        c_tmp.plane = i;
                        c_tmp.pos[0] = x + u;
                        c_tmp.pos[1] = y + v;

                        Point new_center = calc_cell_center(c_tmp);
                        new_center.normalize();

                        c_tmp = get_cell(new_center);
                        d += get_data(c_tmp);
                    }
                }

                //d *= 1.0f / (blur_size * blur_size);
                tmp_buffer.set_data(c, d);
            }
        }
    }

    /*
    std::vector<Cube_cell> const& cells = get_cube_cells();

    for (unsigned int i = 0; i < cells.size(); ++i)
    {
        // set_data(cells[i], tmp_buffer.get_data(cells[i]));
        set_data(cells[i], Data(0.0f));
    }
    */

    return tmp_buffer;
}




__END_YAFRAY

#endif // CUBERASTERBUFFER_H

