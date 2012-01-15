#ifndef CUBERASTERBUFFER_H
#define CUBERASTERBUFFER_H

#include <algorithm>
#include <vector>
#include <map>

#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <utilities/Debug_info.h>
#include <utilities/Frame_buffer.h>
#include <utilities/sample_utils.h>
#include <utilities/mcqmc.h>

__BEGIN_YAFRAY

class Spherical_function;

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

template <class T>
Abstract_frame_buffer* create(int resolution, int plane_index)
{
    return new T(resolution, plane_index);
}

typedef Abstract_frame_buffer* (*Frame_buffer_create_function)(int, int);

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
    Spherical_function * spherical_function;
};

struct Cube_cell
{
    int plane;
    int pos[2];
    float offset[2]; // 0..1

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
    enum Type { Simple = 0, Accumulating, Reweighting, Distance_weighted, Front_stacked, Directional_stacked };
    enum Splat_type { Single_pixel = 0, Disc_tracing, AA_square, SA_tracing, Stocastic_tracing, Stocastic_node_tracing };

    typedef vector3d_t Point;
    typedef color_t Color;

    static std::map<std::string, Splat_type> enum_splat_type_map;
    static std::map<std::string, Type> enum_type_map;
    static std::map<std::string, Frame_buffer_create_function> frame_buffer_create_map;

    static const vector3d_t cube_plane_normals[];

    Cube_raster_buffer();

    Cube_raster_buffer & operator= (Cube_raster_buffer const& cbr);

    ~Cube_raster_buffer();

    void clear();

    void setup_simple(int resolution);
    void setup_surfel_buffer(int resolution);

    void setup(Type fb_type,
               int resolution,
               Splat_type const node_splat_type,
               Splat_type const surfel_far_splat_type,
               Splat_type const surfel_near_splat_type);

    int get_resolution() const;

    Cube_cell get_cell(Point const& dir, bool * ok = NULL) const;

    inline int get_serial_index(Cube_cell const& c) const;

    Point const& get_cell_center(int const serial_index) const;
    Point const& get_cell_center(Cube_cell const& c) const;
    Point        get_cell_direction(int const serial_index) const;
    Point        get_cell_direction(Cube_cell const& c) const;

    std::vector<Cube_cell> const& get_cube_cells() const;

    std::vector<Point> get_cell_corners(Cube_cell const& c) const;

    void add_point_single_pixel(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_stochastic_disc_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_stochastic_node_tracing(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_aa_square(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point_solid_angle_rays(Gi_point_info const& point_info, GiPoint const* node = NULL);
    void add_point(Gi_point_info const& point_info, GiPoint const* gi_point = NULL);

    void set_color(Cube_cell const& c, Color const& color);
    void add_color(Cube_cell const& c, Color const& color);
    Color get_color(Cube_cell const& c, Debug_info * debug_info = NULL) const;
    Color get_color(int const plane, int const x, int const y, Debug_info * debug_info = NULL) const;
    Color get_color_interpolated(Cube_cell const& c, Debug_info * debug_info = NULL) const;

    Abstract_frame_buffer * get_component(int const index) const;
    std::vector<float> component_to_vector(int const component_index, bool const use_first_component_only = false) const;
    void from_component_vector(int const component_index, std::vector<float> const& v, bool const use_first_component_only = false);

    float calc_solid_angle(Cube_cell const& c) const;
    float get_solid_angle (int const serial_index) const;
    float get_solid_angle (Cube_cell const& c) const;
    float get_solid_angle (Point const& dir)   const;

    Color get_diffuse(Point const& normal) const;

    float get_total_energy();
    float get_non_zero_area();

    void get_debug_info(Debug_info & debug_info);

    void blur();

    static int get_corresponding_axis(int const axis, int const plane);

    Point calc_cell_center(Cube_cell const& c) const;
    Point calc_cell_position_with_offset(Cube_cell const& c) const;

private:
    typedef void (Cube_raster_buffer::*Add_point_function_ptr)(Gi_point_info const& point_info, GiPoint const* gi_point);

    Abstract_frame_buffer* buffers[6];
    mutable float _total_energy;
    mutable float _non_zero_area;
    int _resolution;
    int _resolution_2;

    std::vector<Add_point_function_ptr> add_point_function_map;

    mutable bool _get_color_done;

    static int const corresponding_axis_0[];
    static int const corresponding_axis_1[];

    static std::vector< std::vector<float> > solid_angle_luts; // [buffer resolution][cell serial index]
    static std::vector< std::vector<Cube_cell> > cube_cells_luts;
    static std::vector< std::vector<Point> > cell_centers_luts;
};

__END_YAFRAY

#endif // CUBERASTERBUFFER_H

