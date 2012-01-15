#include <map>

#include <utilities/CubeRasterBuffer.h>

#include <utilities/Frame_buffer.h>

#include <utilities/interpolation.h>

#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

std::map<std::string, Cube_raster_buffer::Splat_type> create_splat_type_map()
{
    std::map<std::string, Cube_raster_buffer::Splat_type> m;
    m["Single_pixel"] = Cube_raster_buffer::Single_pixel;
    m["Disc_tracing"] = Cube_raster_buffer::Disc_tracing;
    m["AA_square"] = Cube_raster_buffer::AA_square;
    m["SA_tracing"] = Cube_raster_buffer::SA_tracing;
    m["Stocastic_tracing"] = Cube_raster_buffer::Stocastic_tracing;
    m["Stocastic_node_tracing"] = Cube_raster_buffer::Stocastic_node_tracing;
    return m;
}

std::map<std::string, Cube_raster_buffer::Splat_type> Cube_raster_buffer::enum_splat_type_map = create_splat_type_map();

std::map<std::string, Cube_raster_buffer::Type> create_type_map()
{
    std::map<std::string, Cube_raster_buffer::Type> m;
    m["Simple"] = Cube_raster_buffer::Simple;
    m["Accumulating"] = Cube_raster_buffer::Accumulating;
    m["Reweighting"] = Cube_raster_buffer::Reweighting;
    m["Distance_weighted"] = Cube_raster_buffer::Distance_weighted;
    m["Front_stacked"] = Cube_raster_buffer::Front_stacked;
    m["Directional_stacked"] = Cube_raster_buffer::Directional_stacked;
    return m;
}

std::map<std::string, Cube_raster_buffer::Type> Cube_raster_buffer::enum_type_map = create_type_map();

std::map<std::string, Frame_buffer_create_function> create_fb_create_map()
{
    std::map<std::string, Frame_buffer_create_function> m;
    m["Simple"] = &Simple_frame_buffer::create;
    m["Accumulating"] = &Accumulating_frame_buffer::create;
    m["Reweighting"] = &Reweighting_frame_buffer::create;
    m["Distance_weighted"] = &Distance_weighted_frame_buffer::create;
    m["Front_stacked"] = &Front_stacked_frame_buffer::create;
    m["Directional_stacked"] = &Directional_stacked_frame_buffer::create;
    return m;
}

std::map<std::string, Frame_buffer_create_function> Cube_raster_buffer::frame_buffer_create_map = create_fb_create_map();

const vector3d_t Cube_raster_buffer::cube_plane_normals[] =
{
    yafaray::Cube_raster_buffer::Point( 1.0f,  0.0f,  0.0f),
    yafaray::Cube_raster_buffer::Point(-1.0f,  0.0f,  0.0f),
    yafaray::Cube_raster_buffer::Point( 0.0f,  1.0f,  0.0f),
    yafaray::Cube_raster_buffer::Point( 0.0f, -1.0f,  0.0f),
    yafaray::Cube_raster_buffer::Point( 0.0f,  0.0f,  1.0f),
    yafaray::Cube_raster_buffer::Point( 0.0f,  0.0f, -1.0f)
};



std::vector< std::vector<Cube_cell> > create_cube_cells_luts()
{
    int const max_resolution = 64;

    std::vector< std::vector<Cube_cell> > cube_cells_luts(max_resolution + 1);

    // only actually set up evry even resolution size, other should not be used anyway
    for (unsigned int i = 2; i <= cube_cells_luts.size(); i += 2)
    {
        int const resolution = i;
        int const resolution_2 = resolution / 2;

        Cube_raster_buffer crb;
        crb.setup_simple(resolution);

        std::vector<Cube_cell> cells(6 * resolution * resolution);

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            Cube_cell c;
            c.plane = plane_index;

            for (int u = -resolution_2; u < resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -resolution_2; v < resolution_2; ++v)
                {
                    c.pos[1] = v;

                    cells[crb.get_serial_index(c)] = c;
                }
            }
        }

        cube_cells_luts[resolution] = cells;
    }

    return cube_cells_luts;
}

std::vector< std::vector<Cube_cell> > Cube_raster_buffer::cube_cells_luts = create_cube_cells_luts();



std::vector< std::vector<Point> > create_cell_centers_luts()
{
    int const max_resolution = 64;

    std::vector< std::vector<Point> > cell_centers_luts(max_resolution + 1);

    // only actually set up evry even resolution size, other should not be used anyway
    for (unsigned int i = 2; i <= cell_centers_luts.size(); i += 2)
    {
        int const resolution = i;
        int const resolution_2 = resolution / 2;

        Cube_raster_buffer crb;
        crb.setup_simple(resolution);

        std::vector<Point> cell_centers(6 * resolution * resolution);

        for (int plane_index = 0; plane_index < 6; ++plane_index)
        {
            Cube_cell c;
            c.plane = plane_index;

            for (int u = -resolution_2; u < resolution_2; ++u)
            {
                c.pos[0] = u;

                for (int v = -resolution_2; v < resolution_2; ++v)
                {
                    c.pos[1] = v;

                    cell_centers[crb.get_serial_index(c)] = crb.calc_cell_center(c);
                }
            }
        }

        cell_centers_luts[resolution] = cell_centers;
    }

    return cell_centers_luts;
}

std::vector< std::vector<Point> > Cube_raster_buffer::cell_centers_luts = create_cell_centers_luts();


std::vector< std::vector<float> > create_solid_angle_luts()
{
    int const max_resolution = 64;

    std::vector< std::vector<float> > solid_angle_luts(max_resolution + 1);

    // only actually set up evry even resolution size, other should not be used anyway
    for (unsigned int i = 2; i <= solid_angle_luts.size(); i += 2)
    {
        Cube_raster_buffer crb;
        crb.setup_simple(i);

        std::vector<Cube_cell> const& cells = crb.get_cube_cells();

        solid_angle_luts[i].resize(6 * i * i);

        for (unsigned int j = 0; j < cells.size(); ++j)
        {
            int const serial_index = crb.get_serial_index(cells[j]);

            solid_angle_luts[i][serial_index] = crb.calc_solid_angle(cells[j]);
        }
    }

    return solid_angle_luts;
}

std::vector< std::vector<float> > Cube_raster_buffer::solid_angle_luts = create_solid_angle_luts();



/*
  major axis
  direction     target                             sc     tc    ma
  ----------    -------------------------------    ---    ---   ---
   +rx          TEXTURE_CUBE_MAP_POSITIVE_X_ARB    -rz    -ry   rx
   -rx          TEXTURE_CUBE_MAP_NEGATIVE_X_ARB    +rz    -ry   rx
   +ry          TEXTURE_CUBE_MAP_POSITIVE_Y_ARB    +rx    +rz   ry
   -ry          TEXTURE_CUBE_MAP_NEGATIVE_Y_ARB    +rx    -rz   ry
   +rz          TEXTURE_CUBE_MAP_POSITIVE_Z_ARB    +rx    -ry   rz
   -rz          TEXTURE_CUBE_MAP_NEGATIVE_Z_ARB    -rx    -ry   rz

 Using the sc, tc, and ma determined by the major axis direction as
 specified in the table above, an updated (s,t) is calculated as
 follows

    s   =   ( sc/|ma| + 1 ) / 2
    t   =   ( tc/|ma| + 1 ) / 2
*/

int const Cube_raster_buffer::corresponding_axis_0[] = { 5 /*-rz*/, 4 /*+rz*/, 0 /*+rx*/, 0 /*+rx*/, 0 /*+rx*/, 1 /*-rx*/};
int const Cube_raster_buffer::corresponding_axis_1[] = { 3 /*-ry*/, 3 /*-ry*/, 4 /*+rz*/, 5 /*-rz*/, 3 /*-ry*/, 3 /*-ry*/};


// p0, p1: edge to test
// v0: point on the plane, n: plane normal
// alpha: return value, intersection point: p0 + alpha * (p1 - p0)
bool linePlaneIntersection(Cube_raster_buffer::Point const& p0, Cube_raster_buffer::Point const& p1, Cube_raster_buffer::Point const& v0, Cube_raster_buffer::Point const& n, float & alpha)
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
bool linePlaneIntersection(Cube_raster_buffer::Point const& dir, Cube_raster_buffer::Point const& n, float & alpha)
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
bool sphere_line_intersection(Cube_raster_buffer::Point const& sphere_center, float const radius, Cube_raster_buffer::Point const& line_dir)
{
    float const c_l = sphere_center * line_dir;

    if (c_l <= 0.0f) return false;

    return (c_l * c_l) - sphere_center.lengthSqr() + (radius * radius) >= 0.0f;
}


float round(float const d)
{
    return std::floor(d + 0.5f);
}

float next_raster_in_positive_direction(float u)
{
    return std::floor(u + 1.0f);
}

color_t interpolate(color_t const& i1 , color_t const& i2 , float const offset)
{
    assert(offset >= 0.0f && offset <= 1.0f);
    return i1 * (1.0f - offset) + i2 * offset;
}

int get_longest_axis(Cube_raster_buffer::Point const& dir)
{
    int longest_axis = 0;
    if (std::abs(dir[1]) > std::abs(dir[0]) && std::abs(dir[1]) > std::abs(dir[2])) longest_axis = 1;
    else if (std::abs(dir[2]) > std::abs(dir[0])) longest_axis = 2;
    return longest_axis;
}


float pyramid_solid_angle(std::vector<vector3d_t> const& corners)
{
    int orders[][3] = { {0, 1, 2}, {0, 2, 3} };

    float omega = 0.0f;

    for (int i = 0; i < 2; ++i)
    {
        int const* order = orders[i];

        vector3d_t const& a = corners[order[0]];
        vector3d_t const& b = corners[order[1]];
        vector3d_t const& c = corners[order[2]];

        float const determ = std::abs(a * (b ^ c));

        float const al = a.length();
        float const bl = b.length();
        float const cl = c.length();

        float const div = al * bl * cl + (a * b) * cl + (a * c) * bl + (b * c) * al;
        float at = std::atan2(determ, div);

        if (at < 0) at += M_PI; // If det > 0 and div < 0 arctan2 returns < 0, so add pi.
        omega += 2.0f * at;
    }

    return omega;
}



Cube_raster_buffer::Cube_raster_buffer() :
    _get_color_done(false)
{
    for (int i = 0; i < 6; ++i)
    {
        buffers[i] = NULL;
    }
}

Cube_raster_buffer & Cube_raster_buffer::operator= (Cube_raster_buffer const& cbr)
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
        add_point_function_map = cbr.add_point_function_map;
        _get_color_done        = cbr._get_color_done;
    }

    return *this;
}

Cube_raster_buffer::~Cube_raster_buffer()
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

void Cube_raster_buffer::clear()
{
    for (int i = 0; i < 6; ++i)
    {
        buffers[i]->clear();
    }

    // _cell_centers.clear();
    // _cube_cells.clear();

    add_point_function_map.clear();
}

void Cube_raster_buffer::setup_surfel_buffer(int resolution)
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
        buffers[i] = new Simple_single_value_frame_buffer(_resolution, i);
    }
}

void Cube_raster_buffer::setup_simple(int resolution)
{
    setup(Simple, resolution, Single_pixel, Single_pixel, Single_pixel);
}

void Cube_raster_buffer::setup(Type fb_type,
                               int resolution,
                               Splat_type const node_splat_type,
                               Splat_type const surfel_far_splat_type,
                               Splat_type const surfel_near_splat_type)
{
    _resolution = resolution;
    _resolution_2 = _resolution / 2;

    std::map<Splat_type, Add_point_function_ptr> splat_type_to_function_map;
    splat_type_to_function_map[Single_pixel]           = &Cube_raster_buffer::add_point_single_pixel;
    splat_type_to_function_map[Disc_tracing]           = &Cube_raster_buffer::add_point_disc_tracing;
    splat_type_to_function_map[AA_square]              = &Cube_raster_buffer::add_point_aa_square;
    splat_type_to_function_map[SA_tracing]             = &Cube_raster_buffer::add_point_solid_angle_rays;
    splat_type_to_function_map[Stocastic_tracing]      = &Cube_raster_buffer::add_point_stochastic_disc_tracing;
    splat_type_to_function_map[Stocastic_node_tracing] = &Cube_raster_buffer::add_point_stochastic_node_tracing;

    add_point_function_map.resize(3);
    add_point_function_map[Gi_point_info::Node]        = splat_type_to_function_map[node_splat_type];
    add_point_function_map[Gi_point_info::Far_surfel]  = splat_type_to_function_map[surfel_far_splat_type];
    add_point_function_map[Gi_point_info::Near_surfel] = splat_type_to_function_map[surfel_near_splat_type];

    for (int i = 0; i < 6; ++i)
    {
        // buffers[i] = frame_buffer_create_map[() // FIXME: missing id string, need some restructuring

        if (fb_type == Simple)
        {
            // buffers[i] = new Simple_frame_buffer(_resolution, i);
            buffers[i] = new Simple_frame_buffer_without_queue(_resolution, i);
        }
        else if (fb_type == Accumulating)
        {
            buffers[i] = new Accumulating_frame_buffer_without_queue(_resolution, i);
        }
        else if (fb_type == Reweighting)
        {
            buffers[i] = new Reweighting_frame_buffer(_resolution, i);
        }
        else if (fb_type == Distance_weighted)
        {
            buffers[i] = new Grouping_frame_buffer_without_queue(_resolution, i);
//            buffers[i] = new Forward_distance_weighted_frame_buffer(_resolution, i);
        }
        else if (fb_type == Front_stacked)
        {
            buffers[i] = new Front_stacked_frame_buffer(_resolution, i);
        }
        else if (fb_type == Directional_stacked)
        {
            buffers[i] = new Directional_stacked_frame_buffer(_resolution, i);
        }
        else
        {
            std::cout << "Unknown FB type ..." << std::endl;
            assert(false);
        }
    }
}

int Cube_raster_buffer::get_resolution() const
{
    return _resolution;
}


Cube_cell Cube_raster_buffer::get_cell(Point const& dir, bool * ok) const
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

    int const axis_0 = corresponding_axis_0[plane_index] / 2;
    int const axis_1 = corresponding_axis_1[plane_index] / 2;

    int const sign_axis_0 = (corresponding_axis_0[plane_index] % 2 == 0) ? 1 : -1;
    int const sign_axis_1 = (corresponding_axis_1[plane_index] % 2 == 0) ? 1 : -1;

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


inline int Cube_raster_buffer::get_serial_index(Cube_cell const& c) const
{
    int serial =
            (c.pos[1] + _resolution_2) +
            (c.pos[0] + _resolution_2) * _resolution +
            c.plane * _resolution * _resolution;

    return serial;
}

Cube_raster_buffer::Point const& Cube_raster_buffer::get_cell_center(Cube_cell const& c) const
{
    int serial = get_serial_index(c);

    return get_cell_center(serial);
}

inline Cube_raster_buffer::Point const& Cube_raster_buffer::get_cell_center(int const serial_index) const
{
    return cell_centers_luts[_resolution][serial_index];
}

Cube_raster_buffer::Point Cube_raster_buffer::get_cell_direction(int const serial_index) const
{
    Point dir = get_cell_center(serial_index);
    dir.normalize();
    return dir;
}

Cube_raster_buffer::Point Cube_raster_buffer::get_cell_direction(Cube_cell const& c) const
{
    Point dir = get_cell_center(c);
    dir.normalize();
    return dir;
}

std::vector<Cube_cell> const& Cube_raster_buffer::get_cube_cells() const
{
    //return _cube_cells;
    return cube_cells_luts[_resolution];
}

std::vector<Cube_raster_buffer::Point> Cube_raster_buffer::get_cell_corners(Cube_cell const& c) const
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



void Cube_raster_buffer::add_point_single_pixel(Gi_point_info const& point_info, GiPoint const* node)
{
    Point dir = point_info.direction;
    Color const& color = point_info.color;
    float const depth = point_info.depth;
    float const radius = point_info.radius;
    dir.normalize();

    bool ok;
    Cube_cell c = get_cell(dir, &ok);
    if (!ok) return;

    buffers[c.plane]->add_point(c.pos[0], c.pos[1], color, 1.0f, depth, radius, dir, node);
}

// raytrace a disc through every pixel, assumes receiving point is in the origin and the disc's center relative to it
void Cube_raster_buffer::add_point_disc_tracing(Gi_point_info const& point_info, GiPoint const* node)
{
    // Color const& color       = point_info.color;
    Point const& disc_normal = point_info.disc_normal;
    // Point const& dir         = point_info.direction;
    Point disc_center        = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
    float const disc_radius  = point_info.radius;
    // float const depth = point_info.depth;

    Cube_cell c;

    float disc_radius_sqr = disc_radius * disc_radius;
    float cell_radius_sqr = 1.0f / float(_resolution * _resolution);

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

                    Point const relative_hitpoint = hit_point - disc_center;
                    float dist_sqr = relative_hitpoint.lengthSqr();

                    if (dist_sqr < disc_radius_sqr)
                    {
                        //float const relative_dist = dist_sqr / disc_radius_sqr;
                        //float const gauss_weight = std::exp(-3.0f * relative_dist);
                        float const gauss_weight = 1.0f;

                        if (gauss_weight > 0.05f)
                        {
                            // float const filling_degree = std::min(1.0f, point_info.solid_angle / _pixel_solid_angle);
                            float const filling_degree = 1.0f;

                            Point dir = point_info.receiver_position - (relative_hitpoint + point_info.position);
                            dir.normalize();
                            Color color = point_info.spherical_function->get_color(dir) * gauss_weight;

                            buffers[plane_index]->add_point(u, v, color, filling_degree, alpha, disc_radius, dir, node);
                        }
                    }
                }

                /*

                    float alpha = -1.0f;
                    linePlaneIntersection(Point(0.0f), disc_center, cell_center, -cube_plane_normals[plane_index], alpha);

                    if (alpha > 0.0f)
                    {
                        Point hit_point = alpha * disc_center;

                        // float dist_sqr = (hit_point - cell_center).lengthSqr();

                        Point vec_hp_cell_center_abs = (hit_point - cell_center);
                        vec_hp_cell_center_abs.abs();

                        if (vec_hp_cell_center_abs.x < _resolution_2 && vec_hp_cell_center_abs.y < _resolution_2 && vec_hp_cell_center_abs.z < _resolution_2)
                        // if (dist_sqr < cell_radius_sqr)
                        {
                            float const filling_degree = std::min(1.0f, point_info.solid_angle / _pixel_solid_angle);

                            buffers[plane_index]->add_point(u, v, color, filling_degree, alpha, disc_radius, dir, node);
                        }
                    }
                    */
            }
        }
    }
}


void Cube_raster_buffer::add_point_stochastic_disc_tracing(Gi_point_info const& point_info, GiPoint const* node)
{
    Color const& color       = point_info.color;
    Point const& disc_normal = point_info.disc_normal;
    Point const& dir         = point_info.direction;
    Point const disc_center  = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
    float const disc_radius  = point_info.radius;
    float const disc_solid_angle = point_info.solid_angle;

    std::map<Cube_cell, int> coverage_map;

    int const num_samples = _resolution * _resolution * 0.5f;
    // int const num_samples = 20;

    random_t my_random;

    vector3d_t NU, NV;
    createCS(disc_normal, NU, NV);

//    for (int i = 0; i < num_samples; ++i)
//    {
//        vector3d_t sample_uv = SampleDisc(my_random(), my_random());

//        vector3d_t sample_point = disc_center + NU * sample_uv.x * disc_radius + NV * sample_uv.y * disc_radius;

//        vector3d_t sample_dir = sample_point;
//        sample_dir.normalize();

//        Cube_cell hit_cell = get_cell(sample_dir);

//        vector3d_t const& cell_center = get_cell_center(hit_cell);
//        int const plane_index = hit_cell.plane;

//        float alpha = -1.0f;
//        linePlaneIntersection(Point(0.0f), sample_point, cell_center, -cube_plane_normals[plane_index], alpha);

//        if (alpha > 0.0f)
//        {
//            Point hit_point = alpha * disc_center;

//            // float dist_sqr = (hit_point - cell_center).lengthSqr();

//            Point vec_hp_cell_center_abs = (hit_point - cell_center);
//            vec_hp_cell_center_abs.abs();

//            if (vec_hp_cell_center_abs.x < _resolution_2 && vec_hp_cell_center_abs.y < _resolution_2 && vec_hp_cell_center_abs.z < _resolution_2)
//            {
//                coverage_map[hit_cell] += 1;
//            }
//        }
//    }

    for (int i = 0; i < num_samples; ++i)
    {
        vector3d_t sample_uv = SampleDisc(my_random(), my_random());

        vector3d_t sample_point = disc_center + NU * sample_uv.x * disc_radius + NV * sample_uv.y * disc_radius;

        vector3d_t sample_dir = sample_point;
        sample_dir.normalize();

        Cube_cell hit_cell = get_cell(sample_dir);
        coverage_map[hit_cell] += 1;
    }

    for (std::map<Cube_cell, int>::const_iterator iter = coverage_map.begin(); iter != coverage_map.end(); ++iter)
    {
        Cube_cell const& c = iter->first;

        float const cell_solid_angle = get_solid_angle(c);
        float const sample_ratio = (iter->second / float(num_samples));
        float const solid_angle_ratio = disc_solid_angle / cell_solid_angle;
        float const filling_degree = std::min(1.0f, sample_ratio * solid_angle_ratio);

        buffers[c.plane]->add_point(c.pos[0], c.pos[1], color, filling_degree, point_info.depth, disc_radius, dir, node);
    }
}



void Cube_raster_buffer::add_point_stochastic_node_tracing(Gi_point_info const& point_info, GiPoint const* node)
{
    Color const& color            = point_info.color;
    Point const& node_normal      = point_info.direction;
    Point const& dir              = point_info.direction;
    Point const  node_center      = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
    float const  node_radius      = point_info.radius;
    float const  node_solid_angle = point_info.solid_angle;

    std::map<Cube_cell, int> coverage_map;

    int const num_samples = _resolution * _resolution * 0.5f;
    // int const num_samples = 20;

    random_t my_random;

    vector3d_t NU, NV;
    createCS(node_normal, NU, NV);

    for (int i = 0; i < num_samples; ++i)
    {
        vector3d_t sample_uv = SampleDisc(my_random(), my_random());

//        float s1, s2;
//        hammersley_2(i, num_samples, s1, s2);
//        vector3d_t sample_uv = SampleDisc(s1, s2);

        vector3d_t sample_point = node_center + NU * sample_uv.x * node_radius + NV * sample_uv.y * node_radius;

        vector3d_t sample_dir = sample_point;
        sample_dir.normalize();

        Cube_cell hit_cell = get_cell(sample_dir);
        coverage_map[hit_cell] += 1;
    }

    for (std::map<Cube_cell, int>::const_iterator iter = coverage_map.begin(); iter != coverage_map.end(); ++iter)
    {
        Cube_cell const& c = iter->first;

        float const cell_solid_angle = get_solid_angle(c);
        float const sample_ratio = (iter->second / float(num_samples));
        float const solid_angle_ratio = node_solid_angle / cell_solid_angle;
        float const filling_degree = std::min(1.0f, sample_ratio * solid_angle_ratio);

        buffers[c.plane]->add_point(c.pos[0], c.pos[1], color, filling_degree, point_info.depth, node_radius, dir, node);
    }
}



void Cube_raster_buffer::add_point_aa_square(Gi_point_info const& point_info, GiPoint const* node)
{
    float const depth           = point_info.depth;
    color_t const& color        = point_info.color;
    float const solid_angle     = point_info.solid_angle;
    vector3d_t const& direction = point_info.direction;
    float const radius          = point_info.radius;

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

    float const square_area = solid_angle * 4.0f;
    float const square_width_2 = std::sqrt(square_area) / 2.0f;

    // find intersection point on plane
    for (int i = 0; i < 3; ++i)
    {
        int const normal_axis = i;

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

        //            int const axis_0 = (normal_axis + 1) % 3;
        //            int const axis_1 = (normal_axis + 2) % 3;

        int const axis_0 = corresponding_axis_0[plane_index] / 2;
        int const axis_1 = corresponding_axis_1[plane_index] / 2;

        int const sign_axis_0 = (corresponding_axis_0[plane_index] % 2 == 0) ? 1 : -1;
        int const sign_axis_1 = (corresponding_axis_1[plane_index] % 2 == 0) ? 1 : -1;

        float u_minus = (hit_point[axis_0] - square_width_2) * _resolution_2 * sign_axis_0;
        float u_plus  = (hit_point[axis_0] + square_width_2) * _resolution_2 * sign_axis_0;

        float v_minus = (hit_point[axis_1] - square_width_2) * _resolution_2 * sign_axis_1;
        float v_plus  = (hit_point[axis_1] + square_width_2) * _resolution_2 * sign_axis_1;

        float const u_min = u_minus;
        float const u_max = u_plus;

        float const v_min = v_minus;
        float const v_max = v_plus;

        float const u_center = (u_max + u_min) / 2.0f;
        float const v_center = (v_max + v_min) / 2.0f;

        u_minus = into_range(-float(_resolution_2), float(_resolution_2), u_minus);
        v_minus = into_range(-float(_resolution_2), float(_resolution_2), v_minus);

        u_plus = into_range(-float(_resolution_2), float(_resolution_2), u_plus);
        v_plus = into_range(-float(_resolution_2), float(_resolution_2), v_plus);

        if (u_minus > u_plus) std::swap(u_minus, u_plus);
        if (v_minus > v_plus) std::swap(v_minus, v_plus);

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

                float const sigma = (square_width_2 * _resolution_2) / 2.0f;
                float const dist_from_center = std::sqrt((cell_u + 0.5f - u_center) * (cell_u + 0.5f - u_center) + (cell_v + 0.5f - v_center) * (cell_v + 0.5f - v_center));
                float const x = dist_from_center;
                // float const gauss = std::min(1.0, 1.0f / std::sqrt(M_2PI * sigma * sigma) * std::exp(-(x*x) / (2.0f * sigma * sigma)));
                float const gauss = std::min(1.0f, std::exp(-(x*x) / (2.0f * sigma * sigma)));
                fill_ratio *= gauss;

                buffers[plane_index]->add_point(cell_u, cell_v, color, fill_ratio, depth, radius, dir, node);

                v = next_v;
            }

            u = next_u;
        }
    }
}


void Cube_raster_buffer::add_point_solid_angle_rays(Gi_point_info const& point_info, GiPoint const* node)
{
    Color const& color       = point_info.color;
    Point const& disc_normal = (point_info.receiver_position - point_info.position).normalize();
    Point const& dir         = point_info.direction;
    Point disc_center        = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
    float const disc_radius  = std::sqrt(point_info.solid_angle * point_info.depth * point_info.depth / M_PI);
    float const radius       = point_info.radius;
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
                        buffers[plane_index]->add_point(u, v, color, 1.0f, alpha, radius, dir, node);
                    }
                }
            }
        }
    }
}





void Cube_raster_buffer::add_point(Gi_point_info const& point_info, GiPoint const* gi_point)
{
    (this->*add_point_function_map[point_info.type])(point_info, gi_point);
}

void Cube_raster_buffer::set_color(Cube_cell const& c, Color const& color)
{
    buffers[c.plane]->set_color(c.pos[0], c.pos[1], color);
}

void Cube_raster_buffer::add_color(Cube_cell const& c, Color const& color)
{
    Color current_color = get_color(c);
    set_color(c, current_color + color);
}

Cube_raster_buffer::Color Cube_raster_buffer::get_color(Cube_cell const& c, Debug_info * debug_info) const
{
    return buffers[c.plane]->get_color(c.pos[0], c.pos[1], debug_info);
}

Cube_raster_buffer::Color Cube_raster_buffer::get_color(int const plane, int const x, int const y, Debug_info * debug_info) const
{
    return buffers[plane]->get_color(x, y, debug_info);
}

Cube_raster_buffer::Color Cube_raster_buffer::get_color_interpolated(Cube_cell const& c, Debug_info * debug_info) const
{
    // return buffers[c.plane]->get_color_interpolated(c.pos[0] + c.offset[0], c.pos[1] + c.offset[1], debug_info);

    Cube_cell cell = c;

    Color colors[4];

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
            colors[u + v * 2] = get_color(new_cell);

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

    Color const c1_2 = interpolate(colors[0], colors[1], offset_x);
    Color const c3_4 = interpolate(colors[2], colors[3], offset_x);

    Color const result = interpolate(c1_2, c3_4, offset_y);

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

Abstract_frame_buffer * Cube_raster_buffer::get_component(int const index) const
{
    return buffers[index];
}

std::vector<float> Cube_raster_buffer::component_to_vector(int const component_index, bool const use_first_component_only) const
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

            Color const& color = get_color(c);

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

void Cube_raster_buffer::from_component_vector(int const component_index, std::vector<float> const& v, bool const use_first_component_only)
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

            Color color;

            if (use_first_component_only)
            {
                color = Color(v[i + 0], v[i + 0], v[i + 0]);
                i += 1;
            }
            else
            {
                color = Color(v[i + 0], v[i + 1], v[i + 2]);
                i += 3;
            }

            set_color(c, color);
        }
    }
}


float Cube_raster_buffer::calc_solid_angle(Cube_cell const& c) const
{
    std::vector<vector3d_t> corners = get_cell_corners(c);
    float const solid_angle = pyramid_solid_angle(corners);
    return solid_angle;
}

float Cube_raster_buffer::get_solid_angle(Cube_cell const& c) const
{
    int const serial_index = get_serial_index(c);
    return solid_angle_luts[_resolution][serial_index];
}

float Cube_raster_buffer::get_solid_angle(int const serial_index) const
{
    return solid_angle_luts[_resolution][serial_index];
}

float Cube_raster_buffer::get_solid_angle(Point const& dir) const
{
    Cube_cell c = get_cell(dir);
    return get_solid_angle(c);
}


Cube_raster_buffer::Color Cube_raster_buffer::get_diffuse(Point const& normal) const
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

                Color color = get_color(c);

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

    // std::cout << "get_diffuse(): " << debug_weight_sum << std::endl;

    _total_energy = total.energy();

    return diffuse;
}

float Cube_raster_buffer::get_total_energy()
{
    return _total_energy;
}

float Cube_raster_buffer::get_non_zero_area()
{
    return _non_zero_area;
}


void Cube_raster_buffer::get_debug_info(Debug_info & debug_info)
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

Cube_raster_buffer::Point Cube_raster_buffer::calc_cell_center(Cube_cell const& c) const
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

    int const axis_0 = corresponding_axis_0[c.plane] / 2;
    int const axis_1 = corresponding_axis_1[c.plane] / 2;

    int const sign_axis_0 = (corresponding_axis_0[c.plane] % 2 == 0) ? 1 : -1;
    int const sign_axis_1 = (corresponding_axis_1[c.plane] % 2 == 0) ? 1 : -1;

    result[axis_0] = ((c.pos[0] + 0.5f) * sign_axis_0) / float(_resolution_2);
    result[axis_1] = ((c.pos[1] + 0.5f) * sign_axis_1) / float(_resolution_2);

    return result;
}

Cube_raster_buffer::Point Cube_raster_buffer::calc_cell_position_with_offset(Cube_cell const& c) const
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

    int const axis_0 = corresponding_axis_0[c.plane] / 2;
    int const axis_1 = corresponding_axis_1[c.plane] / 2;

    int const sign_axis_0 = (corresponding_axis_0[c.plane] % 2 == 0) ? 1 : -1;
    int const sign_axis_1 = (corresponding_axis_1[c.plane] % 2 == 0) ? 1 : -1;

    result[axis_0] = ((c.pos[0] + 0.5f + c.offset[0]) * sign_axis_0) / float(_resolution_2);
    result[axis_1] = ((c.pos[1] + 0.5f + c.offset[1]) * sign_axis_1) / float(_resolution_2);

    return result;
}

int Cube_raster_buffer::get_corresponding_axis(int const axis, int const plane)
{
    assert(axis == 0 || axis == 1);

    if (axis == 0)
    {
        return corresponding_axis_0[plane];
    }

    return corresponding_axis_1[plane];
}

__END_YAFRAY
