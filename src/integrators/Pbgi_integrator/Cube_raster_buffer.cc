#include <map>

#include <integrators/Pbgi_integrator/CubeRasterBuffer.h>
#include <utilities/math_utils.h>
#include <integrators/Pbgi_integrator/Frame_buffer.h>
#include <utilities/interpolation.h>
#include <core_api/material.h>


__BEGIN_YAFRAY

std::map<std::string, Splat_cube_raster_buffer::Splat_type> create_splat_type_map()
{
    std::map<std::string, Splat_cube_raster_buffer::Splat_type> m;
    m["Single_pixel"] = Splat_cube_raster_buffer::Single_pixel;
    m["Disc_splat_strategy"] = Splat_cube_raster_buffer::Disc_tracing;
    m["Square_splat_strategy"] = Splat_cube_raster_buffer::AA_square;
    m["Gaussian_splat_strategy"] = Splat_cube_raster_buffer::Gaussian_splat;
    m["SA_tracing"] = Splat_cube_raster_buffer::SA_tracing;
    m["Stocastic_tracing"] = Splat_cube_raster_buffer::Stocastic_tracing;
    m["Stocastic_node_tracing"] = Splat_cube_raster_buffer::Stocastic_node_tracing;
    return m;
}

std::map<std::string, Splat_cube_raster_buffer::Splat_type> Splat_cube_raster_buffer::enum_splat_type_map = create_splat_type_map();

std::map<std::string, Splat_cube_raster_buffer::Buffer_type> create_type_map()
{
    std::map<std::string, Splat_cube_raster_buffer::Buffer_type> m;
    m["Simple"] = Splat_cube_raster_buffer::Simple;
    m["Accumulating_frame_buffer_without_queue"] = Splat_cube_raster_buffer::Accumulating;
    //m["Reweighting"] = Splat_cube_raster_buffer::Reweighting;
    m["Distance_weighted"] = Splat_cube_raster_buffer::Distance_weighted;
    //m["Front_stacked"] = Splat_cube_raster_buffer::Front_stacked;
    //m["Directional_stacked"] = Splat_cube_raster_buffer::Directional_stacked;
    m["Parameter"] = Splat_cube_raster_buffer::Parameter;
    return m;
}

std::map<std::string, Splat_cube_raster_buffer::Buffer_type> Splat_cube_raster_buffer::enum_type_map = create_type_map();


/*
std::map<std::string, Frame_buffer_create_function<Data> > create_fb_create_map()
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
*/

const vector3d_t Cube_static_data::cube_plane_normals[] =
{
    vector3d_t( 1.0f,  0.0f,  0.0f),
    vector3d_t(-1.0f,  0.0f,  0.0f),
    vector3d_t( 0.0f,  1.0f,  0.0f),
    vector3d_t( 0.0f, -1.0f,  0.0f),
    vector3d_t( 0.0f,  0.0f,  1.0f),
    vector3d_t( 0.0f,  0.0f, -1.0f)
};




std::vector< std::vector<Cube_cell> > create_cube_cells_luts()
{
    std::cout << "create_cube_cells_luts()" << std::endl;

    int const max_resolution = 64;

    std::vector< std::vector<Cube_cell> > cube_cells_luts(max_resolution + 1);

    // only actually set up evry even resolution size, other should not be used anyway
    for (unsigned int i = 2; i <= cube_cells_luts.size(); i += 2)
    {
        int const resolution = i;
        int const resolution_2 = resolution / 2;

        Cube_raster_buffer<float> crb;
        crb.setup_surfel_buffer(resolution);

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

std::vector< std::vector<Cube_cell> > Cube_static_data::cube_cells_luts = create_cube_cells_luts();



std::vector< std::vector<vector3d_t> > create_cell_centers_luts()
{
    std::cout << "create_cell_centers_luts()" << std::endl;

    int const max_resolution = 64;

    std::vector< std::vector<vector3d_t> > cell_centers_luts(max_resolution + 1);

    // only actually set up evry even resolution size, other should not be used anyway
    for (unsigned int i = 2; i <= cell_centers_luts.size(); i += 2)
    {
        int const resolution = i;
        int const resolution_2 = resolution / 2;

        Cube_raster_buffer<int> crb;
        crb.setup_surfel_buffer(resolution);

        std::vector<vector3d_t> cell_centers(6 * resolution * resolution);

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

std::vector< std::vector<vector3d_t> > Cube_static_data::cell_centers_luts = create_cell_centers_luts();


std::vector< std::vector<float> > create_solid_angle_luts()
{
    std::cout << "create_solid_angle_luts()" << std::endl;

    int const max_resolution = 64;

    std::vector< std::vector<float> > solid_angle_luts(max_resolution + 1);

    // only actually set up evry even resolution size, other should not be used anyway
    for (unsigned int i = 2; i <= solid_angle_luts.size(); i += 2)
    {
        Cube_raster_buffer<int> crb;
        crb.setup_surfel_buffer(i);

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

std::vector< std::vector<float> > Cube_static_data::solid_angle_luts = create_solid_angle_luts();



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

int const Cube_static_data::corresponding_axis_0[] = { 5 /*-rz*/, 4 /*+rz*/, 0 /*+rx*/, 0 /*+rx*/, 0 /*+rx*/, 1 /*-rx*/};
int const Cube_static_data::corresponding_axis_1[] = { 3 /*-ry*/, 3 /*-ry*/, 4 /*+rz*/, 5 /*-rz*/, 3 /*-ry*/, 3 /*-ry*/};


// p0, p1: edge to test
// v0: point on the plane, n: plane normal
// alpha: return value, intersection point: p0 + alpha * (p1 - p0)
bool linePlaneIntersection(vector3d_t const& p0, vector3d_t const& p1, vector3d_t const& v0, vector3d_t const& n, float & alpha)
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
bool linePlaneIntersection(vector3d_t const& dir, vector3d_t const& n, float & alpha)
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
bool sphere_line_intersection(vector3d_t const& sphere_center, float const radius, vector3d_t const& line_dir)
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


int get_longest_axis(vector3d_t const& dir)
{
    int longest_axis = 0;
    if (std::abs(dir[1]) > std::abs(dir[0]) && std::abs(dir[1]) > std::abs(dir[2])) longest_axis = 1;
    else if (std::abs(dir[2]) > std::abs(dir[0])) longest_axis = 2;
    return longest_axis;
}



float pyramid_solid_angle(std::vector<vector3d_t> const& corners)
{
    assert(corners.size() == 4);

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



void Splat_cube_raster_buffer::setup_empty()
{
    _resolution = 2;
    _resolution_2 = _resolution / 2;

    for (int i = 0; i < 6; ++i)
    {
        // buffers[i] = new Simple_frame_buffer_without_queue<color_t>(_resolution, i);
        buffers[i] = boost::shared_ptr<Abstract_frame_buffer<color_t> >(new Simple_frame_buffer_without_queue<color_t>(_resolution));
    }
}




void Splat_cube_raster_buffer::setup(int const resolution)
{
    _resolution = resolution;
    _resolution_2 = _resolution / 2;

    for (int i = 0; i < 6; ++i)
    {
        buffers[i] = boost::shared_ptr<Abstract_frame_buffer<color_t> >(new Accumulating_frame_buffer_without_queue<color_t>(resolution));
    }

    add_point_function_map.resize(3);
    add_point_function_map[Gi_point_info::Node]        = new Square_splat_strategy;
    add_point_function_map[Gi_point_info::Far_surfel]  = new Square_splat_strategy;
    add_point_function_map[Gi_point_info::Near_surfel] = new Disc_splat_strategy;
}



void Splat_cube_raster_buffer::add_point(Gi_point_info const& point_info)
{
    // (this->*add_point_function_map[point_info.type])(point_info, gi_point);
    add_point_function_map[point_info.type]->splat(*this, point_info);
}


void Splat_cube_raster_buffer::add_background(background_t * background)
{
    Cube_cell c;

    Gi_point_info point_info;

    for (int plane_index = 0; plane_index < 6; ++plane_index)
    {
        c.plane = plane_index;

        for (int u = -_resolution_2; u < _resolution_2; ++u)
        {
            c.pos[0] = u;

            for (int v = -_resolution_2; v < _resolution_2; ++v)
            {
                c.pos[1] = v;

                ray_t ray(point3d_t(0.0f), get_cell_direction(c));
                buffers[plane_index]->add_point(u, v, background->eval(ray), 1.0f, point_info);
            }
        }
    }
}


color_t Splat_cube_raster_buffer::get_brdf_response(renderState_t & state, surfacePoint_t const& receiving_point, vector3d_t const& wo) const
{
    assert(!_get_color_done);
    _get_color_done = true;

    color_t outgoing_color(0.0f);
    color_t total(0.0f);

    _non_zero_area = 0;

    Cube_cell c;

    float debug_weight_sum = 0.0f;
    float weight_sum = 0.0f;

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
                Point const cell_dir = get_cell_direction(serial_index);

                float const cos_sp_normal_cell_dir = cell_dir * receiving_point.N; // assumes the receiver is a lambertian reflector
                assert(cos_sp_normal_cell_dir <= 1.0001f);
                if (cos_sp_normal_cell_dir < 0.001f) continue;

                float const solid_angle = get_solid_angle(serial_index);
                float const weight = solid_angle; // normalization not to 1 but the hemisphere area (2.0f * M_PI) which it already is, FIXME: still not sure about the norm.
                assert(weight >= 0.0f);

                color_t const incoming_color = get_data(c);

                vector3d_t const& light_dir = cell_dir;


//                outgoing_color += pdf * incoming_color * weight * cos_sp_normal_cell_dir;

                color_t const surf_col = receiving_point.material->eval(state, receiving_point, wo, light_dir, BSDF_ALL);
                outgoing_color += surf_col * incoming_color * weight * cos_sp_normal_cell_dir;

                float const pdf = receiving_point.material->pdf(state, receiving_point, wo, light_dir, BSDF_ALL);
                weight_sum += pdf * weight;


//                total += color; // * weight;

//                if (color.energy() > 0.0000001f)
//                {
//                    // std::cout << "non_zero:" << solid_angle << std::endl;
//                    _non_zero_area += solid_angle;
//                }


                debug_weight_sum += weight;

                // outgoing_color += color * cos_sp_normal_cell_dir * weight;
            }
        }
    }

    // not sure about the normalization, the solid angle weights sum up to 2 pi (over the hemisphere)
    // so it seems reasonable to take the factor out again
    // outgoing_color *= 1.0f / (2.0f * M_PI);
    outgoing_color *= 1.0f / weight_sum;

    // std::cout << "get_diffuse(): " << debug_weight_sum << std::endl;

    _total_energy = total.energy();

    return outgoing_color;
}


int Cube_static_data::get_corresponding_axis(int const axis, int const plane)
{
    assert(axis == 0 || axis == 1);

    if (axis == 0)
    {
        return corresponding_axis_0[plane];
    }

    return corresponding_axis_1[plane];
}




Splat_cube_raster_buffer Splat_cube_raster_buffer::blur() const
{
    Splat_cube_raster_buffer tmp_buffer;
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

                Color d(0.0f);

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

                d *= 1.0f / ((blur_size * 2 + 1) * (blur_size * 2 +1));
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



void Square_splat_strategy::splat(Splat_cube_raster_buffer & cube_buffer, Gi_point_info const& point_info)
{
    int const _resolution_2 = cube_buffer.get_resolution() / 2;

    color_t const& color        = point_info.color;
    float const solid_angle     = point_info.solid_angle;
    vector3d_t const& direction = point_info.direction;

    float const square_area = solid_angle;
    float const square_width_2 = std::sqrt(square_area) / 2.0f;


    // find intersection point on plane
    for (int i = 0; i < 3; ++i)
    {
        int const normal_axis = i;

        if (std::abs(direction[normal_axis]) < 1e-10f) continue;

        int plane_index = normal_axis * 2;

        if (direction[normal_axis] < 0.0f) ++plane_index;

        float const cos_hit_dir_normal = direction * Cube_static_data::cube_plane_normals[plane_index];
        assert(cos_hit_dir_normal > 0.0f);

        int const axis_0 = Cube_static_data::corresponding_axis_0[plane_index] / 2;
        int const axis_1 = Cube_static_data::corresponding_axis_1[plane_index] / 2;

        vector3d_t hit_point = direction;
        hit_point[axis_0] /= std::abs(hit_point[normal_axis]);
        hit_point[axis_1] /= std::abs(hit_point[normal_axis]);
        hit_point[normal_axis] = hit_point[normal_axis] > 0.0f ? 1.0f : -1.0f;

        float corrected_width_2[] = { square_width_2, square_width_2 };

        int const sign_axis_0 = (Cube_static_data::corresponding_axis_0[plane_index] % 2 == 0) ? 1 : -1;
        int const sign_axis_1 = (Cube_static_data::corresponding_axis_1[plane_index] % 2 == 0) ? 1 : -1;

        float u_minus = (hit_point[axis_0] - corrected_width_2[0]) * _resolution_2 * sign_axis_0;
        float u_plus  = (hit_point[axis_0] + corrected_width_2[0]) * _resolution_2 * sign_axis_0;

        float v_minus = (hit_point[axis_1] - corrected_width_2[1]) * _resolution_2 * sign_axis_1;
        float v_plus  = (hit_point[axis_1] + corrected_width_2[1]) * _resolution_2 * sign_axis_1;

        u_minus = into_range(-float(_resolution_2), float(_resolution_2), u_minus);
        v_minus = into_range(-float(_resolution_2), float(_resolution_2), v_minus);

        u_plus = into_range(-float(_resolution_2), float(_resolution_2), u_plus);
        v_plus = into_range(-float(_resolution_2), float(_resolution_2), v_plus);

        if (u_minus > u_plus) std::swap(u_minus, u_plus);
        if (v_minus > v_plus) std::swap(v_minus, v_plus);


        float u = u_minus;

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

                float fill_ratio;

                int const cell_u = into_range(-_resolution_2, _resolution_2 - 1, int(std::floor(u)));
                int const cell_v = into_range(-_resolution_2, _resolution_2 - 1, int(std::floor(v)));

                float const covered_pixel_area = (next_u - u) * (next_v - v);
                fill_ratio = std::min(1.0f, covered_pixel_area);

                cube_buffer.get_buffer(plane_index)->add_point(cell_u, cell_v, color, fill_ratio, point_info);

                v = next_v;
            }

            u = next_u;
        }
    }
}

void Disc_splat_strategy::splat(Splat_cube_raster_buffer & cube_buffer, Gi_point_info const& point_info)
{
    int const _resolution_2 = cube_buffer.get_resolution() / 2;

    vector3d_t const& disc_normal = point_info.disc_normal;
    vector3d_t disc_center        = point_info.position - point_info.receiver_position; // -> make receiving point the origin as seen from the disc's center
    float const disc_radius       = point_info.radius;
    color_t const& color          = point_info.color;

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

                vector3d_t cell_center = cube_buffer.get_cell_center(c);

                cell_center.normalize();

                float alpha = -1.0f;
                linePlaneIntersection(vector3d_t(0.0f), cell_center, disc_center, disc_normal, alpha);

                if (alpha > 0.0f)
                {
                    vector3d_t hit_point = alpha * cell_center;

                    vector3d_t const relative_hitpoint = hit_point - disc_center;
                    float dist_sqr = relative_hitpoint.lengthSqr();

                    if (dist_sqr < disc_radius_sqr)
                    {
                        cube_buffer.get_buffer(plane_index)->add_point(u, v, color, 1.0f, point_info);
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

__END_YAFRAY
