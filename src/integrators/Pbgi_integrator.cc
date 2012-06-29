#include <integrators/Pbgi_integrator.h>

#include <sstream>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>


#include <omp.h>

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <core_api/camera.h>
#include <yafraycore/meshtypes.h>
#include <yafraycore/timer.h>
#include <utilities/mcqmc.h>
#include <integrators/Pbgi_integrator/CubeRasterBuffer.h>

#include <utilities/sample_utils.h>

#include <integrators/Pbgi_integrator/Scene_sampler.h>

__BEGIN_YAFRAY



Pbgi_integrator_t::Pbgi_integrator_t(bool transpShad, int shadowDepth, int rayDepth)
{
    type = SURFACE;
    causRadius = 0.25;
    causDepth = 10;
    nCausPhotons = 100000;
    nCausSearch = 100;
    trShad = transpShad;
    usePhotonCaustics = false;
    sDepth = shadowDepth;
    rDepth = rayDepth;
    intpb = 0;
    integratorName = "PB";
    integratorShortName = "PBGI";
}



Cube_raster_buffer< Spherical_harmonics<float> > precalculate_spherical_harmonics(Spherical_harmonics_factory<float> const& spherical_function_area_factory)
{
    std::cout << "precalculate_spherical_harmonics()" << std::endl;

    int const resolution = 32;

    Cube_raster_buffer< Spherical_harmonics<float> > normal_cube;
    normal_cube.setup_surfel_buffer(resolution);

    std::vector<Cube_cell> const& cells = normal_cube.get_cube_cells();

    for (size_t i = 0; i < cells.size(); ++i)
    {
        vector3d_t normal = normal_cube.get_cell_direction(cells[i]);
        Abstract_spherical_function_estimator<float> * area_estimator = new Spherical_function_area_estimator<float>(normal, 1.0f);

        Spherical_harmonics<float> sf = spherical_function_area_factory.create();
        sf.calc_coefficients(area_estimator);

        normal_cube.set_data(cells[i], sf);

        delete area_estimator;
    }

    return normal_cube;
}


void Pbgi_integrator_t::generate_gi_points_data_2(std::vector<MyTree::Tree_node*> const& nodes,
                                             std::tr1::unordered_map<vector3d_t, pbgi_sample_t*, Vector_hash> const& point_to_sampling_point_map,
                                             float const area,
                                             float const radius,
                                             int const treat_depth_as_leaf)
{
    std::cout << "Pbgi_integrator_t::generate_gi_points_data_2(), small tree, nodes: " << nodes.size() << std::endl;

    Gi_point_averager averager(_spherical_function_color_factory, _spherical_function_area_factory);

    renderState_t state;
    unsigned char userdata[USER_DATA_SIZE+7];
    state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
    state.cam = scene->getCamera();

    for (unsigned int node_index = 0; node_index < nodes.size(); ++node_index)
    {
        MyTree::Tree_node* node = nodes[node_index];

        if ((node_index % 1000) == 0)
        {
            std::cout << "." << std::flush;
        }

        if (node->is_leaf())
        {
            MyTree::Tree_leaf_node * leaf_node = node->get_derived<MyTree::Tree_leaf_node>();

            std::vector<Gi_point_surfel> & surfels = leaf_node->get_surfels();
            assert(surfels.size() > 0);

            for (size_t surfel_index = 0; surfel_index < surfels.size(); ++surfel_index)
            {
                Gi_point_surfel & surfel = surfels[surfel_index];

                std::tr1::unordered_map<vector3d_t, pbgi_sample_t*>::const_iterator iter = point_to_sampling_point_map.find(surfel.pos);
                assert(iter != point_to_sampling_point_map.end());

                pbgi_sample_t & sampling_point = *(iter->second);

                point3d_t const& hitPoint = sampling_point.position;
                triangle_t const& tri = *sampling_point.tri_pointer;
                intersectData_t & iData = sampling_point.intersect_data;

                surfacePoint_t sp;

                tri.getSurface(sp, hitPoint, iData);

                assert(!std::isnan(surfel.pos.x) && !std::isinf(surfel.pos.x));

                surfel.normal = sp.N;
                surfel.area   = area;

                // FIXME: allow only diffuse BRDF parts, check if emit is part of eval
                unsigned char userdata[USER_DATA_SIZE];
                state.userdata = (void *) userdata;

                material_t const* material = sp.material;
                BSDF_t bsdfs;
                material->initBSDF(state, sp, bsdfs);

                surfel.color = estimateAllDirectLight(state, sp, sp.N) + sp.material->emit(state, sp, sp.N);
            }

            Gi_point_inner_node const leaf_node_data = averager.average_leaf(surfels, _precalculated_sf);

            node->set_data(leaf_node_data);
        }
        else
        {
            // bool const treat_as_leaf = (node->get_depth() == treat_depth_as_leaf);
            bool const treat_as_leaf = false;
            assert(treat_depth_as_leaf == -1); // FIXME: add get_depth() to node again

            MyTree::Tree_inner_node const* inner_node = node->get_derived<MyTree::Tree_inner_node>();

            if (!treat_as_leaf)
            {
                std::vector< Gi_point_inner_node const* > children_data = inner_node->get_children_data();
                Gi_point_inner_node const inner_node_data = averager.average_node(children_data);
                node->set_data(inner_node_data);
            }
        }
    }
}


// create an enclosing cube shaped bound
bound_t create_fitting_cube(bound_t const& bound)
{
    int const longest_axis = bound.largestAxis();
    float const longest_extent_2 = bound.get_length(longest_axis) / 2.0f;

    vector3d_t center = vector3d_t(bound.center());

    bound_t result;

    result.a = center - vector3d_t(longest_extent_2 * 1.001f);
    result.g = center + vector3d_t(longest_extent_2 * 1.001f);

    return result;
}

void Pbgi_integrator_t::generate_gi_points(MyTree & tree, int const number_of_samples)
{
    std::vector<triangle_t const*> triangles = get_scene_triangles(scene->get_meshes());

    float const scene_area = get_total_scene_area(triangles);

    float const area_per_sample = scene_area / float(number_of_samples);

    std::cout << "generate_gi_points(): scene_area: " << scene_area << " area_per_sample: " << area_per_sample << std::endl;

    float const radius = std::sqrt(area_per_sample / M_PI) / std::sqrt(2.0f);
    std::cout << "generate_gi_points(): best radius: " << radius << std::endl;

    timer_t timer;
    timer.addEvent("sampling");
    timer.start("sampling");

    bool const use_cdf_sampler = true;

    Scene_sampler_cdf sampler(number_of_samples);
    std::vector<pbgi_sample_t> sampling_points = sampler.generate_samples(triangles);


    std::cout << "generate_gi_points(): samples: " << sampling_points.size() << std::endl;

    std::cout << "samples/desired_samples: " << (sampling_points.size() / float(number_of_samples)) << " time: " << timer.getTime("sampling") << std::endl;

    float const real_covered_area = sampling_points.size() * radius * radius * M_PI;
    float const missing_area_factor = real_covered_area / scene_area; // <= 1
    float const needed_area_per_sample = area_per_sample / missing_area_factor;
    float needed_radius = std::sqrt(needed_area_per_sample / M_PI);

    if (use_cdf_sampler)
    {
        needed_radius *= 2.0f; // additional factor of 2.0 when using the simple CDF sampling
    }

    float const area_estimation = M_PI * needed_radius * needed_radius;

    std::cout << "generate_gi_points(): needed radius: " << needed_radius << std::endl;

    std::cout << "generate_gi_points(): creating points" << std::endl;

    std::tr1::unordered_map<vector3d_t, pbgi_sample_t*, Vector_hash> point_to_sampling_point_map;


    std::vector<vector3d_t> point_data;
    point_data.reserve(sampling_points.size());

    timer.addEvent("point_creation");
    timer.start("point_creation");

    // -------------------------------------------
    // Put the GI points into the tree and create a correspondency map between
    // point and its surface point
    for (int i = 0; i < int(sampling_points.size()); ++i)
    {
        point_to_sampling_point_map[vector3d_t(sampling_points[i].position)] = &sampling_points[i];
        point_data.push_back(vector3d_t(sampling_points[i].position));
    }

    std::cout << "generate_gi_points(): finish tree" << std::endl;

    Tree_subdivision_decider subdivision_decider;

    subdivision_decider.max_num_points = 0;
    subdivision_decider.max_size = needed_radius * 1.0f;

    if (MyTree::Arity == 2)
    {
        tree.build_tree(point_data, scene->getSceneBound(), subdivision_decider);
    }
    else
    {
        bound_t cube_bound = create_fitting_cube(scene->getSceneBound());
        tree.build_tree(point_data, cube_bound, subdivision_decider);
    }

    timer.stop("point_creation");
    // ------------------------------------------


    // ------------------------------------------
    // go through all nodes in post order, generate the SF representation in the children,
    // then go up to the parent, use the children to generate its SF and convert it and delete the children's unconverted etc.

    timer.addEvent("node_data_creation");
    timer.start("node_data_creation");

    std::vector<MyTree::Tree_node*> nodes = tree.get_post_order_queue();
    generate_gi_points_data_2(nodes, point_to_sampling_point_map, area_estimation, needed_radius);

    timer.stop("node_data_creation");

    // ------------------------------------------

    std::cout << "surfel count: " << sampling_points.size() <<
                 " timing, points: " << timer.getTime("point_creation") <<
                 " node data: "      << timer.getTime("node_data_creation") <<
                 std::endl;

}




bool Pbgi_integrator_t::preprocess()
{
    Y_INFO << "PBGI Preprocess" << std::endl;

    bool success = true;
    settings = "";

    std::stringstream ss;
    ss << "#s: "    << surfel_samples <<
          " | rr: " << raster_buffer_resolution;

    settings = ss.str();

    background = scene->getBackground();
    lights = scene->lights;

    _precalculated_sf = precalculate_spherical_harmonics(_spherical_function_area_factory);


    generate_gi_points(_point_tree, surfel_samples);

    Y_INFO << settings << std::endl;

    return success;
}


void Pbgi_integrator_t::cleanup()
{
    _point_tree.clear();
}


void process_surfel(
        Gi_point_surfel const& surfel,
        surfacePoint_t const& rp,
        float const surfel_near_threshold,
        float const disc_scale_factor,
        std::vector<Gi_point_info> & point_infos)
{
    yafaray::vector3d_t surfel_to_rp = (vector3d_t(rp.P) - surfel.pos);

    float const distance = surfel_to_rp.length();

    surfel_to_rp.normalize();

    float const cos_rp_n_gip = rp.N * (-surfel_to_rp);
    if (cos_rp_n_gip < -0.1f) return;

    // float const cos_sp_gip = std::max(gi_point.normal * giToSp, 0.0f);
    float const cos_rp_surfel = surfel.normal * surfel_to_rp;
    bool const back_facing = cos_rp_surfel < 0.0f;

    if (back_facing) return;

    float const visible_area = cos_rp_surfel * surfel.area;
    float const max_area = surfel.area;

    float const solid_angle_real = visible_area / (distance * distance); // std::max(0.001f, distance * distance);

    float const disc_radius = std::sqrt(max_area / M_PI);
    // float const visible_radius = std::sqrt(visible_area / M_PI);

    color_t contribution = surfel.color;

    Gi_point_info point_info;

    point_info.type = Gi_point_info::Far_surfel;
    if (distance < disc_radius * surfel_near_threshold)
    {
        point_info.type = Gi_point_info::Near_surfel;
    }

    point_info.color             = contribution;
    point_info.disc_normal       = surfel.normal;
    point_info.direction         = -surfel_to_rp;
    point_info.radius            = disc_radius * disc_scale_factor;
    point_info.depth             = distance;
    point_info.position          = surfel.pos;
    point_info.receiver_position = vector3d_t(rp.P);
    point_info.solid_angle       = solid_angle_real * disc_scale_factor * disc_scale_factor;

    point_infos.push_back(point_info);
}

void process_node(
        Gi_point_inner_node const* gi_point,
        surfacePoint_t const& receiving_point,
        float const distance,
        vector3d_t const& node_to_rp,
        float const node_radius,
        float const disc_scale_factor,
        std::vector<Gi_point_info> & point_infos)
{
    vector3d_t const& position = gi_point->pos;

    float const visible_area = std::max(0.0f, gi_point->area.get_value(node_to_rp));

    float const real_solid_angle = visible_area / (distance * distance);

    color_t cluster_contribution = gi_point->color.get_value(node_to_rp);
    cluster_contribution *= 1.0f / visible_area;

    Gi_point_info point_info;
    point_info.type               = Gi_point_info::Node;
    point_info.color              = cluster_contribution;
    point_info.direction          = -node_to_rp;
    point_info.depth              = distance;
    point_info.position           = position;
    point_info.receiver_position  = vector3d_t(receiving_point.P);
    point_info.solid_angle        = real_solid_angle * disc_scale_factor * disc_scale_factor;
    point_info.radius             = node_radius * disc_scale_factor;

    point_infos.push_back(point_info);
}



struct Compare_point_info_by_depth
{
    bool operator() (Gi_point_info const& p1, Gi_point_info const& p2)
    {
        return (p1.depth < p2.depth);
    }
};


color_t doPointBasedGiTree_sh_fb(
        Pbgi_integrator_t::MyTree const& tree,
        renderState_t & state,
        surfacePoint_t const& receiving_point,
        float const solid_angle_factor,
        background_t * background,
        vector3d_t const& wo,
        int const raster_buffer_resolution,
        float const surfel_near_threshold,
        float const disc_scale_factor)
{
    typedef Pbgi_integrator_t::MyTree::Tree_node Node;

    Splat_cube_raster_buffer frame_buffer;

    frame_buffer.setup(raster_buffer_resolution);


    std::vector<Gi_point_info> point_infos;

    point_infos.reserve(10000);

    // float const fixed_max_solid_angle = frame_buffer.get_solid_angle(vector3d_t(1, 0, 0)) * solid_angle_factor;
    float const fixed_max_solid_angle = solid_angle_factor;

    float const culling_bias = 0.1f;

    // traverse tree, if solid angle of node > max, traverse into the children

    std::queue<Node const*> queue;
    queue.push(tree.get_root());

    while (!queue.empty())
    {
        Node const* node = queue.front();

        queue.pop();


        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->is_leaf())
        {
            Pbgi_integrator_t::MyTree::Tree_leaf_node const* leaf_node = node->get_derived<Pbgi_integrator_t::MyTree::Tree_leaf_node>();

            std::vector<Gi_point_surfel> const& surfels = leaf_node->get_surfels();

            for (size_t i = 0; i < surfels.size(); ++i)
            {
                Gi_point_surfel const& surfel = surfels[i];
                process_surfel(surfel, receiving_point, surfel_near_threshold, disc_scale_factor, point_infos);
            }
        }
        else
        {
            Gi_point_inner_node const& gi_point = node->get_data();

            if (is_bound_behind_plane(gi_point.bounding_box, vector3d_t(receiving_point.P), receiving_point.N, culling_bias))
            {
                continue;
            }


            vector3d_t node_to_rp = (vector3d_t(receiving_point.P) - gi_point.pos);

            float const distance = node_to_rp.length();

            node_to_rp.normalize();

            float const node_radius = gi_point.bounding_box.get_enclosing_radius();
            float const max_visible_area = node_radius * node_radius * M_PI;
            float const max_node_solid_angle = max_visible_area / (distance * distance);

            float const direction_max_solid_angle = fixed_max_solid_angle;



            if (max_node_solid_angle > direction_max_solid_angle /* solid_angle_threshold */ || distance < node_radius * 1.5f)
            {

                Pbgi_integrator_t::MyTree::Tree_inner_node const* inner_node = node->get_derived<Pbgi_integrator_t::MyTree::Tree_inner_node>();
                std::vector<Pbgi_integrator_t::MyTree::Tree_node*> const& children = inner_node->get_children();

                for (size_t i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
                }
            }
            else
            {
                    process_node(&gi_point, receiving_point, distance, node_to_rp, node_radius, disc_scale_factor, point_infos);
            }
        }
    }


    // sort and splat, then integrate
    std::sort(point_infos.begin(), point_infos.end(), Compare_point_info_by_depth());

    for (unsigned int i = 0; i < point_infos.size(); ++i)
    {
        frame_buffer.add_point(point_infos[i]);
    }

    frame_buffer.add_background(background);

    color_t col = frame_buffer.get_brdf_response(state, receiving_point, wo);

    return col;
}


colorA_t Pbgi_integrator_t::integrate(renderState_t &state, diffRay_t &ray) const
{
    color_t col(0.0);
    float alpha = 0.0;
    surfacePoint_t sp;
    void *o_udat = state.userdata;
    bool oldIncludeLights = state.includeLights;

    if(scene->intersect(ray, sp)) // If it hits
    {
        unsigned char userdata[USER_DATA_SIZE];
        const material_t *material = sp.material;
        BSDF_t bsdfs;

        state.userdata = (void *) userdata;
        vector3d_t wo = -ray.dir;
        if(state.raylevel == 0) state.includeLights = true;

        material->initBSDF(state, sp, bsdfs);

        if (!indirectOnly)
        {
            col += estimateAllDirectLight(state, sp, wo);

            if(bsdfs & BSDF_EMIT)
            {
                col += material->emit(state, sp, wo);
            }
        }

        col += doPointBasedGiTree_sh_fb(
                    _point_tree,
                    state,
                    sp,
                    _solid_angle_factor,
                    background,
                    wo,
                    raster_buffer_resolution,
                    surfel_near_threshold,
                    _disc_scale_factor
                    );

        if (!indirectOnly && !(bsdfs & BSDF_GLOSSY))
        {
            recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);
        }

        float m_alpha = material->getAlpha(state, sp, wo);
        alpha = m_alpha + (1.f - m_alpha) * alpha;
    }
    else // Nothing hit, return background if any
    {
        if (background) col += (*background)(ray, state, false);
    }

    state.userdata = o_udat;
    state.includeLights = oldIncludeLights;
    return colorA_t(col, alpha);
}

integrator_t* Pbgi_integrator_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
    bool transpShad=false;
    bool caustics=false;
    bool do_AO=false;
    int shadowDepth=5;
    int raydepth=5, cDepth=10;
    int search=100, photons=500000;
    int AO_samples = 32;
    double cRad = 0.25;
    double AO_dist = 1.0;
    color_t AO_col(1.f);

    int samples = 10;
    bool indirectOnly = false;

    float maxSolidAngle = 1.0f;
    int fb_resolution = 8;
    float surfel_near_threshold = 2.0f;
    float disc_scale_factor = 1.0f;

    params.getParam("raydepth", raydepth);
    params.getParam("transpShad", transpShad);
    params.getParam("shadowDepth", shadowDepth);
    params.getParam("caustics", caustics);
    params.getParam("photons", photons);
    params.getParam("caustic_mix", search);
    params.getParam("caustic_depth", cDepth);
    params.getParam("caustic_radius", cRad);
    params.getParam("do_AO", do_AO);
    params.getParam("AO_samples", AO_samples);
    params.getParam("AO_distance", AO_dist);
    params.getParam("AO_color", AO_col);

    params.getParam("samplesPerArea", samples);
    params.getParam("indirectOnly", indirectOnly);
    params.getParam("maxSolidAngle", maxSolidAngle);
    params.getParam("fb_resolution", fb_resolution);
    params.getParam("disc_scale_factor", disc_scale_factor);



    Pbgi_integrator_t *inte = new Pbgi_integrator_t(transpShad, shadowDepth, raydepth);
    // caustic settings
    inte->usePhotonCaustics = caustics;
    inte->nCausPhotons = photons;
    inte->nCausSearch = search;
    inte->causDepth = cDepth;
    inte->causRadius = cRad;
    // AO settings
    inte->useAmbientOcclusion = do_AO;
    inte->aoSamples = AO_samples;
    inte->aoDist = AO_dist;
    inte->aoCol = AO_col;

    inte->surfel_samples = samples;
    inte->indirectOnly = indirectOnly;
    inte->raster_buffer_resolution = fb_resolution;
    inte->surfel_near_threshold    = surfel_near_threshold;
    inte->_disc_scale_factor = disc_scale_factor;
    inte->_solid_angle_factor = maxSolidAngle;

    bool use_exact_sh = true;

    std::string spherical_function_factory_type = "SH";
    params.getParam("spherical_function_type", spherical_function_factory_type);


    inte->_spherical_function_color_factory = Spherical_harmonics_factory<color_t>(use_exact_sh);
    inte->_spherical_function_area_factory =  Spherical_harmonics_factory<float>(use_exact_sh);


    Y_INFO <<
              "maxSolidAngle: " << inte->_solid_angle_factor << " " <<
              "fb_resolution: " << inte->raster_buffer_resolution << " " <<
              "surfel_near_threshold: " << inte->surfel_near_threshold << " " <<
              std::endl;


    return inte;
}

extern "C"
{

YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
{
    render.registerFactory("pbgi_simple", Pbgi_integrator_t::factory);
}

}

__END_YAFRAY
