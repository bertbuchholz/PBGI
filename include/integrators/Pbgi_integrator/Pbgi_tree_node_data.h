#ifndef PBGI_TREE_NODE_H
#define PBGI_TREE_NODE_H

#include "spherical_harmonics.h"

#include <core_api/vector3d.h>

__BEGIN_YAFRAY

class Gi_point_surfel
{
public:
    Gi_point_surfel() :
        pos(vector3d_t(0.0f)),
        normal(vector3d_t(0.0f)),
        color(color_t(0.0f)),
        area(0.0f)
    { }

    vector3d_t const& get_normal()
    {
        return normal;
    }

    color_t           get_color()
    {
        return color;
    }

    float             get_area ()
    {
        return area;
    }

    vector3d_t pos;
    vector3d_t normal;
    color_t    color;
    float      area;
};


class Gi_point_inner_node
{
public:
    Gi_point_inner_node() :
        pos(vector3d_t(0.0f))
    { }

    color_t           get_color(vector3d_t const& dir)
    {
        return color.get_value(dir);
    }

    float             get_area (vector3d_t const& dir)
    {
        return area.get_value(dir);
    }

    bound_t    const& get_bound()
    {
        return bounding_box;
    }

// private:
    vector3d_t pos;
    Spherical_harmonics<float>   area;
    Spherical_harmonics<color_t> color;
    bound_t bounding_box;
};



inline
bound_t generate_disc_bounding_box(vector3d_t const& p, vector3d_t const& nu, vector3d_t const& nv, float const radius)
{
    bound_t bounding_box;

    vector3d_t disc_bound[4];
    disc_bound[0] = p + nu * radius + nv * radius;
    disc_bound[1] = p + nu * radius - nv * radius;
    disc_bound[2] = p - nu * radius - nv * radius;
    disc_bound[3] = p - nu * radius + nv * radius;

    bounding_box = bound_t(disc_bound[0], disc_bound[1]);
    bounding_box.include(disc_bound[2]);
    bounding_box.include(disc_bound[3]);

    return bounding_box;
}


class Gi_point_averager
{
    public:
    Gi_point_averager(Spherical_harmonics_factory<color_t> const& sf_color_factory, Spherical_harmonics_factory<float> const& sf_area_factory) :
        _sh_color_factory(sf_color_factory),
        _sh_area_factory(sf_area_factory)
    { }


    Gi_point_inner_node average_leaf(std::vector<Gi_point_surfel> const& points, Cube_raster_buffer< Spherical_harmonics<float> > const& precalculated_sf)
    {
        assert(points.size() > 0);

        Gi_point_inner_node result;
        result.color = _sh_color_factory.create();
        result.area  = _sh_area_factory.create();

        result.bounding_box = bound_t(point3d_t(1e10f, 1e10f, 1e10f), point3d_t(-1e10f, -1e10f, -1e10f));

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            Gi_point_surfel const& surfel = points[i];

            result.pos += surfel.pos;

            vector3d_t NU, NV;
            createCS(surfel.normal, NU, NV);
            float const radius = std::sqrt(surfel.area / M_PI);

            bound_t bounding_box = generate_disc_bounding_box(surfel.pos, NU, NV, radius);

            for (int j = 0; j < 3; ++j)
            {
                result.bounding_box.a[j] = std::min(result.bounding_box.a[j], bounding_box.a[j]);
                result.bounding_box.g[j] = std::max(result.bounding_box.g[j], bounding_box.g[j]);
            }

             Spherical_harmonics<color_t> color_sh = _sh_color_factory.create();
             Spherical_harmonics<float>   area_sh  = _sh_area_factory.create();

             area_sh.get_precalculated_coefficients(precalculated_sf, surfel.area, surfel.normal);

             color_t const scale = surfel.color * surfel.area;
             color_sh.get_precalculated_coefficients(precalculated_sf, scale, surfel.normal);

             result.area.add(area_sh);
             result.color.add(color_sh);


            // result.area += get_precalculated_sh(precalculated_sf, surfel.normal) * surfel.area;
            // result.color += get_precalculated_sh(precalculated_sf, surfel.normal) * (surfel.color * surfel.area);
        }

        result.pos *= 1.0f / float(points.size());

        assert (!std::isnan(result.pos.x) && !std::isinf(result.pos.x));

        /*
        // -------------------- sanity check

        vector3d_t dir(0.4f, 0.3f, 0.5f);
        dir.normalize();

        float surfel_area_sum = 0.0f;
        color_t surfel_color_sum(0.0f);

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            Gi_point_surfel const& p = points[i];

            surfel_area_sum +=  std::max(0.0f, dir * p.normal) * p.area;   // cos() * area
            surfel_color_sum += std::max(0.0f, dir * p.normal) * p.area * p.color; // cos() * area
        }

        float const averaged_area = result.area.get_value(dir);
        color_t const averaged_color = result.color.get_value(dir);

        if (averaged_area > 0.0f && !is_in_range(0.7f, 1.3f, averaged_area / surfel_area_sum))
        {
            std::cout << "Gi_point_averager::average(): not in range: " << surfel_area_sum << " " << averaged_area << std::endl;
        }

        if (surfel_color_sum.energy() > 0.0f && !is_in_range(0.7f, 1.3f, averaged_color.energy() / surfel_color_sum.energy()))
        {
            std::cout << "Gi_point_averager::average(): not in range: " << surfel_area_sum << " " << result.area.get_value(dir) << std::endl;
        }

        // -------------------------
        */

        return result;
    }



    // for new Gi_point_base
    Gi_point_inner_node average_node(std::vector<Gi_point_inner_node const*> const& points) const
    {
        assert(points.size() > 0);

        Gi_point_inner_node result;
        result.color = _sh_color_factory.create();
        result.area  = _sh_area_factory.create();

        result.bounding_box = bound_t(point3d_t(1e10f, 1e10f, 1e10f), point3d_t(-1e10f, -1e10f, -1e10f));

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            Gi_point_inner_node const* p = points[i];
            result.pos    += p->pos;

            for (int j = 0; j < 3; ++j)
            {
                result.bounding_box.a[j] = std::min(result.bounding_box.a[j], p->bounding_box.a[j]);
                result.bounding_box.g[j] = std::max(result.bounding_box.g[j], p->bounding_box.g[j]);
            }

            result.area.add(p->area);
            result.color.add(p->color);
        }

        result.pos *= 1.0f / float(points.size());


        assert (!std::isnan(result.pos.x) && !std::isinf(result.pos.x));


        /*
        // -------------------- test

        vector3d_t dir(0.4f, 0.3f, 0.5f);
        dir.normalize();

        float area_sum = 0;

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            Gi_point_inner_node const& p = *points[i];

            area_sum += p.area.get_value(dir);
        }

        if (result.area.get_value(dir) > 0.0f && !is_in_range(0.9f, 1.1f, result.area.get_value(dir) / area_sum))
        {
            std::cout << "Gi_point_averager::average(): not in range: " << area_sum << " " << result.area.get_value(dir) << std::endl;
        }

        // -------------------------
        */

        return result;
    }

private:
    Spherical_harmonics_factory<color_t> const& _sh_color_factory;
    Spherical_harmonics_factory<float> const& _sh_area_factory;
};

__END_YAFRAY

#endif // PBGI_TREE_NODE_H
