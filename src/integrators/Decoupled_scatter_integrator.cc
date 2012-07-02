#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/integrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <integrators/integr_utils.h>
#include <yafraycore/photon.h>
#include <utilities/mcqmc.h>
#include <yafraycore/scr_halton.h>
#include <vector>
#include <cassert>

__BEGIN_YAFRAY

class YAFRAYPLUGIN_EXPORT Decoupled_scatter_integrator : public volumeIntegrator_t
{
public:
    struct Ray_volume_sample
    {
        Ray_volume_sample() :
            distance(0.0f),
            delta(0.0f),
            sigma_s(0.0f),
            sigma_t(0.0f),
            transmittance(0.0f),
            pdf(0.0)
        {}

        float distance, delta, sigma_s, sigma_t, transmittance, pdf;
    };

    struct Pixel_ray_data_cache
    {
        Pixel_ray_data_cache() : initialized(false)
        {}

        std::vector<Ray_volume_sample> density_samples;
        std::vector<float> density_samples_cdf;
        float t0, t1;
        bool initialized;
    };

    Decoupled_scatter_integrator(float sSize)
    {
        _step_size = sSize;

        Y_INFO << "Decoupled_scatter_integrator: stepSize: " << _step_size << yendl;
    }

    virtual bool preprocess()
    {
        Y_INFO << "Decoupled_scatter_integrator: Preprocessing..." << yendl;

        _lights = scene->lights;
        _volumes = scene->getVolumes();

        // FIXME: should be the number of pixels of the image, 1e6 should be larger than any actual image
        _pixel_ray_data_caches.resize(1e6);

        return true;
    }

    color_t getInScatter(renderState_t& state, ray_t& stepRay, float currentStep) const
    {
        color_t inScatter(0.f);
        surfacePoint_t sp;
        sp.P = stepRay.from;

        ray_t lightRay;
        lightRay.from = sp.P;

        float const inv_num_volumes = 1.0f / float(_volumes.size());

        for(std::vector<light_t *>::const_iterator l=_lights.begin(); l!=_lights.end(); ++l)
        {
            color_t lcol(0.0);

            // handle lights with delta distribution, e.g. point and directional lights
            if( (*l)->diracLight() )
            {
                if( (*l)->illuminate(sp, lcol, lightRay) )
                {
                    // ...shadowed...
                    lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                    if (lightRay.tmax < 0.f) lightRay.tmax = 1e10; // infinitely distant light
                    bool shadowed = scene->isShadowed(state, lightRay);
                    if (!shadowed)
                    {
                        float lightTr = 0.0f;
                        // replace lightTr with precalculated attenuation

                        // replaced by
                        color_t lightstepTau(0.f);
                        for (unsigned int i = 0; i < _volumes.size(); i++)
                        {
                            VolumeRegion* vr = _volumes.at(i);
                            float t0Tmp = -1, t1Tmp = -1;
                            if (_volumes.at(i)->intersect(lightRay, t0Tmp, t1Tmp))
                            {
                                lightstepTau += vr->tau(lightRay, currentStep, 0.f) * inv_num_volumes;
                            }
                        }
                        // transmittance from the point p in the volume to the light (i.e. how much light reaches p)
                        lightTr = fExp(-lightstepTau.energy());

                        lightTr *= inv_num_volumes;

                        inScatter += lightTr * lcol;
                    }
                }
            }
            else // area light and suchlike
            {
                int n = (*l)->nSamples() >> 2; // samples / 4
                if (n < 1) n = 1;
                float iN = 1.f / (float)n; // inverse of n
                color_t ccol(0.0);
                float lightTr = 0.0f;
                lSample_t ls;

                for(int i=0; i<n; ++i)
                {
                    // ...get sample val...
                    ls.s1 = (*state.prng)();
                    ls.s2 = (*state.prng)();

                    if((*l)->illumSample(sp, ls, lightRay))
                    {
                        // ...shadowed...
                        lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                        if (lightRay.tmax < 0.f) lightRay.tmax = 1e10; // infinitely distant light
                        bool shadowed = scene->isShadowed(state, lightRay);
                        if(!shadowed) {
                            ccol += ls.col / ls.pdf;

                            // replaced by
                            color_t lightstepTau(0.f);
                            for (unsigned int i = 0; i < _volumes.size(); i++)
                            {
                                VolumeRegion* vr = _volumes.at(i);
                                float t0Tmp = -1, t1Tmp = -1;
                                if (_volumes.at(i)->intersect(lightRay, t0Tmp, t1Tmp))
                                {
                                    lightstepTau += vr->tau(lightRay, currentStep * 4.f, 0.0f);
                                }
                            }
                            // transmittance from the point p in the volume to the light (i.e. how much light reaches p)
                            lightTr += fExp(-lightstepTau.energy()) * inv_num_volumes;
                        }
                    }
                    lightTr *= inv_num_volumes;
                } // end of area light sample loop

                lightTr *= iN;

                ccol = ccol * iN;
                inScatter += lightTr * ccol;
            } // end of area lights loop
        }

        return inScatter;
    }


    // optical thickness, absorption, attenuation, extinction
    virtual colorA_t transmittance(renderState_t &state, ray_t &ray) const
    {
        colorA_t Tr(1.f);
        //return Tr;

        if (_volumes.size() == 0) return Tr;

        for (unsigned int i = 0; i < _volumes.size(); i++)
        {
            VolumeRegion const* vr = _volumes[i];
            float t0 = -1, t1 = -1;
            if (vr->intersect(ray, t0, t1))
            {
                float random = (*state.prng)();
                color_t opticalThickness = vr->tau(ray, _step_size, random);
                Tr *= colorA_t(fExp(-opticalThickness.energy()));
            }
        }

        return Tr;
    }

    void get_sigmas_from_all_volumes(std::vector<VolumeRegion*> const& volumes, vector3d_t const& pos, float & sigma_s, float & sigma_t)
    {
        sigma_s = 0.0f;
        sigma_t = 0.0f;

        for (unsigned int i = 0; i < _volumes.size(); i++)
        {
            VolumeRegion const* vr = _volumes[i];

            sigma_s += vr->sigma_s(pos, vector3d_t()).energy();
            sigma_t += vr->sigma_t(pos, vector3d_t()).energy();
        }
    }



    std::vector<Ray_volume_sample> calc_ray_samples(ray_t const& ray, float & t0, float & t1) const
    {
        std::vector<Ray_volume_sample> samples;

        float const delta = _step_size;

        bool const hit = (ray.tmax > 0.f);


        t0 = 1e10f, t1 = -1e10f;

        // find min t0 and max t1
        for (unsigned int i = 0; i < _volumes.size(); i++)
        {
            float t0_tmp = 0.f, t1_tmp = 0.f;
            VolumeRegion const* vr = _volumes[i];

            if (!vr->intersect(ray, t0_tmp, t1_tmp)) continue;

            if (hit && ray.tmax < t0_tmp) continue;

            if (t0_tmp < 0.f) t0_tmp = 0.f;

            if (hit && ray.tmax < t1_tmp) t1_tmp = ray.tmax;

            if (t1_tmp > t1) t1 = t1_tmp;
            if (t0_tmp < t0) t0 = t0_tmp;
        }

        float dist = t1-t0;
        if (dist < 1e-3f) return samples;

        float t = t0 + delta;

        Ray_volume_sample sample_0;
        sample_0.delta = delta;
        sample_0.sigma_s = 0.0f;
        sample_0.sigma_t = 0.0f;
        sample_0.transmittance = 1.0f;
        sample_0.distance = t0;
        samples.push_back(sample_0);

        while (t < t1)
        {
            vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * t;

            Ray_volume_sample sample;

            sample.delta = delta;
            sample.distance = t;

            for (unsigned int i = 0; i < _volumes.size(); i++)
            {
                VolumeRegion const* vr = _volumes[i];

                sample.sigma_s += vr->sigma_s(sample_pos, vector3d_t()).energy();
                sample.sigma_t += vr->sigma_t(sample_pos, vector3d_t()).energy();
            }

            Ray_volume_sample const& prev_sample = samples.back();
            sample.transmittance = prev_sample.transmittance * fExp(-prev_sample.delta * prev_sample.sigma_t);
            sample.pdf = sample.delta * sample.sigma_s * sample.transmittance;

            samples.push_back(sample);

            t += delta;
        }

        return samples;
    }



    // emission and in-scattering
    colorA_t integrate_discrete_density(renderState_t &state, ray_t &ray, float & final_pdf) const
    {
        colorA_t result(0.f, 1.0f);

        if (_volumes.size() == 0) return result;

        float t0, t1;

        if (!_pixel_ray_data_caches[state.pixelNumber].initialized)
        {
            _pixel_ray_data_caches[state.pixelNumber].density_samples = calc_ray_samples(ray, t0, t1);
        }

        std::vector<Ray_volume_sample> const& samples = _pixel_ray_data_caches[state.pixelNumber].density_samples;

        if (samples.size() < 1) return result;

//        return colorA_t(1.0f) * samples.back().transmittance;

        std::vector<float> samples_cdf(samples.size());
        samples_cdf[0] = samples[0].pdf;

        for (unsigned int i = 1; i < samples_cdf.size(); ++i)
        {
            samples_cdf[i] = samples_cdf[i - 1] + samples[i].pdf;
        }

        float const pdf_sum = samples_cdf.back();


        Halton halton2(2);
        halton2.setStart(state.pixelSample + state.samplingOffs);

        float const ksi_0 = halton2.getNext();

        int sample_index = std::lower_bound(samples_cdf.begin(), samples_cdf.end(), ksi_0 * pdf_sum) - samples_cdf.begin();

        Ray_volume_sample const& sample = samples[sample_index];




        Halton halton3(3);
        halton3.setStart(state.pixelSample + state.samplingOffs);

        float const ksi_1 = halton3.getNext();

        float const a = sample.distance;
        float const b = sample.distance + sample.delta;
        float const t = a - (1.0f / sample.sigma_t) * fLog(1.0f - ksi_1 * (1.0f - fExp(-(b - a) * sample.sigma_t)));

//        float const pdf = sample.sigma_t / (fExp(-(t - a) * sample.sigma_t) - fExp(-(t - b) * sample.sigma_t));


        // vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * sample.distance;
        vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * t;


        // randomly choose light
        light_t const* light = _lights[0];

        vector3d_t dummy_wo;
        lSample_t light_sample;
        light_sample.sp = new surfacePoint_t;
        color_t const light_color = light->emitSample(dummy_wo, light_sample);
        vector3d_t const& light_pos = vector3d_t(light_sample.sp->P);

        vector3d_t sample_to_light = (light_pos - sample_pos);
        float const dist_sample_to_light = sample_to_light.length();
        sample_to_light.normalize();

        ray_t lightRay(sample_pos, sample_to_light, 0.0f, dist_sample_to_light);
        bool const shadowed = scene->isShadowed(state, lightRay);

        if (!shadowed)
        {
            float const pdf = sample.pdf / pdf_sum;

            color_t const s_s = sample.sigma_s;
            // color_t const s_t = sample.sigma_t;

            float const transmittance_to_light = transmittance(state, lightRay).energy();

            assert(t >= sample.distance);
            float total_transmittance = sample.transmittance * fExp(t - sample.distance) * transmittance_to_light;

            result = s_s * total_transmittance * light_color / (dist_sample_to_light * dist_sample_to_light);
            result *= 1.0f / pdf;

            final_pdf = pdf;
        }
        else
        {
            final_pdf = 0.0f;
        }

        delete light_sample.sp;

        return result;
    }


    float project_point_on_ray(vector3d_t const& point, ray_t const& ray) const
    {
        vector3d_t const& p0 = vector3d_t(ray.from);
        vector3d_t const& p1 = vector3d_t(ray.from) + ray.dir;

        vector3d_t v = p1 - p0;
        vector3d_t w = point - p0;

        float const c1 = w * v;
        float const c2 = v * v;

        float const u = c1 / c2;

        return u;
    }





    // emission and in-scattering
    colorA_t integrate_equi_angular(renderState_t &state, ray_t &ray, float & final_pdf) const
    {
        float t0 = 1e10f, t1 = -1e10f;

        colorA_t result(0.f, 1.0f);

        if (_volumes.size() == 0) return result;

        bool hit = (ray.tmax > 0.f);

        // find min t0 and max t1
        for (unsigned int i = 0; i < _volumes.size(); i++)
        {
            float t0Tmp = 0.f, t1Tmp = 0.f;
            VolumeRegion const* vr = _volumes[i];

            if (!vr->intersect(ray, t0Tmp, t1Tmp)) continue;

            if (hit && ray.tmax < t0Tmp) continue;

            if (t0Tmp < 0.f) t0Tmp = 0.f;

            if (hit && ray.tmax < t1Tmp) t1Tmp = ray.tmax;

            if (t1Tmp > t1) t1 = t1Tmp;
            if (t0Tmp < t0) t0 = t0Tmp;
        }

        float dist = t1-t0;
        if (dist < 1e-3f) return result;

        // randomly choose light
        light_t const* light = _lights[0];

        VolumeRegion const* vr = _volumes[0];

        vector3d_t dummy_wo;
        lSample_t light_sample;
        light_sample.sp = new surfacePoint_t;
        color_t const light_color = light->emitSample(dummy_wo, light_sample);
        vector3d_t const& light_pos = vector3d_t(light_sample.sp->P);
        float const closest_t = project_point_on_ray(light_pos, ray);
        vector3d_t const closest_pos = vector3d_t(ray.from) + ray.dir * closest_t;

//        std::cout << "lpos: " << light_pos << " "
//                  << "closest_pos: " << closest_pos << " "
//                  << std::endl;


        float const D = (closest_pos - light_pos).length();
        float const a = t0 - closest_t;
        float const b = t1 - closest_t;

        float const theta_a = std::atan(a / D);
        float const theta_b = std::atan(b / D);

        Halton halton2(2);
        halton2.setStart(state.pixelSample + state.samplingOffs);

//        float const ksi = (*state.prng)();
        float const ksi = halton2.getNext();
        float const t = D * std::tan((1.0f - ksi) * theta_a + ksi * theta_b);
        vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * (closest_t + t);

        vector3d_t sample_to_light = (light_pos - sample_pos);
        float const dist_sample_to_light = std::sqrt(D * D + t * t);
        sample_to_light.normalize();
        ray_t lightRay(sample_pos, sample_to_light, 0.0f, dist_sample_to_light);
        bool const shadowed = scene->isShadowed(state, lightRay);

        if (!shadowed)
        {
            float const pdf = D / ((theta_b - theta_a) * (D * D + t * t));

            color_t const s_s = vr->sigma_s(vector3d_t(sample_pos), vector3d_t());
            color_t const s_t = vr->sigma_t(vector3d_t(sample_pos), vector3d_t());

            assert(t + closest_t + std::sqrt(D * D + t * t) >= 0.0f);
            result = s_s * fExp(-1.0f * s_t.energy() * (t + closest_t + std::sqrt(D * D + t * t))) * light_color;
            result *= 1.0f / (D * D + t * t);
            result *= 1.0f / pdf;

            final_pdf = pdf;
        }
        else
        {
            final_pdf = 0.0f;
        }

        delete light_sample.sp;

        return result;
    }



    colorA_t integrate_equi_angular_homogeneous(renderState_t &state, ray_t &ray, float & final_pdf) const
    {
        float t0 = 1e10f, t1 = -1e10f;

        colorA_t result(0.f, 1.0f);

        if (_volumes.size() == 0) return result;

        bool hit = (ray.tmax > 0.f);

        // find min t0 and max t1
        for (unsigned int i = 0; i < _volumes.size(); i++)
        {
            float t0Tmp = 0.f, t1Tmp = 0.f;
            VolumeRegion const* vr = _volumes[i];

            if (!vr->intersect(ray, t0Tmp, t1Tmp)) continue;

            if (hit && ray.tmax < t0Tmp) continue;

            if (t0Tmp < 0.f) t0Tmp = 0.f;

            if (hit && ray.tmax < t1Tmp) t1Tmp = ray.tmax;

            if (t1Tmp > t1) t1 = t1Tmp;
            if (t0Tmp < t0) t0 = t0Tmp;
        }

        float dist = t1-t0;
        if (dist < 1e-3f) return result;

        // randomly choose light
        light_t const* light = _lights[0];

        VolumeRegion const* vr = _volumes[0];

        vector3d_t dummy_wo;
        lSample_t light_sample;
        light_sample.sp = new surfacePoint_t;
        color_t const light_color = light->emitSample(dummy_wo, light_sample);
        vector3d_t const& light_pos = vector3d_t(light_sample.sp->P);
        float const closest_t = project_point_on_ray(light_pos, ray);
        vector3d_t const closest_pos = vector3d_t(ray.from) + ray.dir * closest_t;

//        std::cout << "lpos: " << light_pos << " "
//                  << "closest_pos: " << closest_pos << " "
//                  << std::endl;


        float const D = (closest_pos - light_pos).length();
        float const a = t0 - closest_t;
        float const b = t1 - closest_t;

        float const theta_a = std::atan(a / D);
        float const theta_b = std::atan(b / D);

        Halton halton2(2);
        halton2.setStart(state.pixelSample + state.samplingOffs);

//        float const ksi = (*state.prng)();
        float const ksi = halton2.getNext();
        float const t = D * std::tan((1.0f - ksi) * theta_a + ksi * theta_b);
        vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * (closest_t + t);

        vector3d_t sample_to_light = (light_pos - sample_pos);
        float const dist_sample_to_light = std::sqrt(D * D + t * t);
        sample_to_light.normalize();
        ray_t lightRay(sample_pos, sample_to_light, 0.0f, dist_sample_to_light);
        bool const shadowed = scene->isShadowed(state, lightRay);

        if (!shadowed)
        {
            float const pdf = D / ((theta_b - theta_a) * (D * D + t * t));

            color_t const s_s = vr->sigma_s(vector3d_t(sample_pos), vector3d_t());
            color_t const s_t = vr->sigma_t(vector3d_t(sample_pos), vector3d_t());

            assert(t + closest_t + std::sqrt(D * D + t * t) >= 0.0f);
            result = s_s * fExp(-1.0f * s_t.energy() * (t + closest_t + std::sqrt(D * D + t * t))) * light_color;
            result *= 1.0f / (D * D + t * t);
            result *= 1.0f / pdf;

            final_pdf = pdf;
        }
        else
        {
            final_pdf = 0.0f;
        }

        delete light_sample.sp;

        return result;
    }


    float sample_from_density(renderState_t const& state, std::vector<Ray_volume_sample> const& samples, std::vector<float> const& samples_cdf) const
    {
        float const pdf_sum = samples_cdf.back();

        Halton halton2(2);
        halton2.setStart(state.pixelSample + state.samplingOffs);

        float const ksi_0 = halton2.getNext();

        int const sample_index = std::min(int(samples.size() - 1), std::lower_bound(samples_cdf.begin(), samples_cdf.end(), ksi_0 * pdf_sum) - samples_cdf.begin());

        Ray_volume_sample const& sample = samples[sample_index];


        Halton halton3(3);
        halton3.setStart(state.pixelSample + state.samplingOffs);

        float const ksi_1 = halton3.getNext();

        float const a = sample.distance;
        float const b = sample.distance + sample.delta;
        float t = a - (1.0f / sample.sigma_t) * std::log(1.0f - ksi_1 * (1.0f - std::exp(-(b - a) * sample.sigma_t)));

        t = std::max(a + 0.001f * sample.delta, std::min(b - 0.001f * sample.delta, t));


        { // debug sanity
            assert(t >= a && t <= b);

            int const sample_index_new = (t - samples[0].distance) / samples[0].delta;

            // std::cout << sample_index << " " << sample_index_new << " " << t << " " << samples[0].distance << std::endl;

            if (!(t >= samples[sample_index_new].distance && sample_index_new >= 0 && sample_index_new < samples.size()))
            {
                std::cout << sample_index << " " << sample_index_new << " " << t << " " << samples[0].distance << " " << a << " " << b << " " << sample.delta << " " << samples.size() << std::endl;
            }

            assert(t >= samples[sample_index_new].distance && sample_index_new >= 0 && sample_index_new < samples.size());

            float const t_offset = (t - samples[sample_index_new].distance) / samples[0].delta;

            assert(t_offset >= 0.0f && t_offset <= 1.0f);
        }

        return t;
    }

    float pdf_from_density(float const t, std::vector<Ray_volume_sample> const& samples) const
    {
        // int const sample_index = std::min(int((t - samples[0].distance) / samples[0].delta), int(samples.size() - 1));
        int const sample_index = (t - samples[0].distance) / samples[0].delta;
        float const t_offset = std::max(0.0f, std::min(1.0f, (t - samples[sample_index].distance) / samples[0].delta));

        // assert(t >= samples[sample_index].distance);
        assert(t_offset >= 0.0f && t_offset <= 1.0f);

        return samples[sample_index].pdf;

        //assert(sample_index + 1 < samples.size());
        // return samples[sample_index].pdf * (1.0f - t_offset) + samples[sample_index + 1].pdf * t_offset;
    }


    std::vector<float> calc_density_samples_pdf(std::vector<Ray_volume_sample> const& samples) const
    {
        if (samples.size() < 2) return std::vector<float>();

        std::vector<float> samples_cdf(samples.size());
        samples_cdf[0] = samples[0].pdf;

        for (unsigned int i = 1; i < samples_cdf.size(); ++i)
        {
            samples_cdf[i] = samples_cdf[i - 1] + samples[i].pdf;
        }

        return samples_cdf;
    }


    float pdf_from_light(ray_t const& ray, light_t const* light, float const t0, float const t1, vector3d_t const sample_pos) const
    {
        lSample_t light_sample;
        vector3d_t dummy_wo;
        light_sample.sp = new surfacePoint_t;
        light->emitSample(dummy_wo, light_sample);
        vector3d_t const& light_pos = vector3d_t(light_sample.sp->P);

        float const closest_t = project_point_on_ray(light_pos, ray);
        vector3d_t const closest_pos = vector3d_t(ray.from) + ray.dir * closest_t;


        float const D = (closest_pos - light_pos).length();
        float const a = t0 - closest_t;
        float const b = t1 - closest_t;

        float const theta_a = std::atan(a / D);
        float const theta_b = std::atan(b / D);

        float const dist_sample_to_light_sqr = (light_pos - sample_pos).lengthSqr();

        float const pdf = D / ((theta_b - theta_a) * dist_sample_to_light_sqr);

        delete light_sample.sp;

        return pdf;
    }


    float sample_from_light(renderState_t const& state, ray_t const& ray, light_t const* light, float const t0, float const t1) const
    {
        vector3d_t dummy_wo;
        lSample_t light_sample;
        light_sample.sp = new surfacePoint_t;
        light->emitSample(dummy_wo, light_sample);
        vector3d_t const& light_pos = vector3d_t(light_sample.sp->P);
        float const closest_t = project_point_on_ray(light_pos, ray);
        vector3d_t const closest_pos = vector3d_t(ray.from) + ray.dir * closest_t;

//        std::cout << "lpos: " << light_pos << " "
//                  << "closest_pos: " << closest_pos << " "
//                  << std::endl;


        float const D = (closest_pos - light_pos).length();
        float const a = t0 - closest_t;
        float const b = t1 - closest_t;

        float const theta_a = std::atan(a / D);
        float const theta_b = std::atan(b / D);

        Halton halton2(2);
        halton2.setStart(state.pixelSample + state.samplingOffs);

        float const ksi = halton2.getNext();
        float const t = D * std::tan((1.0f - ksi) * theta_a + ksi * theta_b);

        delete light_sample.sp;

        return closest_t + t;
    }


    colorA_t evaluate_light(renderState_t & state, light_t const* light, vector3d_t const& sample_pos, std::vector<Ray_volume_sample> const& samples, float const t) const
    {
        colorA_t color(0.0f);

        // vector3d_t dummy_wo;
        // lSample_t light_sample;
        // light_sample.sp = new surfacePoint_t;
        color_t light_color; // = light->emitSample(dummy_wo, light_sample);

        // vector3d_t const& light_pos = vector3d_t(light_sample.sp->P);

//        vector3d_t sample_to_light = (light_pos - sample_pos);
//        float const dist_sample_to_light = sample_to_light.length();
//        sample_to_light.normalize();

        // ray_t lightRay(sample_pos, sample_to_light, 0.0f, dist_sample_to_light);

        surfacePoint_t sp;
        sp.P = sample_pos;

        ray_t lightRay;
        lightRay.from = sp.P;

        if(light->diracLight())
        {
            if(light->illuminate(sp, light_color, lightRay))
            {
                // ...shadowed...
                lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                if (lightRay.tmax < 0.f) lightRay.tmax = 1e10; // infinitely distant light
            }
        }
        else // area light and suchlike
        {
            int n = 1;
            lSample_t ls;

            for(int i = 0; i < n; ++i)
            {
                // ...get sample val...
                Halton halton2(2);
                Halton halton3(2);

                halton2.setStart(state.pixelSample + state.samplingOffs);
                halton3.setStart(state.pixelSample + state.samplingOffs);

//                ls.s1 = (*state.prng)();
//                ls.s2 = (*state.prng)();

                ls.s1 = halton2.getNext();
                ls.s2 = halton3.getNext();

                if(light->illumSample(sp, ls, lightRay))
                {
                    // ...shadowed...
                    lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                    if (lightRay.tmax < 0.f) lightRay.tmax = 1e10; // infinitely distant light
                    light_color = ls.col;
                }
            }
        }


        vector3d_t light_pos = vector3d_t(lightRay.from) + lightRay.dir * lightRay.tmax;
        vector3d_t sample_to_light = (light_pos - sample_pos);
        float const dist_sample_to_light = sample_to_light.length();


        bool const shadowed = scene->isShadowed(state, lightRay);

        if (!shadowed)
        {
            int const sample_index = std::min(int((t - samples[0].distance) / samples[0].delta), int(samples.size() - 1));

            // std::cout << "eval light: " << sample_index << std::endl;

            Ray_volume_sample const& sample_0 = samples[sample_index];

            color_t const s_s = sample_0.sigma_s;
            // color_t const s_t = sample.sigma_t;

            float const transmittance_to_light = transmittance(state, lightRay).energy();

            // assert(t >= sample_0.distance); // FIXME: this fails sometimes, t is then very close to sample_0.distance
            float total_transmittance = sample_0.transmittance * fExp(t - sample_0.distance) * transmittance_to_light;

            color = s_s * total_transmittance * light_color / (dist_sample_to_light * dist_sample_to_light);

            /*
            assert(sample_index < samples.size());

            float const t_offset = std::max(0.0f, std::min(1.0f, (t - samples[sample_index].distance) / samples[0].delta));

            assert(t_offset >= 0.0f && t_offset <= 1.0f);

            color_t const s_s = samples[sample_index].sigma_s * (1.0f - t_offset) + samples[sample_index + 1].sigma_s * t_offset;

            // color_t const s_s = sample.sigma_s;
            // color_t const s_t = sample.sigma_t;

            float const transmittance_to_light = transmittance(state, lightRay).energy();

            // assert(t >= samples[sample_index].distance);

            float const sample_transmittance = samples[sample_index].transmittance * (1.0f - t_offset) + samples[sample_index + 1].transmittance * t_offset;

            float total_transmittance = sample_transmittance * transmittance_to_light;

            color = s_s * total_transmittance * light_color / (dist_sample_to_light * dist_sample_to_light);
            */
        }

        // delete light_sample.sp;

        return color;
    }


    colorA_t integrate_mis(renderState_t &state, ray_t &ray) const
    {
        colorA_t result(0.f, 1.0f);

        if (_volumes.size() == 0) return result;


        assert(state.pixelNumber < _pixel_ray_data_caches.size());

        float t0, t1;

        if (!_pixel_ray_data_caches[state.pixelNumber].initialized)
        {
            _pixel_ray_data_caches[state.pixelNumber].density_samples = calc_ray_samples(ray, t0, t1);
            _pixel_ray_data_caches[state.pixelNumber].density_samples_cdf = calc_density_samples_pdf(_pixel_ray_data_caches[state.pixelNumber].density_samples);
            _pixel_ray_data_caches[state.pixelNumber].initialized = true;
            _pixel_ray_data_caches[state.pixelNumber].t0 = t0;
            _pixel_ray_data_caches[state.pixelNumber].t1 = t1;
        }

        std::vector<Ray_volume_sample> const& samples = _pixel_ray_data_caches[state.pixelNumber].density_samples;

        if (samples.size() < 2) return result;

        std::vector<float>             const& samples_cdf = _pixel_ray_data_caches[state.pixelNumber].density_samples_cdf;

        t0 = _pixel_ray_data_caches[state.pixelNumber].t0;
        t1 = _pixel_ray_data_caches[state.pixelNumber].t1;


//        std::vector<Ray_volume_sample> samples     = calc_ray_samples(ray, t0, t1);

//        if (samples.size() < 2) return result;

//        std::vector<float>             samples_cdf = calc_density_samples_pdf(samples);


        // iterate over PDFs
        //     1. draw sample from PDF (n + 1: n <- number of lights, one sample from density and 1 from each light's PDF)
        //     2. calc sample position and calculate inscatter from all lights
        //     3. divide by the average of all PDFs at the sample position


        int const num_pdfs = 1 + _lights.size();

        // sample from density's pdf
        {
            // Take sample
            float const t = sample_from_density(state, samples, samples_cdf);
            vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * t;

            // calc average PDF
            float const pdf_sum = samples_cdf.back();
            float average_pdf = pdf_from_density(t, samples) / pdf_sum;

            for (size_t i = 0; i < _lights.size(); ++i)
            {
                average_pdf += pdf_from_light(ray, _lights[i], t0, t1, sample_pos);
            }

            average_pdf *= 1.0f / float(num_pdfs);

            // calc estimator (evaluate all lights at sample_pos)

            color_t inscatter(0.0f);

            for (size_t i = 0; i < _lights.size(); ++i)
            {
                inscatter += evaluate_light(state, _lights[i], sample_pos, samples, t);
            }

            // update result

            result += inscatter / average_pdf;
        }

        // sample from lights' pdfs
        for (size_t i = 0; i < _lights.size(); ++i)
        {
            // Take sample
            float const t = sample_from_light(state, ray, _lights[i], t0, t1);
            vector3d_t const sample_pos = vector3d_t(ray.from) + ray.dir * t;

            // calc average PDF
            float const pdf_sum = samples_cdf.back();
            float average_pdf = pdf_from_density(t, samples) / pdf_sum;

            for (size_t j = 0; j < _lights.size(); ++j)
            {
                average_pdf += pdf_from_light(ray, _lights[j], t0, t1, sample_pos);
            }

            average_pdf *= 1.0f / float(num_pdfs);

            // calc estimator (evaluate all lights at sample_pos)

            color_t inscatter(0.0f);

            for (size_t j = 0; j < _lights.size(); ++j)
            {
                inscatter += evaluate_light(state, _lights[j], sample_pos, samples, t);
            }

            // update result

            result += inscatter / average_pdf;
        }

        return result * (1.0f / float(num_pdfs));
    }


    virtual colorA_t integrate(renderState_t &state, ray_t &ray) const
    {
        // return integrate_equi_angular_homogeneous(state, ray, distance_pdf);
        // return integrate_discrete_density(state, ray, distance_pdf);
        return integrate_mis(state, ray);
    }


    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render)
    {
        float step_size = 1.f;
        params.getParam("stepSize", step_size);
        Decoupled_scatter_integrator * inte = new Decoupled_scatter_integrator(step_size);
        return inte;
    }

    private:
    std::vector<VolumeRegion*> _volumes;
    std::vector<light_t*> _lights;

    //mutable std::vector<std::vector<Ray_volume_sample> > _density_samples;
    //mutable std::vector<std::vector<float> > _density_samples_cdf;
    mutable std::vector<Pixel_ray_data_cache> _pixel_ray_data_caches;
    float _step_size;

};

extern "C"
{

    YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
    {
        render.registerFactory("Decoupled_scatter_integrator", Decoupled_scatter_integrator::factory);
    }

}

__END_YAFRAY
