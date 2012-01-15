#ifndef DICTIONARY_CONVERTER_H
#define DICTIONARY_CONVERTER_H

#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

class Spherical_function_converter
{
    public:
    virtual Spherical_function * convert(Spherical_function const* sf) const = 0;
};

class Cube_to_mises_fisher_converter : public Spherical_function_converter
{
public:
    Spherical_function * convert(Spherical_function const* sf) const
    {
        Cube_spherical_function const* csf = static_cast<Cube_spherical_function const*>(sf);

        assert(csf);

        Mises_Fisher_spherical_function * mf_sf = new Mises_Fisher_spherical_function;
        mf_sf->generate_from_cube_spherical_function(csf, 2, 1);

        return mf_sf;
    }
};

class Dictionary_converter : public Spherical_function_converter
{
public:
    Dictionary_converter(std::vector<Spherical_function*> const* dictionary) :
        _dictionary(dictionary)
    {
        for (unsigned int i = 0; i < _dictionary->size(); ++i)
        {
            _unpacked_dictionary.push_back((*_dictionary)[i]->components_to_vectors()[0]);
        }
    }

    Spherical_function * convert(Spherical_function const* sf) const
    {
        Cube_spherical_function const* csf = static_cast<Cube_spherical_function const*>(sf);
        Indexed_Cube_Spherical_function * indexed_csf = new Indexed_Cube_Spherical_function(_unpacked_dictionary, _dictionary, csf);
        return indexed_csf;
    }

    Spherical_function * reconvert(Spherical_function const* sf) const
    {
        Indexed_Cube_Spherical_function const* icsf = static_cast<Indexed_Cube_Spherical_function const*>(sf);

        Cube_spherical_function * csf = new Cube_spherical_function(icsf->get_resolution(), true);

        std::vector< std::vector<float> > const& data = icsf->components_to_vectors();

        csf->from_component_vectors(data);

        return csf;
    }

    static void test()
    {

        /*
        Dictionary_converter c;

        for (int i = 0; i < 16; ++i)
        {
            std::vector< std::vector<float> > center_components(6);

            center_components[0][i] = 1.0f;

            Spherical_function * sf = spherical_function_factory->create();
            sf->from_component_vectors(center_components);
            dictionary.push_back(sf);
        }
        */

        Spherical_function_factory * spherical_function_factory = new Cube_spherical_function_factory(4, false);

        std::vector< std::vector<float> > sf_components(6);

        for (int i = 0; i < 6; ++i)
        {
            for (unsigned int j = 0; j < 48; ++j)
            {
                // sf_components[i].push_back(i / 11.0f);
                sf_components[i].push_back(j);
            }

            for (unsigned int j = 48; j < 64; ++j)
            {
                sf_components[i].push_back(0.5f + i / 11.0f);
            }

        }

        Spherical_function * surfel_sf = spherical_function_factory->create();
        surfel_sf->from_component_vectors(sf_components);

        std::vector<Spherical_function*> dictionary;

        for (int i = 0; i < 6; ++i)
        {
            std::vector< std::vector<float> > center_components(6);
            center_components[0] = sf_components[i];

            Spherical_function * sf = spherical_function_factory->create();
            sf->from_component_vectors(center_components);
            dictionary.push_back(sf);
        }

        Dictionary_converter c(&dictionary);

        Spherical_function * converted_sf = c.convert(surfel_sf);

        Spherical_function * reconverted_sf = c.reconvert(converted_sf);

        Spherical_function * rereconverted_sf = c.convert(reconverted_sf);
    }

private:
    std::vector<Spherical_function*> const* _dictionary;
    std::vector< std::vector<float> > _unpacked_dictionary;
};

__END_YAFRAY

#endif // DICTIONARY_CONVERTER_H
