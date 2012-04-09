#ifndef DICTIONARY_CONVERTER_H
#define DICTIONARY_CONVERTER_H

#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

template <class Data>
class Spherical_function_converter
{
    public:
    virtual Spherical_function<Data> * convert(Spherical_function<Data> const* sf) const = 0;
    virtual std::string get_name() const = 0;
};

template <class Data>
class Cube_to_mises_fisher_converter : public Spherical_function_converter<Data>
{
public:
    Cube_to_mises_fisher_converter(int num_lobes) :
        _num_lobes(num_lobes)
    {}

    Spherical_function<Data> * convert(Spherical_function<Data> const* sf) const
    {
        Cube_spherical_function<Data> const* csf = static_cast<Cube_spherical_function<Data> const*>(sf);

        assert(csf);

        Mises_Fisher_spherical_function<Data> * mf_sf = new Mises_Fisher_spherical_function<Data>;
        mf_sf->generate_from_cube_spherical_function(csf, _num_lobes);

        return mf_sf;
    }

    virtual std::string get_name() const
    {
        return "C_MF";
    }

private:
    int _num_lobes;
};

template <class Data>
class Mises_fisher_to_cube_converter : public Spherical_function_converter<Data>
{
public:
    Mises_fisher_to_cube_converter(int const resolution) :
        _resolution(resolution)
    {}

    Spherical_function<Data> * convert(Spherical_function<Data> const* sf) const
    {
        Mises_Fisher_spherical_function<Data> const* mf_sf = static_cast<Mises_Fisher_spherical_function<Data> const*>(sf);

        assert(mf_sf);

        Cube_spherical_function<Data> * cube_sf = new Cube_spherical_function<Data>(_resolution, false);

        assert(mf_sf->get_lobes().size() == 1);

        cube_sf->sample_cube_cell_area(mf_sf->get_lobes()[0]);

        return cube_sf;
    }

    virtual std::string get_name() const
    {
        return "C_MF";
    }

private:
    int _resolution;
};


template <class Data>
class Spherical_function_indexed_converter : public Spherical_function_converter<Data>
{
public:
    Spherical_function_indexed_converter(std::vector<Spherical_function<Data> *> const* dictionary, bool const keep_original) :
        _dictionary(dictionary),
        _keep_original(keep_original)
    {
        for (unsigned int i = 0; i < _dictionary->size(); ++i)
        {
            _unpacked_dictionary.push_back((*_dictionary)[i]->to_vector());
        }
    }

    Spherical_function<Data> * convert(Spherical_function<Data> const* sf) const
    {
        Indexed_spherical_function<Data> * indexed_sf = new Indexed_spherical_function<Data>(_unpacked_dictionary, _dictionary, sf, _keep_original);
        return indexed_sf;
    }

//    Spherical_function<Data> * convert(Spherical_function<Data> const* sf) const
//    {
//        Cube_spherical_function<Data> const* csf = static_cast<Cube_spherical_function<Data> const*>(sf);
//        Indexed_Cube_Spherical_function<Data> * indexed_csf = new Indexed_Cube_Spherical_function<Data>(_unpacked_dictionary, _dictionary, csf);
//        return indexed_csf;
//    }

    Spherical_function<Data> * reconvert(Spherical_function<Data> const* sf) const
    {
        Indexed_Cube_Spherical_function<Data> const* icsf = static_cast<Indexed_Cube_Spherical_function<Data> const*>(sf);

        Cube_spherical_function<Data> * csf = new Cube_spherical_function<Data>(icsf->get_resolution(), true);

        std::vector< std::vector<float> > const& data = icsf->components_to_vectors();

        csf->from_component_vectors(data);

        return csf;
    }


    virtual std::string get_name() const
    {
        return "C_D";
    }

    static void test()
    {
        Spherical_function_factory<Data> * spherical_function_factory = new Cube_spherical_function_factory<Data>(4, false);

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

        Spherical_function<Data> * surfel_sf = spherical_function_factory->create();
        surfel_sf->from_component_vectors(sf_components);

        std::vector<Spherical_function<Data>*> dictionary;

        for (int i = 0; i < 6; ++i)
        {
            std::vector< std::vector<float> > center_components(6);
            center_components[0] = sf_components[i];

            Spherical_function<Data> * sf = spherical_function_factory->create();
            sf->from_component_vectors(center_components);
            dictionary.push_back(sf);
        }

        Spherical_function_indexed_converter<Data> c(&dictionary);

        Spherical_function<Data> * converted_sf = c.convert(surfel_sf);

        Spherical_function<Data> * reconverted_sf = c.reconvert(converted_sf);

        Spherical_function<Data> * rereconverted_sf = c.convert(reconverted_sf);
    }

private:
    std::vector< Spherical_function<Data> *> const* _dictionary;
    std::vector< std::vector<float> > _unpacked_dictionary;
    bool _keep_original;
};


__END_YAFRAY

#endif // DICTIONARY_CONVERTER_H
