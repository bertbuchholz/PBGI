#ifndef KMEANS_HPP
#define KMEANS_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>

#include <boost/random.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "distance.hpp"

typedef std::vector<float> Word;
//typedef imdb::l1norm<Word> Distance_function;
typedef imdb::l2norm<Word> Distance_function;

template <class index_t, class collection_t>
void kmeans_init_random(std::vector<index_t>& centers, const collection_t& collection, std::size_t numclusters)
{
    assert(collection.size() >= numclusters);

    centers.resize(collection.size());
    for (std::size_t i = 0; i < centers.size(); i++) centers[i] = i;
    std::random_shuffle(centers.begin(), centers.end());
    centers.resize(numclusters);
}

template <class index_t, class collection_t, class dist_fn>
void kmeans_init_plusplus(std::vector<index_t>& result, const collection_t& collection, std::size_t numclusters, const dist_fn& distfn)
{
    assert(numclusters > 0);
    assert(collection.size() >= numclusters);

    typedef typename collection_t::value_type item_t;
    typedef boost::mt19937                    rng_t;
    typedef boost::uniform_real<double>       unirand_t;

    rng_t rng;

    boost::variate_generator<rng_t&, unirand_t> unirand(rng, unirand_t(0.0, 1.0));

    std::size_t numtrials = 2 + std::log(numclusters);

    // add first cluster, randomly chosen
    std::set<index_t> centers;
    index_t first = unirand() * collection.size();
    centers.insert(first);

    // compute distance between first cluster center and all others
    // and accumulate the distances that gives the current potential
    std::vector<double> dists(collection.size());
    double potential = 0.0;
    for (std::size_t i = 0; i < collection.size(); i++)
    {
        double d = distfn(collection[first], collection[i]);
        dists[i] = d*d;
        potential += dists[i];
    }
    std::cout << "kmeans++ init: numclusters=" << numclusters << " numtrials=" << numtrials << " collection.size=" << collection.size() << " init pot=" << potential << std::endl;

    // iteratively add centers
    for (std::size_t c = 1; c < numclusters; c++)
    {
        double min_potential = std::numeric_limits<double>::max();
        std::size_t best_index = 0;

        for (std::size_t i = 0; i < numtrials; i++)
        {
            std::size_t index;

            // get new center
            double r = unirand() * potential;
            for (index = 0; index < collection.size()-1 && r > dists[index]; index++)
            {
                r -= dists[index];
            }

            while (centers.count(index) > 0) index = (index + 1) % collection.size();

            // recompute potential
            double p = 0.0;
            for (std::size_t k = 0; k < collection.size(); k++)
            {
                double d = distfn(collection[index], collection[k]);
                p += std::min(dists[k], d*d);
            }

            if (p < min_potential)
            {
                min_potential = p;
                best_index = index;
            }
        }

        for (std::size_t i = 0; i < collection.size(); i++)
        {
            double d = distfn(collection[best_index], collection[i]);
            dists[i] = d*d;
        }

        potential = min_potential;

        centers.insert(best_index);

        std::cout << "new center " << c << ": potential=" << potential << " index=" << best_index << std::endl;
    }

    std::copy(centers.begin(), centers.end(), std::back_inserter(result));
}

enum KmeansInitAlgorithm
{
    KmeansInitRandom,
    KmeansInitPlusPlus
};

template <class collection_t, class dist_fn>
class kmeans
{
    typedef boost::mutex               mutex_t;
    typedef boost::lock_guard<mutex_t> locker_t;

    typedef typename collection_t::value_type item_t;

public:

    kmeans(const collection_t& collection, std::size_t numclusters, KmeansInitAlgorithm initalgorithm = KmeansInitRandom, const dist_fn& distfn = dist_fn())
        : _collection(collection), _distfn(distfn), _centers(numclusters), _clusters(collection.size())
    {
        // get initial centers
        std::vector<std::size_t> initindices;
        if (initalgorithm == KmeansInitPlusPlus)
        {
            kmeans_init_plusplus(initindices, collection, numclusters, distfn);
        }
        else
        {
            kmeans_init_random(initindices, collection, numclusters);
        }

        for (std::size_t i = 0; i < initindices.size(); i++) _centers[i] = collection[initindices[i]];
    }

    void run(std::size_t maxiteration, double minchangesfraction)
    {
        using namespace boost;

        // main iteration
        std::size_t iteration = 0;
        for (;;)
        {
            if (maxiteration > 0 && iteration == maxiteration) break;

            //QTime time;
            //time.start();
            std::size_t changes = 0;

            // distribute items on clusters in parallel
            thread_group pool;
            std::size_t idx = 0;
            mutex_t mtx;
            for (std::size_t i = 0; i < thread::hardware_concurrency(); i++)
            {
                pool.create_thread(bind(&kmeans::distribute_samples, this, ref(idx), ref(changes), ref(mtx)));
            }
            pool.join_all();

            iteration++;

            //std::cout << "changes: " << changes << " distribution time: " << time.elapsed() << std::endl;

            if (changes <= std::ceil(_collection.size() * minchangesfraction)) break;

            // compute new centers
            std::vector<std::size_t> clustersize(_centers.size(), 0);
            for (std::size_t i = 0; i < _collection.size(); i++)
            {
                std::size_t k = _clusters[i];

                // if it is the first assignment for that cluster, then zero the center
                if (clustersize[k] == 0) std::fill(_centers[k].begin(), _centers[k].end(), 0);

                add_operation(_centers[k], _collection[i]);
                clustersize[k]++;
            }

            // assign new centers
            std::vector<std::size_t> invalid, valid;
            for (std::size_t i = 0; i < _centers.size(); i++)
            {
                if (clustersize[i] > 0)
                {
                    div_operation(_centers[i], clustersize[i]);
                    valid.push_back(i);
                }
                else
                {
                    invalid.push_back(i);
                }
                //std::cout << "clustersize " << i << ": " << clustersize[i] << std::endl;
            }

            // fix invalid centers, i.e. those with no members
            while (!invalid.empty() && !valid.empty())
            {
                std::size_t current = invalid.back();
                std::cout << "handle invalid clusters: " << invalid.size() << std::endl;
                // compute for each valid cluster the variance
                // of distances to all members and get the
                // most distant member
                std::vector<double> maxdist(_centers.size(), 0.0);
                std::vector<double> variance(_centers.size(), 0.0);
                std::vector<std::size_t> farthest(_centers.size());
                for (std::size_t i = 0; i < valid.size(); i++)
                {
                    std::size_t c = valid[i];


                    for (std::size_t k = 0; k < _collection.size(); k++)
                    {
                        if (_clusters[k] != c) continue;

                        double d = _distfn(_collection[k], _centers[c]);
                        if (d > maxdist[c])
                        {
                            maxdist[c] = d;
                            farthest[c] = k;
                        }
                        variance[c] += d*d;
                    }
                    variance[c] /= clustersize[c];
                }

                // get cluster with highest variance and make
                // the farthest member of that cluster the
                // new center
                std::size_t c = std::distance(variance.begin(), std::max_element(variance.begin(), variance.end()));
                _centers[current] = _collection[farthest[c]];
                _clusters[farthest[c]] = current;

                std::cout << "reassign " << current << " to sample " << farthest[c] << " of cluster " << c << std::endl;

                valid.pop_back();
                invalid.pop_back();
            }

            //std::cout << "iteration " << iteration << " time: " << time.elapsed() << std::endl;
        }

        std::cout << "kmeans iterations: " << iteration << std::endl;
    }

    void run_default()
    {
        this->run(std::numeric_limits<std::size_t>::max(), 0.01);
    }

    const std::vector<std::size_t>& clusters() const
    {
        return _clusters;
    }

    const std::vector<item_t>& centers() const
    {
        return _centers;
    }

    template <class index_t>
    void make_cluster_table(std::vector<std::vector<index_t> >& table) const
    {
        table.resize(_centers.size());
        for (std::size_t i = 0; i < _clusters.size(); i++) table[_clusters[i]].push_back(i);
    }

    template <class T>
    static void add_operation(T& lhs, const T& rhs)
    {
        for (std::size_t i = 0; i < lhs.size(); i++) lhs[i] += rhs[i];
    }

    template <class T>
    static void div_operation(T& lhs, double rhs)
    {
        for (std::size_t i = 0; i < lhs.size(); i++) lhs[i] /= rhs;
    }

private:

    void distribute_samples(std::size_t& index, std::size_t& changes, mutex_t& mutex)
    {
        std::vector<double> dists(_centers.size());
        std::size_t currentchanges = 0;

        for (;;)
        {
            std::size_t i;

            {
                locker_t locker(mutex);
                if (index == _collection.size()) break;
                i = index++;
            }

            // compute distance of current point to every center
            std::transform(_centers.begin(), _centers.end(), dists.begin(), boost::bind(_distfn, boost::ref(_collection[i]), boost::arg<1>()));

            // find the minimum distance, i.e. the nearest center
            std::size_t c = std::distance(dists.begin(), std::min_element(dists.begin(), dists.end()));

            // update cluster membership
            {
                locker_t locker(_mutex);

                if (_clusters[i] != c)
                {
                    _clusters[i] = c;
                    currentchanges++;
                }
                //if (i % 1000 == 0) std::cout << "kmeans distribute: " << i << std::endl;
            }
        }

        {
            locker_t locker(mutex);
            changes += currentchanges;
        }
    }

    const collection_t& _collection;
    const dist_fn&      _distfn;

    std::vector<item_t>      _centers;
    std::vector<std::size_t> _clusters;

    boost::mutex        _mutex;
};




template <class collection_t, class dist_fn>
class kmedian
{
    typedef boost::mutex               mutex_t;
    typedef boost::lock_guard<mutex_t> locker_t;

    typedef typename collection_t::value_type item_t;

public:

    kmedian(const collection_t& collection, std::size_t numclusters, KmeansInitAlgorithm initalgorithm = KmeansInitRandom, const dist_fn& distfn = dist_fn())
        : _collection(collection), _distfn(distfn), _centers(numclusters), _clusters(collection.size())
    {
        // get initial centers
        std::vector<std::size_t> initindices;
        if (initalgorithm == KmeansInitPlusPlus)
        {
            kmeans_init_plusplus(initindices, collection, numclusters, distfn);
        }
        else
        {
            kmeans_init_random(initindices, collection, numclusters);
        }

        for (std::size_t i = 0; i < initindices.size(); i++) _centers[i] = collection[initindices[i]];
    }

    void run(std::size_t maxiteration, double minchangesfraction)
    {
        using namespace boost;

        // main iteration
        std::size_t iteration = 0;
        for (;;)
        {
            if (maxiteration > 0 && iteration == maxiteration) break;

            //QTime time;
            //time.start();
            std::size_t changes = 0;

            // distribute items on clusters in parallel
            thread_group pool;
            std::size_t idx = 0;
            mutex_t mtx;
            for (std::size_t i = 0; i < thread::hardware_concurrency(); i++)
            {
                pool.create_thread(bind(&kmedian::distribute_samples, this, ref(idx), ref(changes), ref(mtx)));
            }
            pool.join_all();

            iteration++;

            //std::cout << "changes: " << changes << " distribution time: " << time.elapsed() << std::endl;

            if (changes <= std::ceil(_collection.size() * minchangesfraction)) break;

            // compute new centers
            std::vector<std::size_t> clustersize(_centers.size(), 0);
            std::vector<std::vector<std::vector<double> > > medians(_centers.size());
            for (std::size_t i = 0; i < _collection.size(); i++)
            {
                std::size_t k = _clusters[i];

                // if it is the first assignment for that cluster, then zero the center
                if (clustersize[k] == 0)
                {
                    std::fill(_centers[k].begin(), _centers[k].end(), 0);
                    // this could be done in the initializer of medians, but
                    // there we don't have explicit item size information
                    medians[k].resize(_collection[i].size());
                }

                for (std::size_t j = 0; j < _collection[i].size(); j++)
                {
                    medians[k][j].push_back(_collection[i][j]);
                }

                clustersize[k]++;
            }

            // assign new centers
            std::vector<std::size_t> invalid, valid;
            for (std::size_t i = 0; i < _centers.size(); i++)
            {
                if (clustersize[i] > 0)
                {
                    // get the median of ith cluster by sorting each
                    // component of the assigned items separately
                    for (std::size_t k = 0; k < medians[i].size(); k++)
                    {
                        std::vector<double>& components = medians[i][k];
                        std::sort(components.begin(), components.end());
                        _centers[i][k] = components[components.size()/2];
                    }

                    valid.push_back(i);
                }
                else
                {
                    invalid.push_back(i);
                }
                //std::cout << "clustersize " << i << ": " << clustersize[i] << std::endl;
            }

            // fix invalid centers, i.e. those with no members
            while (!invalid.empty() && !valid.empty())
            {
                std::size_t current = invalid.back();

                // compute for each valid cluster the variance
                // of distances to all members and get the
                // most distant member
                std::vector<double> maxdist(_centers.size(), 0.0);
                std::vector<double> variance(_centers.size(), 0.0);
                std::vector<std::size_t> farthest(_centers.size());
                for (std::size_t i = 0; i < valid.size(); i++)
                {
                    std::size_t c = valid[i];
                    for (std::size_t k = 0; k < _collection.size(); k++)
                    {
                        if (_clusters[k] != c) continue;

                        double d = _distfn(_collection[k], _centers[c]);
                        if (d > maxdist[c])
                        {
                            maxdist[c] = d;
                            farthest[c] = k;
                        }
                        variance[c] += d*d;
                    }
                    variance[c] /= clustersize[c];
                }

                // get cluster with highest variance and make
                // the farthest member of that cluster the
                // new center
                std::size_t c = std::distance(variance.begin(), std::max_element(variance.begin(), variance.end()));
                _centers[current] = _collection[farthest[c]];
                _clusters[farthest[c]] = current;

                //std::cout << "reassign " << current << " to sample " << farthest[c] << " of cluster " << c << std::endl;

                valid.pop_back();
                invalid.pop_back();
            }

            //std::cout << "iteration " << iteration << " time: " << time.elapsed() << std::endl;
        }

        std::cout << "kmeans iterations: " << iteration << std::endl;
    }

    void run_default()
    {
        this->run(std::numeric_limits<std::size_t>::max(), 0.01);
    }

    const std::vector<std::size_t>& clusters() const
    {
        return _clusters;
    }

    const std::vector<item_t>& centers() const
    {
        return _centers;
    }

    template <class index_t>
    void make_cluster_table(std::vector<std::vector<index_t> >& table)
    {
        table.resize(_centers.size());
        for (std::size_t i = 0; i < _clusters.size(); i++) table[_clusters[i]].push_back(i);
    }

private:

    void distribute_samples(std::size_t& index, std::size_t& changes, mutex_t& mutex)
    {
        std::vector<double> dists(_centers.size());
        std::size_t currentchanges = 0;

        for (;;)
        {
            std::size_t i;

            {
                locker_t locker(mutex);
                if (index == _collection.size()) break;
                i = index++;
            }

            // compute distance of current point to every center
            std::transform(_centers.begin(), _centers.end(), dists.begin(), std::bind1st(_distfn, _collection[i]));

            // find the minimum distance, i.e. the nearest center
            std::size_t c = std::distance(dists.begin(), std::min_element(dists.begin(), dists.end()));

            // update cluster membership
            {
                locker_t locker(_mutex);

                if (_clusters[i] != c)
                {
                    _clusters[i] = c;
                    currentchanges++;
                }
                //if (i % 1000 == 0) std::cout << "kmeans distribute: " << i << std::endl;
            }
        }

        {
            locker_t locker(mutex);
            changes += currentchanges;
        }
    }

    const collection_t& _collection;
    const dist_fn&      _distfn;

    std::vector<item_t>      _centers;
    std::vector<std::size_t> _clusters;

    boost::mutex        _mutex;
};

#endif // KMEANS_HPP
