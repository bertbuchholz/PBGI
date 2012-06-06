#include <utilities/ANN_wrapper_functions.h>

void ANN_wrapper::generate_tree_from_centers(std::vector<Word> const& centers)
{
    int const dim = centers[0].size(); // dimension
    int const maxPts = centers.size();

    _dataPts = annAllocPts(maxPts, dim); // allocate data points

    for (size_t i = 0; i < centers.size(); ++i)
    {
        Word const& center = centers[i];

        for (size_t j = 0; j < center.size(); ++j)
        {
            _dataPts[i][j] = center[j];
        }

    }

    _tree = new ANNkd_tree(_dataPts, // the data points
                            maxPts, // number of points
                            dim); // dimension of spac
}

int ANN_wrapper::find_closest_center_ann(Word const& word) const
{
    int const k = 1; // number of nearest neighbors
    int const dim = word.size(); // dimension
    double const eps = 0; // error boun

    ANNpoint queryPt; // query point
    ANNidxArray nnIdx; // near neighbor indices
    ANNdistArray dists; // near neighbor distances

    queryPt = annAllocPt(dim); // allocate query point

    for (size_t i = 0; i < word.size(); ++i)
    {
        queryPt[i] = word[i];
    }

    nnIdx = new ANNidx[k]; // allocate near neigh indices
    dists = new ANNdist[k]; // allocate near neighbor dists

    _tree->annkSearch(queryPt, // query point
                      k, // number of near neighbors
                      nnIdx, // nearest neighbors (returned)
                      dists, // distance (returned)
                      eps); // error boun


    int const closest_center = nnIdx[0];

    delete [] nnIdx;
    delete [] dists;

    annDeallocPt(queryPt);

    return closest_center;
}
