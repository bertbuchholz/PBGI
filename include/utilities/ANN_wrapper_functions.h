#ifndef ANN_WRAPPER_H
#define ANN_WRAPPER_H

#include <vector>
#include <cassert>
#include <ANN/ANN.h>

typedef std::vector<float> Word;

class ANN_wrapper
{
public:
    ANN_wrapper() : _tree(NULL)
    { }

    ~ANN_wrapper()
    {
        if (_tree)
        {
            annDeallocPts(_dataPts);
            delete _tree;
        }
    }

    void generate_tree_from_centers(std::vector<Word> const& centers);
    int find_closest_center_ann(Word const& word) const;

private:
    ANNkd_tree * _tree;
    ANNpointArray _dataPts;
};

#endif
