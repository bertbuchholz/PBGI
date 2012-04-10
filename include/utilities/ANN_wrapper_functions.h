#ifndef ANN_WRAPPER_H
#define ANN_WRAPPER_H

#include <vector>
#include <ANN/ANN.h>

typedef std::vector<float> Word;

ANNkd_tree * generate_kd_tree_from_centers(std::vector<Word> const& centers);
int find_closest_center_ann(Word const& word, ANNkd_tree* kdTree);

#endif
