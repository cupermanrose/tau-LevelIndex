
#ifndef K_LEVEL_UTILS_H
#define K_LEVEL_UTILS_H
#include <vector>
#include "rtree.h"
void kskyband(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k,
        bool rtree=true,
        const Rtree *rtree_rt= nullptr);

void kskyband_rtree(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k,
        const Rtree *rtree_rt,
        std::unordered_map<long int, RtreeNode*> &ramTree);

void kskyband_nortree(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k);

#endif //K_LEVEL_UTILS_H


