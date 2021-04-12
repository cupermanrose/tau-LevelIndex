//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_KCELL_H
#define K_LEVEL_KCELL_H

#include "region.h"
#include <algorithm>
#include <unordered_set>
#include <boost/functional/hash.hpp>

class kcell {
public:
    int curk, objID;
    size_t hash_value;
    region r;
    unordered_set<int> topk; // the top-kth in this region
    unordered_set<int> Stau; // top-tau candidates set

public:
    kcell();
    ~kcell();

    void Get_HashValue();
    void TobeRoot(vector<int>& candidates, int dim);
};


#endif //K_LEVEL_KCELL_H
