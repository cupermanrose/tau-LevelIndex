//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_KCELL_H
#define K_LEVEL_KCELL_H

#include "region.h"
#include <set>

class kcell {
public:
    int curk, objID;
    region r;
    set<int> topk; // the top-kth in this region
    set<int> Stau; // top-tau candidates set

public:
    kcell();
    ~kcell();

    void TobeRoot(vector<int>& candidates, int dim);
};


#endif //K_LEVEL_KCELL_H
