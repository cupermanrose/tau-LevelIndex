//
// Created by jiahao on 2021/3/30.
//

#ifndef K_LEVEL_REGION_H
#define K_LEVEL_REGION_H

#include <vector>
#include "lp_adapter.h"
#include "qhull_adapter.h"

using namespace std;

class region {
public:

    vector<halfspace> H; // halfspace representation
    vector<point> V; // vertice representation
    vector<float> innerPoint;

public:
    region();
    ~region();

    void ToBeRoot(int dim);
    void ComputeHP(vector<float>& w, vector<float>& o1, vector<float>& o2, int dim);
};


#endif //K_LEVEL_REGION_H
