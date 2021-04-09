//
// Created by jiahao on 2021/3/30.
//

#ifndef K_LEVEL_QHULLADAPTER_H
#define K_LEVEL_QHULLADAPTER_H

#include "lp_adapter.h"
//#include "libqhull_r/libqhull_r.h"
#include "libqhull_r/qhull_ra.h"

class qhull_adapter{
public:
    qhull_adapter();
    ~qhull_adapter();
    static void ComputeVertex(vector<halfspace>& H, vector<point>& V, vector<float>& innerPoint, int& dim);
};

#endif //K_LEVEL_QHULLADAPTER_H