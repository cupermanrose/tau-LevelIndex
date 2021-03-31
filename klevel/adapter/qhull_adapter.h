//
// Created by jiahao on 2021/3/30.
//

#ifndef K_LEVEL_QHULLADAPTER_H
#define K_LEVEL_QHULLADAPTER_H

#include "lp_adapter.h"

class qhull_adapter{

public:
    static void ComputeVertex(vector<halfspace>& H, point& innerPoint, int dim);
};

#endif //K_LEVEL_QHULLADAPTER_H