//
// Created by jiahao on 2021/3/30.
//

#include "region.h"

region::region() {
    H.clear();
    V.clear();
    innerPoint.clear();
}

region::~region() {
    H.clear();
    vector<halfspace>().swap(H);
    V.clear();
    vector<point>().swap(V);
    innerPoint.clear();
    vector<float>().swap(innerPoint);
}

void region::ToBeRoot(int dim) {
    H.clear();
    V.clear();
    //Generate vertices for the whole space
    point origin; origin.w.clear();
    for (int d = 0; d < dim - 1; d++) origin.w.push_back(0);
    V.push_back(origin);
    for (int d = 0; d < dim - 1; d++) {
        origin.w[d] = 1.0;
        V.push_back(origin);
        origin.w[d] = 0.0;
    }
    innerPoint.clear();
    for (int d = 0; d < dim - 1; d++) innerPoint.push_back(1.0 / (float)dim);
    innerPoint.push_back(0.0);
    return;
}

void region::ComputeHP(vector<float> &w, vector<float> &o1, vector<float> &o2, int dim) {
    w.clear();
    float o1_d = o1[dim - 1];
    float o2_d = o2[dim - 1];
    for (int d = 0; d < dim - 1; d++) {
        w.push_back((o1[d] - o1_d) - (o2[d] - o2_d));
    }
    w.push_back(o2_d - o1_d);
    return;
}