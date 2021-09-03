
#ifndef K_LEVEL_LPADAPTER_H
#define K_LEVEL_LPADAPTER_H

#include <vector>
#include <iostream>
#include "lp_lib.h"

#define Max_Dimension 8
using namespace std;

struct halfspace {
    vector<float> w;
    bool side;
};

class lp_adapter{

public:
    static void lpModel(vector<halfspace>& H, lprec* lp, int dim); // the lp model of this region
    static void addHP(lprec* lp, int dim, vector<float>& HP, bool sideindicator); // add a hyperplane into the lp model;
    static bool is_Feasible(vector<halfspace>& H, vector<float>& innnerPoint, int dim); // r is non-empty
    static void ComputeHP(vector<float>& o1, vector<float>& o2, vector<float>& w, int dim);

};

#endif //K_LEVEL_LPADAPTER_H




