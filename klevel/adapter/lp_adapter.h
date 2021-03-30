//
// Created by jiahao on 2021/3/30.
//

#ifndef K_LEVEL_LPADAPTER_H
#define K_LEVEL_LPADAPTER_H

#include <vector>
#include <iostream>
#include "lp_lib.h"

#define Max_Dimension 8
#define EPS 0.000001

using namespace std;

struct halfspace {
    vector<float> w;
    bool side;
};

struct point{
    vector<float> w;
};

class lp_adapter{

public:
    static void lpModel(vector<halfspace>& H, lprec* lp, int dim); // the lp model of this region
    static void addHP(lprec* lp, int dim, vector<float>& HP, bool sideindicator); // add a hyperplane into the lp model;
    static bool is_Feasible(vector<halfspace>& H, point& innnerPoint, int dim); // r is non-empty
    static void ComputeHP(vector<float>& o1, vector<float>& o2, vector<float>& w, int dim);
};

#endif //K_LEVEL_LPADAPTER_H




bool Halfspace::is_Feasible(Halfspace & r, int dimm) {

    double row[Max_Dimension];

    int dim = dimm - 1;

    lprec* lp = make_lp(0, dim);

    r.lpModel(lp, dim);

    /*for (unordered_map<int, bool>::iterator it = r.HPs.begin(); it != r.HPs.end(); it++) {
        r.addHP(lp, dim, it->first, it->second);
    }*/


    //set_verbose(lp, IMPORTANT);
    set_verbose(lp, SEVERE);
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
    //set_presolve(lp, PRESOLVE_ROWDOMINATE, get_presolveloops(lp));

    double var[Max_Dimension], var1[Max_Dimension];
    for (int i = 0; i < dim + 1; i++) row[i] = 0.0;
    row[1] = 1.0;
    set_maxim(lp);
    set_obj_fn(lp, row);
    set_timeout(lp,1);
    int ret = solve(lp);

    get_variables(lp, var);

    // for reduced space

    if (ret == 0)
    {
        if (dim > 1) {
            row[1] = 0.0; row[2] = 1.0;
            set_obj_fn(lp, row);
            ret = solve(lp);
            get_variables(lp, var1);
            for (int i = 0; i < dim; i++) var[i] = (var[i] + var1[i]) / 2.0;
        }
        if (dim > 2) {
            row[2] = 0.0; row[3] = 1.0;
            set_obj_fn(lp, row);
            ret = solve(lp);
            get_variables(lp, var1);
            for (int i = 0; i < dim; i++) var[i] = (var[i] + var1[i]) / 2.0;
        }
        for (int i = 0; i < dim; i++) r.innerPoint[i] = var[i];
        delete_lp(lp);
        return true;
    }
    else
    {
        if (ret==7) cout <<"lpsolver time out" << endl;
        delete_lp(lp);
        return false;
    }
}