#include "lp_adapter.h"

void lp_adapter::lpModel(vector<halfspace>& H, lprec * lp, int dim) {
    double row[Max_Dimension];

    if (lp == NULL)
    {
        fprintf(stderr, "Unable to create new LP model\n");
        exit(0);
    }

    // constraints in each dimension //???
    for (int d = 1; d <= dim; d++)
    {
        row[0] = 0;
        for (int d_s = 1; d_s <= dim; d_s++)
        {
            if (d_s == d)
            {
                row[d_s] = 1.0;
            }
            else
            {
                row[d_s] = 0.0;
            }
        }

        add_constraint(lp, row, GE, 0.0);
        add_constraint(lp, row, LE, 1.0);
    }

    // in reduced space, sum_{w_i} should less than 1 // add it later discuss!!!
    row[0] = 0.0;
    for (int d = 1; d <= dim; d++)
    {
        row[d] = 1.0;
    }
    add_constraint(lp, row, GE, 0);
    add_constraint(lp, row, LE, 1);

    // constraints in intersected hyperplanes
    for (vector<halfspace>::iterator iter = H.begin(); iter != H.end(); iter++)
    {
        addHP(lp, dim, iter->w, iter->side);
    }

    // set scale
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

void lp_adapter::addHP(lprec* lp, int dim, vector<float>& HP, bool sideindicator) {
    double row[Max_Dimension];
    float EPS=0;
//    float EPS=0.000001;
    row[0] = 0;
    for (int d = 1; d <= dim; d++)
    {
        row[d] = HP[d - 1];
    }
    if (sideindicator == false) // o1 <= o2
    {
        add_constraint(lp, row, LE, HP[dim] - EPS);
        //add_constraint(lp, row, LE, HP[dim]);
    }
    else if (sideindicator == true) /// o1>=o2
    {
        add_constraint(lp, row, GE, HP[dim] + EPS);
        //add_constraint(lp, row, GE, HP[dim]);
    }
    else
    {
        std::cout << "Unable to detect half plane direction!!!" << endl;
    }
}

bool lp_adapter::is_Feasible(vector<halfspace>& H, vector<float>& innerPoint, int dimm) {

    float EPS=0.000001;

    double row[Max_Dimension];

    int dim = dimm - 1;

    lprec* lp = make_lp(0, dim);

    lpModel(H, lp, dim);

    //set_verbose(lp, IMPORTANT);
    set_verbose(lp, SEVERE);
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
    //set_presolve(lp, PRESOLVE_ROWDOMINATE, get_presolveloops(lp));

    double var[Max_Dimension], var1[Max_Dimension];

    for (int i = 0; i < dim + 1; i++) row[i] = 1.0;
    set_obj_fn(lp, row);
    set_maxim(lp);
    set_timeout(lp,10);
    int ret = solve(lp);
    get_variables(lp, var);

    // for reduced space

    if (ret == 0)
    {
        set_obj_fn(lp,row);
        set_minim(lp);
        set_timeout(lp,10);
        int ret = solve(lp);
        get_variables(lp, var1);
        for (int i = 0; i < dim; i++) {
            var[i] = (var[i] + var1[i])/2.0;
            if (dimm>=7){
                if (var[i]<EPS) var[i]=EPS; // for dimension 8
                if (var[i]>1.0-EPS) var[i]=1-EPS;
            }
        }
        innerPoint.clear();
        for (int i = 0; i < dim; i++) innerPoint.push_back(var[i]);
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

void lp_adapter::ComputeHP(vector<float> &o1, vector<float> &o2, vector<float> &w, int dim) {
    w.clear();
    float o1_d = o1[dim - 1];
    float o2_d = o2[dim - 1];
    for (int d = 0; d < dim - 1; d++) {
        w.push_back((o1[d] - o1_d) - (o2[d] - o2_d));
    }
    w.push_back(o2_d - o1_d);
}