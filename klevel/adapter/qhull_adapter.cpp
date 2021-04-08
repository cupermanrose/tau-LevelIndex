#include "qhull_adapter.h"

qhull_adapter::qhull_adapter() {}

qhull_adapter::~qhull_adapter() {}

void qhull_adapter::ComputeVertex(vector<halfspace> &H, vector<point>& V, vector<float>& innerPoint, int& dim) {
    const int numpoints = H.size() + (dim - 1) * 2 + 1;            /* number of points */
    //coordT* halfspaces = (coordT*)malloc(dim * numpoints * sizeof(coordT)); /* array of coordinates for each point */
    coordT* halfspaces = new coordT[dim*numpoints];

    // Ax+B<=0
    int res = 0;
    for (int d = 0; d < dim; d++) {
        if (d < dim-1) {
            halfspaces[res] = 1.0;
            res++;
        }
        else {
            halfspaces[res] = -1.0;
            res++;
        }
    }
    for (int i = 0; i < dim-1; i++) {
        for (int d = 0; d < dim - 1; d++) {
            if (i == d) halfspaces[res] = -1.0;
            else halfspaces[res] = 0.0;
            res++;
        }
        halfspaces[res] = 0.0; res++;
        for (int d = 0; d < dim - 1; d++) {
            if (i == d) halfspaces[res] = 1.0;
            else halfspaces[res] = 0.0;
            res++;
        }
        halfspaces[res] = -1.0; res++;
    }
    for (auto it = H.begin(); it != H.end(); it++) {
        for (int d = 0; d < dim - 1; d++) {
            if (it->side) halfspaces[res] = -it->w[d];
            else halfspaces[res] = it->w[d];
            res++;
        }
        if (it->side) halfspaces[res] = it->w[dim - 1];
        else halfspaces[res] = -it->w[dim - 1];
        res++;
    }

    //coordT* rows[TOTpoints];
    boolT ismalloc = True;    /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[250];          /* option flags for qhull, see qh-quick.htm */
    //FILE* outfile = stdout;    /* output from qh_produce_output()
    // use NULL to skip qh_produce_output() */
    //FILE* errfile = stderr;    /* error messages from qhull code */

    int exitcode;             /* 0 if no error from qhull */
    facetT* facet;            /* set by FORALLfacets */
    int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
    //int i;

    qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
    qhT* qh = &qh_qh;

    QHULL_LIB_CHECK
            qh_zero(qh, NULL);
    fflush(NULL);

    /* use qh_sethalfspace_all to transform the halfspaces yourself.
       If so, set 'qh->feasible_point and do not use option 'Hn,...' [it would retransform the halfspaces]
    */
    //char qhull_cmd[] = "qhull H0 s Tcv Fp";
    //char qhull_cmd[] = "qhull H0 Fp";

    string s="qhull H";
    for (int i = 0; i < dim - 2; i++) {
        s += to_string(innerPoint[i]) + ",";
    }
    s += to_string(innerPoint[dim-2])+" Fp";
    char qhull_cmd[s.length()+1];
    strcpy(qhull_cmd, s.c_str());

    //coordT* feasible_point = new coordT[dim - 1];
    /*for (int i = 0; i < dim - 1; i++) {
        feasible_point[i] = innerPoint[i];
    }*/

    //exitcode = qh_new_qhull_klevel(qh, dim, numpoints, halfspaces, ismalloc, qhull_cmd, feasible_point, NULL, NULL);

    exitcode = qh_new_qhull(qh, dim, numpoints, halfspaces, ismalloc, qhull_cmd, NULL, NULL);

    V.clear();

    boolT zerodiv;
    FORALLfacets{
            point tmp; tmp.w.clear();
            for (int d = 0; d < qh->hull_dim; d++) {
                if (facet->offset < -qh->MINdenom) {
                    tmp.w.push_back((facet->normal[d] / -facet->offset) + qh->feasible_point[d]);
                }
                else {
                    tmp.w.push_back(qh_divzero(facet->normal[d], facet->offset, qh->MINdenom_1,
                                                    &zerodiv) + qh->feasible_point[d]);
                }
            }
            V.emplace_back(tmp);
    }

#ifdef qh_NOmem
    qh_freeqhull(qh, qh_ALL);
#else
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
    if (curlong || totlong)  /* could also check previous runs */
        fprintf(stderr, "qhull internal warning (user_eg, #3): did not free %d bytes of long memory (%d pieces)\n",
                totlong, curlong);
#endif
    //qh_freeqhull(qh, !qh_ALL);

    // memory free

    return;
}
