#include "qhull_adapter.h"

void qhull_adapter::ComputeVertex(vector<halfspace> &H, point &innerPoint, int dim) {

    return;
}

//void Halfspace::ComputeVertexByQhull(int dim, Halfspace& r) {
//    //int dim = DIM;             /* dimension of points */
//    const int numpoints = r.HPs.size() + (dim - 1) * 2 + 1;            /* number of points */
//    //coordT* halfspaces = (coordT*)malloc(dim * numpoints * sizeof(coordT)); /* array of coordinates for each point */
//    coordT* halfspaces = new coordT[dim*numpoints];
//
//    //const int numpoints = (dim - 1) * 2;            /* number of points */
//    //coordT* halfspaces = (coordT*)malloc(dim * numpoints * sizeof(coordT)); /* array of coordinates for each point */
//
//    // Ax+B<=0
//    int res = 0;
//    for (int d = 0; d < dim; d++) {
//        if (d < dim-1) {
//            halfspaces[res] = 1.0;
//            res++;
//        }
//        else {
//            halfspaces[res] = -1.0;
//            res++;
//        }
//    }
//    for (int i = 0; i < dim-1; i++) {
//        for (int d = 0; d < dim - 1; d++) {
//            if (i == d) halfspaces[res] = -1.0;
//            else halfspaces[res] = 0.0;
//            res++;
//        }
//        halfspaces[res] = 0.0; res++;
//        for (int d = 0; d < dim - 1; d++) {
//            if (i == d) halfspaces[res] = 1.0;
//            else halfspaces[res] = 0.0;
//            res++;
//        }
//        halfspaces[res] = -1.0; res++;
//    }
//    for (auto it = r.HPs.begin(); it != r.HPs.end(); it++) {
//        for (int d = 0; d < dim - 1; d++) {
//            if (it->second) halfspaces[res] = -AllHP[it->first].w[d];
//            else halfspaces[res] = AllHP[it->first].w[d];
//            res++;
//        }
//        if (it->second) halfspaces[res] = AllHP[it->first].w[dim - 1];
//        else halfspaces[res] = -AllHP[it->first].w[dim - 1];
//        res++;
//    }
//
//    //coordT* rows[TOTpoints];
//    boolT ismalloc = True;    /* True if qhull should free points in qh_freeqhull() or reallocation */
//    char flags[250];          /* option flags for qhull, see qh-quick.htm */
//    //FILE* outfile = stdout;    /* output from qh_produce_output()
//    // use NULL to skip qh_produce_output() */
//    //FILE* errfile = stderr;    /* error messages from qhull code */
//
//    int exitcode;             /* 0 if no error from qhull */
//    facetT* facet;            /* set by FORALLfacets */
//    int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
//    //int i;
//
//    qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
//    qhT* qh = &qh_qh;
//
//    QHULL_LIB_CHECK
//            qh_zero(qh, NULL);
//    fflush(NULL);
//
//    /* use qh_sethalfspace_all to transform the halfspaces yourself.
//       If so, set 'qh->feasible_point and do not use option 'Hn,...' [it would retransform the halfspaces]
//    */
//    //char qhull_cmd[] = "qhull H0 s Tcv Fp";
//    char qhull_cmd[] = "qhull H0 Fp";
//
//    coordT* feasible_point = new coordT[dim - 1];
//    for (int i = 0; i < dim - 1; i++) feasible_point[i] = r.innerPoint[i];
//    exitcode = qh_new_qhull_klevel(qh, dim, numpoints, halfspaces, ismalloc, qhull_cmd, feasible_point, NULL, NULL);
//
//    //exitcode = qh_new_qhull(qh, dim, numpoints, halfspaces, ismalloc, qhull_cmd, NULL, NULL);
//
//    r.vertices.clear();
//
//    boolT zerodiv;
//    FORALLfacets{
//            vector<float> tmp_vertex;
//            tmp_vertex.clear();
//            for (int d = 0; d < qh->hull_dim; d++) {
//                if (facet->offset < -qh->MINdenom) {
//                    tmp_vertex.push_back((facet->normal[d] / -facet->offset) + qh->feasible_point[d]);
//                }
//                else {
//                    tmp_vertex.push_back(qh_divzero(facet->normal[d], facet->offset, qh->MINdenom_1,
//                                                    &zerodiv) + qh->feasible_point[d]);
//                }
//            }
//            r.vertices.emplace_back(tmp_vertex);
//    }
//
//#ifdef qh_NOmem
//    qh_freeqhull(qh, qh_ALL);
//#else
//    qh_freeqhull(qh, !qh_ALL);
//    qh_memfreeshort(qh, &curlong, &totlong);
//    if (curlong || totlong)  /* could also check previous runs */
//        fprintf(stderr, "qhull internal warning (user_eg, #3): did not free %d bytes of long memory (%d pieces)\n",
//                totlong, curlong);
//#endif
//    //qh_freeqhull(qh, !qh_ALL);
//
//    // memory free
//
//    return;
//}