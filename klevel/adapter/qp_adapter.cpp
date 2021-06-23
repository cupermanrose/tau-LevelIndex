#include <vector>

#include "osqp.h"
#include <cassert>
#include <chrono>
#include  <numeric>
#include <cmath>
using namespace std;


double getDistance(const vector<float>& w, const vector<vector<float>>& H){
    /*
     * if   w is dim
     * then H is dim+1
     * minimize (x-w)^T \cdot (x-w) subject to H_{i, 1:dim} \cdot x <= H_{i, dim+1}
     */
    /*
     * csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p)
     * in OSQP,
     * QP problem's "P" and "A" are represented as csc format sparse matrix
     * for more detail about transform dense_matrix into csc_matrix
     * see https://en.wikipedia.org/wiki/Sparse_matrix
     *
     * minimize 0.5 x^T P x + q^T x subject to l <= Ax <= u
     *
     * The OSQP example can be seen from here:
     * https://osqp.org/docs/examples/setup-and-solve.html
     */
    /*
     * the detail document can be reached from /doc/For k-level Distance from a point to a convex hull and QP problem.docx
     */
    assert(!w.empty());
    int dim=w.size();
    double w_L2=0;
    for(float i:w){
        w_L2+=i*i;
    }
    OSQPWorkspace *work;
    OSQPSettings  *settings;
    OSQPData      *data;
    vector<c_float> P_x;
    vector<c_int> P_i;
    vector<c_int> P_p;
    vector<c_float> A_x;
    vector<c_int> A_i;
    vector<c_int> A_p;
    vector<c_float> q;
    vector<c_float> l;
    vector<c_float> u;

    P_x=vector<c_float>(dim, 1.0);
    c_int P_nnz = dim;
    P_i=vector<c_int>(dim);
    iota(P_i.begin(), P_i.end(), 0); //{0, 1, 2, ..., dim-1}
    P_p=vector<c_int>(dim+1);
    iota(P_p.begin(), P_p.end(), 0); //{0, 1, 2, ..., dim-1}

    q=vector<c_float>(w.begin(), w.end());
    for (c_float &qi:q) {
        qi=-qi;
    }

    A_x=vector<c_float>((2+H.size())*dim);
    c_int A_x_idx=0, h_idx=0;
    while(A_x_idx<A_x.size()){
        A_x[A_x_idx++]=1.0;
        for (int i = 0; i < H.size(); ++i) {
            A_x[A_x_idx++]=H[i][h_idx];
        }
        h_idx++;
        A_x[A_x_idx++]=1.0;
    }
    c_int A_nnz = A_x.size();
    A_i=vector<c_int>((2+H.size())*dim);
    c_int A_i_idx=0, h_iidx=0;
    while(A_i_idx<A_i.size()){
        A_i[A_i_idx++]=0;
        for (int i = 0; i < H.size(); ++i) {
            A_i[A_i_idx++]=i+1;
        }
        A_i[A_i_idx++]=h_iidx+1+H.size();
        ++h_iidx;
    }
    A_p=vector<c_int>(dim+1);
    c_int A_p_idx=0, h_pidx=0;
    while(A_p_idx+1<A_p.size()){
        A_p[A_p_idx+1]=A_p[A_p_idx]+2+H.size();
        ++A_p_idx;
    }

    l=vector<c_float>(dim+1+H.size(), 0);
    l[0]=-INFINITY;
    for (int j = 0; j <H.size() ; ++j) {
        l[j+1]=-INFINITY;
    }
    u=vector<c_float>(dim+1+H.size(), INFINITY);
    u[0]=1;
    for (int j = 0; j <H.size() ; ++j) {
        u[j+1]=H[j][dim];
    }

    settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    data     = (OSQPData *)c_malloc(sizeof(OSQPData));
    data->n = dim;
    data->m = dim+1+H.size(); //constrain number
    data->P = csc_matrix(data->n, data->n, P_nnz, P_x.data(), P_i.data(), P_p.data());
    data->q = q.data();
    data->A = csc_matrix(data->m, data->n, A_nnz, A_x.data(), A_i.data(), A_p.data());
    data->l = l.data();
    data->u = u.data();

    double ret=INFINITY;
    if (settings) {
        osqp_set_default_settings(settings);
        settings->verbose=false; // keep quiet
        settings->alpha = 1.0; // Change alpha parameter
        osqp_setup(&work, data, settings);

        int flag = osqp_solve(work);
        if(flag!=0){
            ret= INFINITY;
        }
        ret=2*work->info->obj_val+w_L2;
        ret>=0?sqrt(ret):0;
    }

    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings) {
        c_free(settings);
    }

    return ret;
}

