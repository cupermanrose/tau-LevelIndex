#include "qhull_adapter.h"
#include <string>
qhull_adapter::qhull_adapter() {}

qhull_adapter::~qhull_adapter() {}

void qhull_adapter::ComputeVertex2D(vector<halfspace> &H, vector<vector<float>> &V, vector<float> &innerPoint) {
    V.clear();
    float vmin=0.0,vmax=1.0;
    for (auto it=H.begin();it!=H.end();it++){
        if (it->w[0]==0) continue;
        float res=(it->w[1])/(it->w[0]);
        if (it->side) {
            if (it->w[0]>0.0) vmin=max(vmin,res);
            else vmax=min(vmax,res);
        }
        else {
            if (it->w[0]>0.0) vmax=min(vmax,res);
            else vmin=max(vmin,res);
        }
        if (vmax<vmin) cout << "Vertex Error!" << endl;
    }
    vector<float> tmp;
    tmp.clear(); tmp.push_back(vmin); V.push_back(tmp);
    tmp.clear(); tmp.push_back(vmax); V.push_back(tmp);
}

void qhull_adapter::ComputeVertex(vector<halfspace> &H, vector<vector<float>>& V, vector<float>& innerPoint, int& dim) {
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
        halfspaces[res] = 0.0-1e-4; res++;
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

    qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
    qhT* qh = &qh_qh;

    QHULL_LIB_CHECK
            qh_zero(qh, NULL);
    fflush(NULL);

    /* use qh_sethalfspace_all to transform the halfspaces yourself.
       If so, set 'qh->feasible_point and do not use option 'Hn,...' [it would retransform the halfspaces]
    */

    vector<vector<float>> ret;
    qhull_user::points_at_half_inter(ret, halfspaces, numpoints, innerPoint);
    for (auto &i1:ret) {
        vector<float> tmp;
        tmp=i1;
        V.push_back(tmp);
    }
    delete [] (halfspaces);
}

#define POINT_ID int
#define REGION std::vector<std::vector<double>>
qhull_user::qhull_user(){}

void qhull_user::points_at_half_inter(std::vector<std::vector<float>> &ret, const realT *pointCoordinates, int p_num, vector<float> &innerPoint){
    int dim=innerPoint.size();
    Qhull q;
    Coordinates feasible;
    for (float i:innerPoint) {
        feasible << i;
    }
//        cout<<innerPoint.size()<<endl;
    q.setFeasiblePoint(feasible);

//        q.setOutputStream(&cout);

    try {
        std::stringstream output;
        q.setOutputStream(&output);
        q.runQhull("normals of square", dim+1, p_num, pointCoordinates, "QJ H Pp"); // halfspace intersect

        q.outputQhull("Fp");
        int UNUSED;
        output>>UNUSED;
        int num_out;
        output>>num_out;
        for (int i = 0; i < num_out; ++i){
            std::vector<float> vs(dim);
            for (int j = 0; j < dim; ++j){
                output>>vs[j];
            }
            ret.push_back(vs);
        }
    } catch (exception &e) {// catch by ref

    }
}

void qhull_user::get_neiVT_of_VT(Qhull &q, const std::vector<int>&pd_ids, std::unordered_map<int, std::vector<int>> &ret){
    // given a solved Qhull obj q with pd_ids, return neigh vertex of id in pd_ids
    // return pdt_id
    auto pt_neiF = get_neiFacets_of_points(q, pd_ids);  // require index
    auto f_pt = get_points_of_facets(q, pd_ids); // require pdt_id
    for (int i = 0; i < pd_ids.size(); ++i) {
        if(!pt_neiF[i].empty()){
            std::unordered_set<int> nei_f_s;
            for (int j = 0; j < pt_neiF[i].size(); ++j) {
                int f=pt_neiF[i][j];
                for (int k = 0; k < f_pt[f].size(); ++k) {
                    nei_f_s.insert(f_pt[f][k]);
                }
            }
            nei_f_s.erase(pd_ids[i]);
            ret[pd_ids[i]]=std::vector<int>(nei_f_s.begin(), nei_f_s.end());
        }
    }
}

std::unordered_map<int, std::vector<int>> qhull_user::get_neiVT_of_VT(Qhull &q, const std::vector<int>&pd_ids){
    // return pdt_id
    std::unordered_map<int, std::vector<int>> ret;
    get_neiVT_of_VT(q, pd_ids, ret);
    return ret;
}

std::vector<std::vector<int>> qhull_user::get_points_of_facets(Qhull &q, const std::vector<int> &pd_ids){
    // return each facet's options
    // return in form of pdt_id
    std::vector<std::vector<int>> ret;
    std::stringstream output;
    q.setOutputStream(&output);
    q.outputQhull("Fv");
    int num_f;
    output>>num_f;
    int num_v;
    int tmp;
    for (int i = 0; i < num_f; ++i) {
        output>>num_v;
        std::vector<int> vs;
        for (int j = 0; j <num_v ; ++j) {
            output>>tmp;
            if(tmp<pd_ids.size()){
                vs.push_back(pd_ids[tmp]);
            }
        }
        ret.push_back(vs);
    }
    return ret;
}

std::vector<int> qhull_user::get_CH_pointID(Qhull &q, const std::vector<int> &pd_ids){
    // return the upper ch made by pd_ids
    std::stringstream output;
    q.setOutputStream(&output);
    q.outputQhull("Fx");
    std::vector<int> CH;
    int num;
    output>>num;
    int tmp;
    for (int i = 0; i < num; ++i) {
        output>>tmp;
        if(tmp<pd_ids.size()){
            CH.push_back(pd_ids[tmp]);
        }
    }
    return CH;
}

std::unordered_map<POINT_ID, REGION> qhull_user::get_neiFacetsNorm_of_point(Qhull &q, const std::vector<int> &pd_ids){
    std::unordered_map<POINT_ID, REGION> pt_r;
    get_neiFacetsNorm_of_point(q, pd_ids, pt_r);
    return pt_r;
}

void qhull_user::get_neiFacetsNorm_of_point(Qhull &q, const std::vector<int> &pd_ids, std::unordered_map<POINT_ID, REGION> &ret){
    // get the neighbor facets' norm for each point in pd_ids
    // in form of index
    std::vector<std::vector<double>> facets_norms=get_norm_of_all_facets(q);
    std::vector<std::vector<int>> pt_neif=get_neiFacets_of_points(q, pd_ids);// in form with index
    assert(pt_neif.size()>=pd_ids.size());
    for (int i = 0; i <pd_ids.size() ; ++i) {
        if(pt_neif[i].empty()){ // not a vertex

        }else{ // is a vertex
            std::vector<std::vector<double>> cone;
            for (int j = 0; j < pt_neif[i].size(); ++j) {
                cone.push_back(facets_norms[pt_neif[i][j]]);
            }
            ret[pd_ids[i]]=cone;
        }
    }
}

std::vector<std::vector<double>> qhull_user::get_norm_of_all_facets(Qhull &q){
    // get the norm of each facet, norm is pointing out ch
    std::stringstream outer_normstr;
    q.setOutputStream(&outer_normstr);
    q.outputQhull("n");
    int dim_p1;
    outer_normstr >> dim_p1;
    int dim = dim_p1-1;
    int num_facets;
    outer_normstr >> num_facets;
    std::vector<std::vector<double>> cone;
    double offset;
    for (int i = 0; i < num_facets; ++i) {
        std::vector<double> n(dim);
        for (int j = 0; j < dim; ++j) {
            outer_normstr >> n[j];
        }
        outer_normstr >> offset; // offset
        cone.push_back(n);
    }
    return cone;
}

std::vector<std::vector<double>> qhull_user::get_cone_norms(Qhull &q, std::vector<std::vector<double>> &points){
    // trans the cone defined by points(norms) into halfspace form
    // make sure the first of Qhull input is \vec{0}_{d}
    std::stringstream opt_neibor_facets;
    q.setOutputStream(&opt_neibor_facets);
    q.outputQhull("FN");
    int num_point;
    int UNUSED;
    int USED;
    int f_cnt;
    opt_neibor_facets >> num_point;
    opt_neibor_facets>>f_cnt;
    std::vector<int> cone_facets(f_cnt);
    for (int j = 0; j < f_cnt; ++j) {
        opt_neibor_facets>>cone_facets[j];
    }
    std::vector<std::vector<double>> norm_of_facets=get_norm_of_all_facets(q);
    std::vector<std::vector<double>> norm_of_cone;
    for(int f_id: cone_facets){
        if(f_id<0){
            continue;
        }
//            assert(f_id<norm_of_facets.size());
        norm_of_cone.push_back(norm_of_facets[f_id]);
    }
//        assert(!norm_of_cone.empty());
    return norm_of_cone;
}

std::vector<std::vector<int>> qhull_user::get_neiFacets_of_points(Qhull &q, const std::vector<int> &pd_ids){
        // in the order of index
        // get neighboring facets for each point
        std::stringstream opt_neibor_facets;
        q.setOutputStream(&opt_neibor_facets);
        q.outputQhull("FN");
        int num_point;
        int UNUSED;
        int USED;
        opt_neibor_facets >> num_point;
        std::vector<std::vector<int>> p_nei_fs;
        for (int i = 0; i < num_point; ++i) {
            if(i>=pd_ids.size()){
                break;
            }
            int f_cnt;
            opt_neibor_facets>>f_cnt;
            std::vector<int> fs;
            if(f_cnt<=1){
                // interior point or coplanar or vertices belong to not good facet
                if(f_cnt==1){
                    opt_neibor_facets>>UNUSED;
                }
            }else{
                for (int j = 0; j < f_cnt; ++j) {
                    opt_neibor_facets>>USED;
                    if(USED<0){
                        continue;
                    }
                    fs.push_back(USED);
                }
            }
            p_nei_fs.push_back(fs);
        }
        return p_nei_fs;
    }

