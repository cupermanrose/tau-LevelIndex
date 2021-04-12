
#ifndef IPREF_QHULL_USER_H
#define IPREF_QHULL_USER_H
#include "RboxPoints.h"
#include "QhullError.h"
#include "QhullQh.h"
#include "QhullFacet.h"
#include "QhullFacetList.h"
#include "QhullFacetSet.h"
#include "QhullLinkedList.h"
#include "QhullPoint.h"
#include "QhullUser.h"
#include "QhullVertex.h"
#include "QhullVertexSet.h"
#include "Qhull.h"
#include <iostream>
#include <iomanip> // setw
#include <unordered_map>
#include <unordered_set>
//TODO merge to qhull_adapter
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::cout;
using std::endl;
using std::string;
using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetSet;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexSet;
using orgQhull::RboxPoints;
using orgQhull::Coordinates;

class qhull_user{

#define POINT_ID int
#define REGION std::vector<std::vector<double>>
public:
    qhull_user(){
    }
    static void points_at_half_inter(std::vector<std::vector<float>> &ret, const realT *pointCoordinates, int p_num, vector<float> &innerPoint){
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
        } catch (exception) {

        }
    }

    void get_neiVT_of_VT(Qhull &q, const std::vector<int>&pd_ids, std::unordered_map<int, std::vector<int>> &ret){
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

    std::unordered_map<int, std::vector<int>> get_neiVT_of_VT(Qhull &q, const std::vector<int>&pd_ids){
        // return pdt_id
        std::unordered_map<int, std::vector<int>> ret;
        get_neiVT_of_VT(q, pd_ids, ret);
        return ret;
    }

    std::vector<std::vector<int>> get_points_of_facets(Qhull &q, const std::vector<int> &pd_ids){
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

    std::vector<int> get_CH_pointID(Qhull &q, const std::vector<int> &pd_ids){
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

    std::unordered_map<POINT_ID, REGION> get_neiFacetsNorm_of_point(Qhull &q, const std::vector<int> &pd_ids){
        std::unordered_map<POINT_ID, REGION> pt_r;
        get_neiFacetsNorm_of_point(q, pd_ids, pt_r);
        return pt_r;
    }

    void get_neiFacetsNorm_of_point(Qhull &q, const std::vector<int> &pd_ids, std::unordered_map<POINT_ID, REGION> &ret){
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

    std::vector<std::vector<double>> get_norm_of_all_facets(Qhull &q){
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

    std::vector<std::vector<double>> get_cone_norms(Qhull &q, std::vector<std::vector<double>> &points){
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
//        for (int i = 0; i < num_point-1; ++i) {
//            int f_cnt;
//            opt_neibor_facets>>f_cnt;
//            if(f_cnt<=1){
//                // interior point or coplanar or vertices belong to not good facet
//                if(f_cnt==1){
//                    opt_neibor_facets>>UNUSED;
//                }
//            }else{
//                for (int j = 0; j < f_cnt; ++j) {
//                    opt_neibor_facets>>UNUSED;
//                }
//            }
//        }
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

    std::vector<std::vector<int>> get_neiFacets_of_points(Qhull &q, const std::vector<int> &pd_ids){
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
};

#endif //IPREF_QHULL_USER_H
