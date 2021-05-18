
#ifndef K_LEVEL_UTILS_H
#define K_LEVEL_UTILS_H
#include <vector>
#include <string>
#include "rtree.h"
using std::unordered_set;
using std::unordered_map;
using std::vector;

void kskyband(
        vector<int> &ret,
        const vector<vector<float>> &data,
        const int k,
        bool rtree=true,
        const Rtree *rtree_rt= nullptr);

void kskyband_rtree(
        vector<int> &ret,
        const vector<vector<float>> &data,
        const int k,
        const Rtree *rtree_rt,
        unordered_map<long int, RtreeNode*> &ramTree);

void kskyband_nortree(
        vector<int> &ret,
        const vector<vector<float>> &data,
        const int k);

void onionlayer(vector<int> &ret, vector<int>& ret_layer, vector<int>& candidate, vector<vector<float>>& data, int k);

float GetScore(vector<float> &w, vector<float> &p, int dim); // w[dim-1]=1.0-sigma(w[0] to w[dim-2])


// Return true: oj dominates oi in this region;
template<typename POINT>
bool RegionDominate(vector<POINT> &V, vector<float> &oi, vector<float> &oj, int dim) {
    for (auto it=V.begin();it!=V.end();it++){
        if (GetScore(*it,oi,dim)>GetScore(*it,oj,dim)) return false;
    }
    return true;
}

class ch{
    unordered_set<int> rest;
    unordered_map<int, int> pdtid_layer;
    unordered_map<int, int> dominated_cnt;// only use for build k-convex-hull

    vector<vector<int>> chs;

    vector<int> EMPTY;
    vector<vector<float>>& pointSet;
    int d;
public:
    unordered_set<int> all;
    vector<int> rskyband;
    unordered_map<int, vector<int>> A_p;
    unordered_map<int, unordered_set<int>> do_map;
    unordered_map<int, unordered_set<int>> dominated_map;
    ch(vector<int> &idxes, vector<vector<float>> &point_set, int dim);

    void fast_non_dominate_sort(
            const unordered_map<int, unordered_set<int>> &do_map,
            unordered_map<int, int>& dominated_cnt,
            const vector<int> &last_layer);

    void build_do_re(vector<int> &idxes, vector<vector<float>> &point_set, int dim);

    const vector<int>& get_next_layer();

    int get_option_layer(int option);

    const vector<int>& get_neighbor_vertex(int option);

    const vector<int>& get_layer(int which_layer);

    ~ch();
};

void build_onion();

void read_onion(const std::string &filename, vector<vector<int>> &ret);


#endif //K_LEVEL_UTILS_H


