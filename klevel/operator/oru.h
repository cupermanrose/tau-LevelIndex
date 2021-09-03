

#ifndef K_LEVEL_ORU_H
#define K_LEVEL_ORU_H

#include "algorithm/level.h"
#include "adapter/qp_adapter.h"
#include "k_level_lib.h"

class oru {
public:
    //
    static void generate_query(level& idx, int q_num,  vector<vector<float>>& q_list);
    static bool isIn(vector<float>& v, vector<halfspace>& H, int dim);
    static float GetDistance(vector<float>& q, region& r, int dim);
    static float kcell_filter(vector<kcell>& L, vector<bool>& filter, int ret_size, vector<float>& q, int dim);
    static float single_query(level& idx, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree,
                              int k, int ret_size, vector<float>& q, fstream& log);
    //static float single_query_largek(level& idx, int k, int ret_size, vector<float>& q, fstream& log);
    static float single_query(level &idx, int k, int ret_size, vector<float>& q, fstream &log, int& push_cnt, int& pop_cnt);
    static void multiple_query(level &idx, int k, int ret_size, int q_num, fstream &log);

};


#endif //K_LEVEL_ORU_H
