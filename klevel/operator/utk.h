
#ifndef K_LEVEL_UTK_H
#define K_LEVEL_UTK_H

#include "algorithm/level.h"
#include "k_level_lib.h"

class utk {
public:
    static bool isIn(vector<float>& v, vector<float>& Qregion, int dim);
    static bool isIn(vector<float>& v, vector<halfspace>& H, int dim);
    static bool Intersect(vector<float>& Qregion,region& r, int dim);
    static void AddQregion(vector<float>& Qregion, region& r, int dim);
    static void multiple_query(level& idx, int k, int q_num, float utk_side_length, fstream& log);
    static void generate_query(level& idx, int q_num, float utk_side_length, vector<vector<float>>& q_list);
    static void single_query(level& idx, int k, vector<float>& Qregion, int& visit_sum, int& result_sum, fstream& log);


    /*static int single_query(level& idx, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree,
                            int k, vector<float>& Qregion, fstream& log);
    static int single_query_largek(level& idx, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree,
                                   int k, vector<float>& Qregion, fstream& log);
    static void SplitDFS(level& idx, kcell& cell, vector<float>& Qregion, int k, unordered_set<int>& results);*/
};


#endif //K_LEVEL_UTK_H
