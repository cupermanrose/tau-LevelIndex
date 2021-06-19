//
// Created by cuper on 2021/4/16.
//

#ifndef K_LEVEL_KSPR_H
#define K_LEVEL_KSPR_H

#include "algorithm/level.h"

class kspr{
public:
    static void multiple_query(level& idx, int k, int q_num, fstream& log);
    static void generate_query(level& idx, int q_num, vector<int>& q_list);
    static int single_query(level& idx, int k, int q_id, fstream& log);
    static bool Find_qid_topk(kcell& this_cell, int qid);
    static bool Find_qid_Stau(kcell& this_cell, int qid);
    static int single_query_largek(level& idx, int k, int q_id, fstream& log);
    //static void BuildInvertedIdx(vector<vector<int>>& inverted_idx, level& idx, vector<kcell>& L);
};



#endif //K_LEVEL_KSPR_H
