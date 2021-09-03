

#ifndef K_LEVEL_KSPR_H
#define K_LEVEL_KSPR_H

#include "algorithm/level.h"

class kspr{
public:
    static void multiple_query(level& idx, int k, int q_num, fstream& log);
    static void generate_query(level& idx, int q_num, vector<int>& q_list);
    static void single_query(level& idx, int k, int q_id, int& visit_sum, int& result_sum, fstream& log);
    static bool Find_qid_topk(kcell& this_cell, int qid);
    static bool Find_qid_Stau(kcell& this_cell, int qid);
    //static void SplitDFS(level& idx, kcell& cell, int pq, int k, int& cnt);
};



#endif //K_LEVEL_KSPR_H
