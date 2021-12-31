//
// Created by 12859 on 2021/12/24.
//

#ifndef K_LEVEL_TOPK_H
#define K_LEVEL_TOPK_H
#include "algorithm/level.h"

class dtopk{
public:
    static void multiple_query(level& idx, int k, int q_num, fstream& log);
    static void generate_query(level& idx, int q_num, vector<vector<float>>& ws_ret);
    static void single_query(level& idx, int k, vector<float> &w, int& visit_sum, int& result_sum, fstream& log);
};


#endif //K_LEVEL_TOPK_H
