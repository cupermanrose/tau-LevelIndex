//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_LEVEL_H
#define K_LEVEL_LEVEL_H

#include "dominateG.h"
#include "kcell.h"
#include "qhull_adapter.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <time.h>

#define TEST 1

class level {
public:
    int dim, tau;
    vector<vector<kcell>> idx;
    vector<vector<float>> Allobj, OriginD;
    unordered_map<int, dominateG> Grid;

public:
    level(int a_dim, int a_tau);
    ~level();

    void LoadData(char* datafile);
    void GlobalFilter(fstream& log, vector<int>& candidate);

    void GridFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void rskyband(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void NoFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);

    void initIdx(fstream& log);
    void Build(fstream& log);
    bool VerifyDuplicate(int p, kcell& cur_cell, vector<int>& Sk, vector<kcell>& this_level);
    void CreateNewCell(int p, vector<int>& S1, vector<int>& Sk, kcell& cur_cell,kcell& newcell);
    void AddHS(int o1, int o2, bool side, vector<halfspace>& H);
    void UpdateH(kcell& cur_cell);
    void UpdateV(kcell& cur_cell);

    void print_info(int k, clock_t& cur_time, int ave_S1, int ave_Sk, float ave_vertex, set<int>& utk_set, fstream& log);
};


#endif //K_LEVEL_LEVEL_H
