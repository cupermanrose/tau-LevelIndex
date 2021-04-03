//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_LEVEL_H
#define K_LEVEL_LEVEL_H

#include "dominateG.h"
#include "kcell.h"
#include "qhull_adapter.h"
#include <iostream>
#include <fstream>

#define TEST 1

class level {
public:
    int dim, tau;
    vector<vector<kcell>> idx;
    vector<vector<float>> Allobj, OriginD;
    vector<dominateG> Grid;

public:
    level(int a_dim, int a_tau);
    ~level();

    void LoadData(char* datafile);
    void GlobalFilter(vector<int>& candidate);
    //void LocalFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void GridFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void rskyband(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void initIdx(fstream& log);
    void Build(fstream& log);
    bool VerifyDuplicate(int p, kcell& cur_cell, vector<int>& Sk, vector<kcell>& this_level);
    void CreateNewCell(int p, vector<int>& S1, vector<int>& Sk, kcell& cur_cell,kcell& newcell);
    void AddHS(int o1, int o2, bool side, vector<halfspace>& H);
    void UpdateH(kcell& cur_cell);
    void UpdateV(kcell& cur_cell);
};


#endif //K_LEVEL_LEVEL_H
