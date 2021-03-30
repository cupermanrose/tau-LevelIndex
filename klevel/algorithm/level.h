//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_LEVEL_H
#define K_LEVEL_LEVEL_H

#include "kcell.h"
//#include "build.h"

class level {
public:
    int dim, tau;
    vector<vector<kcell>> idx;
    vector<vector<float>> Allobj;

public:
    level(int a_dim, int a_tau);
    ~level();

    void LoadData(string datafile);
    void GlobalFilter(vector<int>& candidate, ALlobj, int tau);
    void LocalFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void initIdx();
    void Build();
    bool VerifyDuplicate(int p, kcell& cur_cell, vector<kcell>& this_level);
    void CreateNewCell(int p, vector<int>& S1, vector<int>& Sk, kcell& cur_cell,vector<kcell>& this_level);
    void UpdateH(kcell& cur_cell);
    void UpdateV(kcell& cur_cell);
};


#endif //K_LEVEL_LEVEL_H
