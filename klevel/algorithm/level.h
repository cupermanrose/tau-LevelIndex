//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_LEVEL_H
#define K_LEVEL_LEVEL_H

//#include "dominateG.h"
#include "kcell.h"
#include "qhull_adapter.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <time.h>


//for system info
#include <unistd.h>

#define TEST 1

class level {
public:
    int dim, tau;
    deque<deque<kcell>> idx;
    unordered_map<size_t,int> region_map;
    vector<vector<float>> Allobj, OriginD;
    vector<int> global_layer;
    //unordered_map<int, dominateG> Grid;

public:
    level(int a_dim, int a_tau);
    ~level();

    void LoadData(char* datafile);
    void GlobalFilter(fstream& log, vector<int>& candidate);
    void FreeMem(int k);

    void LocalFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell, int& ave_S1, int& ave_Sk);
    //void GridFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void rskyband(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    void NoFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);

    void initIdx(fstream& log);
    void Build(fstream& log, ofstream& idxout);
    bool VerifyDuplicate(kcell& newcell, deque<kcell>& this_level); // hash version
    //bool VerifyDuplicate(int p, kcell& cur_cell, vector<int>& Sk, vector<kcell>& this_level); // for-loop version
    void CreateNewCell(int p, vector<int>& S1, vector<int>& Sk, kcell& cur_cell,kcell& newcell);
    void AddHS(int o1, int o2, bool side, vector<halfspace>& H);
    void UpdateH(kcell& cur_cell);
    void UpdateV(kcell& cur_cell, int& ave_vertex);

    void profiling(int k, clock_t& level_zero_time, double& rskyband_time, double& verify_time, double& isFeasible_time,double& updateV_time, fstream& log);
    void print_info(int k, clock_t& level_zero_time, clock_t& level_k_time, int& ave_S1, int& ave_Sk, int& ave_vertex, set<int>& utk_set, fstream& log);
    void print_system_info(fstream& log);

    void WriteToDisk(int k, ofstream& idxout);
    void ReadFromDisk(int k, ifstream& idxin);
};


#endif //K_LEVEL_LEVEL_H
