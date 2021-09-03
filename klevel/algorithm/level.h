

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
    int dim, tau, ik;
    vector<vector<kcell>> idx;
    unordered_map<size_t,int> region_map;
    vector<vector<float>> Allobj, OriginD;
    vector<int> global_layer;
    unordered_map<int, int> levelId_2_dataId;
    //unordered_map<int, dominateG> Grid;

public:
    level(int a_dim, int a_tau, int a_ik);
    ~level();

    // common functions
    /////////////////////////////////////////////////
    void LoadData(string datafile);

    void GlobalFilter(fstream& log, vector<int>& candidate);

    void FreeMem(int k);
    void profiling(int k, clock_t& level_zero_time, double& rskyband_time, double& verify_time, double& isFeasible_time,double& updateV_time, fstream& log);
    void print_info(int k, int cellsum, int valid_cell, clock_t& level_zero_time, clock_t& level_k_time,
                    int& ave_S1, int& ave_Sk, int& ave_vertex, int& ave_next, int& suc_split, int& HS_size,
                    set<int>& utk_set, fstream& log);
    void print_system_info(fstream& log);
    void WriteToDisk(int k, ofstream& idxout);
    void ReadFromDisk(int k, ifstream& idxin);
    int EdgeComputation(int k);

    // advanced building
    /////////////////////////////////////////////////
    void initIdx(fstream& log);
    void Build(fstream& log, ofstream& idxout);
    void Build_nofilter(fstream& log, ofstream& idxout);
    bool VerifyDuplicate(kcell& newcell, vector<kcell>& this_level); // hash version
    void CreateNewCell(int p, vector<int>& S1, vector<int>& Sk, kcell& cur_cell,kcell& newcell);
    void AddHS(int o1, int o2, bool side, vector<halfspace>& H);
    void UpdateH(kcell& cur_cell);
    void UpdateH_S1(kcell &cur_cell, vector<int> &S1);
    void UpdateV(kcell& cur_cell, int& ave_vertex);
    void LocalFilter(int k, vector<int>& S1, vector<int>& Sk, kcell& cur_cell, int& ave_S1, int& ave_Sk);
    void rskyband(vector<int>& S1, vector<int>& Sk, kcell& cur_cell, int k);
    void NoFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);

    //void GridFilter(vector<int>& S1, vector<int>& Sk, kcell& cur_cell);
    //bool VerifyDuplicate(int p, kcell& cur_cell, vector<int>& Sk, vector<kcell>& this_level); // for-loop version

    // Insertion-based building
    void MergeCell(vector<kcell>& L_NoMerge, vector<kcell>& L_Merge);
    void IncBuild(fstream& log, ofstream& idxout);
    void SplitCell(int p, int i, vector<kcell>& L);

    //void InsBuild(fstream& log, ofstream& idxout);
    //void SplitCell_Ins(int p, kcell& cur_cell, vector<kcell>& L);


    // BulkLoad DFS
    //void DFSBuild(fstream& log, ofstream& idxout);
    //void SplitDFS(kcell& cell , vector<kcell>& L, ofstream& idxout, int& kcell_num);

    void SingleCellSplit(int qk, kcell& cur_cell, vector<kcell>& this_level);
};


#endif //K_LEVEL_LEVEL_H
