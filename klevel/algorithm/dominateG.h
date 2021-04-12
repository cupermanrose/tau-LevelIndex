//
// Created by cuper on 2021/4/1.
//

#ifndef K_LEVEL_DOMINATEG_H
#define K_LEVEL_DOMINATEG_H

#include <iostream>
#include <set>
#include <unordered_map>
#include <algorithm>
#include "region.h"

using namespace std;

#define div_num 10
#define EPS 0.000001

class dominateG {
public:
    int Gid;
    unordered_map<int, set<int>> G;
public:
    static int offsetToID(vector<int>& offset, int dim);
    static float GetScore(vector<float>& w, vector<float>& p, int dim);
    static bool RegionDominate(vector<point>& V, vector<float>& oi, vector<float>& oj, int dim); //true if oi<oj always in region
    static void GetCube(vector<point>& cube, vector<int>& offset, int pos, int dim);// compute vertices of the grid
    void GetGraph(vector<int>& offset, vector<vector<float>>& Allobj, int dim);
    static void EnumerateGrid(vector<int>&cur_offset, int cur_dim, int remains, int dim, vector<vector<float>>& Allobj, unordered_map<int, dominateG>& Grid);
    static int FindCube(point& v, int dim);
    static void MergeG(unordered_map<int, set<int>>& G, unordered_map<int, set<int>>& G_cube, set<int>& S);
};


#endif //K_LEVEL_DOMINATEG_H
