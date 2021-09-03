

#ifndef K_LEVEL_K_LEVEL_LIB_H
#define K_LEVEL_K_LEVEL_LIB_H

#include "algorithm/level.h"
#include "rtree_adapter.h"
#include "hypercube.h"
#include <iostream>
#include <fstream>

void PTBuild(level& idx, string datafile, fstream& log, string idxfile);
void PTFBuild(level& idx, string datafile, fstream& log, string idxfile);
void INSBuild(level& idx, string datafile, fstream& log, string idxfile);

void LoadIndex(level& idx, string datafile, fstream& log, string idxfile);
void BuildRtree(const vector<kcell> &L, Rtree* &rt, unordered_map<long int, RtreeNode*>& ramTree, int dim);
void Vertex2BOX(const vector<kcell> &L, vector<vector<float>>& MBRs, int dim);
void RangeQueryFromRtree(Rtree* rt, unordered_map<long int, RtreeNode*> &ramTree,
                         vector<float> &ql, vector<float> &qu, vector<int> &ret_ids);
#endif //K_LEVEL_K_LEVEL_LIB_H
