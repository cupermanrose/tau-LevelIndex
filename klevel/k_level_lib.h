//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_K_LEVEL_LIB_H
#define K_LEVEL_K_LEVEL_LIB_H

#include "algorithm/level.h"
#include <iostream>
#include <fstream>

void BuildIndex(level& idx, string datafile, fstream& log, string idxfile);
void LoadIndex(level& idx, string datafile, fstream& log, string idxfile);
void IncBuildIndex(level& idx, string datafile, fstream& log, string idxfile);
void DFSBuildIndex(level& idx, string datafile, fstream& log, string idxfile);

#endif //K_LEVEL_K_LEVEL_LIB_H
