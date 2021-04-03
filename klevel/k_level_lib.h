//
// Created by cuper on 2021/3/30.
//

#ifndef K_LEVEL_K_LEVEL_LIB_H
#define K_LEVEL_K_LEVEL_LIB_H

#include "algorithm/level.h"
#include <iostream>
#include <fstream>

class k_level_lib {
public:
    k_level_lib();
    ~k_level_lib();
    static void BuildIdx(int dim, int tau, char* datafile, char* logfile);
};


#endif //K_LEVEL_K_LEVEL_LIB_H
