//
// Created by cuper on 2021/3/30.
//

#include "k_level_lib.h"

void k_level_lib::BuildIdx(int dim, int tau, char* datafile, char* logfile) {
    fstream log(logfile, ios::out);
    level idx(dim,tau);
    idx.LoadData(datafile);
    idx.Build(log);
    log.close();
}
