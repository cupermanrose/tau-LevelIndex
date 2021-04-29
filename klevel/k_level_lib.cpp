//
// Created by cuper on 2021/3/30.
//

#include "k_level_lib.h"

void BuildIndex(level& idx, string datafile, fstream& log, string idxfile) {
    ofstream idxout(idxfile);
    idx.LoadData(datafile);
    idx.Build(log,idxout);
    log.close();
    idxout.close();
}

void LoadIndex(level& idx, string datafile, fstream& log, string idxfile) {
    ifstream idxin(idxfile);
    idx.LoadData(datafile);
    vector<int> candidate;
    idx.GlobalFilter(log,candidate); //TODO remove this and load from index file

    vector<kcell> empty_level; empty_level.clear();
    idx.idx.emplace_back(empty_level);
    for (int k=1;k<=idx.tau;k++){
        idx.idx.push_back(empty_level);
        idx.ReadFromDisk(k,idxin);
        cout << "The size of level " << k << ": "<<idx.idx[k].size() << endl;
        log << "The size of level " << k << ": "<<idx.idx[k].size() << endl;
    }
    log.close();
    idxin.close();
    cout << "Load Index Done!" << endl;
    log << "Load Index Done!" << endl;
    idx.print_system_info(log);
}
