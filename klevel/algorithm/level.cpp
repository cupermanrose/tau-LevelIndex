//
// Created by cuper on 2021/3/30.
//

#include "level.h"
#include <fstream>

level::level(int a_dim, int a_tau){
    dim=a_dim;
    tau=a_tau;
    idx.clear();
    Allobj.clear();
}

level::~level() {
    idx.clear();
    vector<vector<kcell>>().swap(idx);
    Allobj.clear();
    vector<vector<float>>().swap(Allobj);
}

void level::LoadData(string datafile) {
    fstream fin(datafile, ios::in);
    Allobj.clear();
    //OriginData.clear();
    while (true) {
        int id;
        float *cl = new float[dim];
        float *cu = new float[dim];
        fin >> id;
        if (fin.eof())
            break;

        vector<float> tmp; tmp.clear();
        for (int d = 0; d < dim; d++) fin >> cl[d];
        for (int d = 0; d < dim; d++) fin >> cu[d];
        for (int d = 0; d < dim; d++) tmp.push_back((cl[d]+cu[d])/2.0);
        Allobj.push_back(tmp);

        /*if (TEST) {
            if (objCnt >= 50) break;
        }*/
        //log information
        if (Allobj.size() % 1000 == 0)
            cout << ".";
        if (Allobj.size() % 10000 == 0)
            cout << Allobj.size() << " objects loaded" << endl;
    }

    cout << "Total number of objects: " << Allobj.size() << endl;
    fin.close();
    return;
}

void level::initIdx(){
    vector<int> candidate;
    GlobalFilter(candidate);

    idx.clear();
    kcell rootcell; rootcell.TobeRoot(candidate, dim);
    vector<kcell> Lzero;
    Lzero={rootcell};
    idx.push_back(Lzero);
}

void level::Build() {
    vector<int> S1,Sk;
    initIdx();
    for (int k=1;k<=tau;k++){
        vector<kcell> this_level;  this_level.clear();
        for (auto cur_cell=idx[k-1].begin(); cur_cell!=idx[k-1].end(); cur_cell++){
            LocalFilter(S1,Sk,*cur_cell); // rskyband or Gridfilter
            for (auto p=S1.begin();p!=S1.end();p++){
                if (!VerifyDuplicate(*p,*cur_cell,this_level)){ // merge s_tau
                    CreateNewCell(*p, S1, Sk, *cur_cell, this_level);
                }
            }
        }

        //Compute V for each cell
        //discuss why we need recompute after all
        for (auto cur_cell=this_level.begin();cur_cell!=this_level.end();cur_cell++){
            UpdateH(*cur_cell);
            UpdateV(*cur_cell);
        }
        idx.emplace_back(this_level);
    }
}