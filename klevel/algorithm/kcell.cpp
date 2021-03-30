//
// Created by cuper on 2021/3/30.
//

#include "kcell.h"

kcell::kcell() {
    curk=-1;
    objID=-1;
    topk.clear();
    Sg.clear();
}

kcell::~kcell() {
}

void kcell::TobeRoot(vector<int> &candidates, int dim) {
    curk=0;
    objID=-1;
    topk.clear();
    for (auto it=candidates.begin();it!=candidates.end();it++) Sg.insert(*it);
    r.ToBeRoot(dim);
}