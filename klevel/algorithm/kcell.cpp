//
// Created by cuper on 2021/3/30.
//

#include "kcell.h"

kcell::kcell() {
    curk=-1;
    objID=-1;
    topk.clear();
    Stau.clear();
}

kcell::~kcell() {
}

void kcell::TobeRoot(vector<int> &candidates, int dim) {
    curk=0;
    objID=-1;
    topk.clear();
    for (int i=0;i!=candidates.size();i++) Stau.insert(i);
    r.ToBeRoot(dim);
}