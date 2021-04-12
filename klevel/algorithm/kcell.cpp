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

void kcell::Get_HashValue() {
    vector<int> tmp; tmp.clear();
    for (auto it=topk.begin();it!=topk.end();it++) tmp.push_back(*it);
    std::sort(tmp.begin(),tmp.end()); // make it in order to output the unique hash_value for a set
    std::size_t seed = 0;
    for (auto it=tmp.begin();it!=tmp.end();it++){
        boost::hash_combine(seed, *it);
    }

    boost::hash_combine(seed, objID); // objID will be computed twice, to make sure kcell shared the same hash value have the same objID
    hash_value=seed;
    return;
}