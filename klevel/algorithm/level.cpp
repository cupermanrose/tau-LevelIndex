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
    //GlobalFilter(candidate);

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
            //LocalFilter(S1,Sk,*cur_cell); // rskyband or Gridfilter
            for (auto p=S1.begin();p!=S1.end();p++){
                if (!VerifyDuplicate(*p,*cur_cell,Sk, this_level)){ // merge s_tau
                    kcell newcell;
                    CreateNewCell(*p, S1, Sk, *cur_cell, newcell);
                    if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim)) // compute innerPoint
                        this_level.emplace_back(newcell);
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

bool level::VerifyDuplicate(int p, kcell &cur_cell, vector<int>& Sk, vector<kcell> &this_level) {
    bool flag = false;

    for (auto r=this_level.begin();r!=this_level.end();r++){
        if (r->objID!=p) continue;
        if (r->topk.find(p)==r->topk.end()) continue;
        bool found=true;
        for (auto it=r->topk.begin();it!=r->topk.end();it++){
            if(cur_cell.topk.find(*it)==cur_cell.topk.end()){
                found=false;
                break;
            }
        }
        if (found) {
            flag = true;
            for (auto it = Sk.begin(); it != Sk.end(); it++) {
                if (*it != p) r->Stau.insert(*it);
            }
            break;
        }
    }

    return flag;
}

void level::CreateNewCell(int p, vector<int> &S1, vector<int> &Sk, kcell &cur_cell, kcell &newcell) {
    newcell.curk=cur_cell.curk+1;
    newcell.objID=p;
    newcell.topk=cur_cell.topk; newcell.topk.insert(p);
    newcell.Stau.clear();
    for (auto it=Sk.begin();it!=Sk.end();it++){
        if (*it!=p) newcell.Stau.insert(*it);
    }

    // just for verification whether p can top-1 in cur_cell
    newcell.r.V.clear();
    newcell.r.H=cur_cell.r.H;

    for (auto it = S1.begin(); it != S1.end(); it++) {
        if (*it != p) AddHS(p,*it,true,newcell.r.H);
    }
    return;
}

void level::AddHS(int o1, int o2, bool side, vector<halfspace> &H) {
    halfspace hp;
    lp_adapter::ComputeHP(Allobj[o1],Allobj[o2],hp.w,dim);
    hp.side=side; // the halfspace score(p)>score(*it)
    H.emplace_back(hp);
}

void level::UpdateH(kcell &cur_cell) {
    // Replace halfspaces to compute the vertices of new region(merged)
    int p=cur_cell.objID;
    cur_cell.r.H.clear();
    for (auto it = cur_cell.topk.begin(); it != cur_cell.topk.end(); it++) {
        // the halfspace score(p)<score(*it)
        if (*it != p) AddHS(p,*it,false,cur_cell.r.H);
    }

    for (auto it = cur_cell.Stau.begin(); it != cur_cell.Stau.end(); it++) {
        // the halfspace score(p)>score(*it)
        AddHS(p,*it,true,cur_cell.r.H);
    }
}

void level::UpdateV(kcell &cur_cell) {
    qhull_adapter::ComputeVertex(cur_cell.r.H,cur_cell.r.innerPoint,dim);
}