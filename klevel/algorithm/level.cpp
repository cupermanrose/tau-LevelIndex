//
// Created by cuper on 2021/3/30.
//

#include "level.h"

level::level(int a_dim, int a_tau){
    dim=a_dim;
    tau=a_tau;
    idx.clear();
    Allobj.clear();
    Grid.clear();
}

level::~level() {
    idx.clear();
    vector<vector<kcell>>().swap(idx);
    Allobj.clear();
    vector<vector<float>>().swap(Allobj);
}

void level::LoadData(char* datafile) {
    fstream fin(datafile, ios::in);
    OriginD.clear();
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
        OriginD.push_back(tmp);

        if (TEST) {
            if (OriginD.size() >= 50) break;
        }
        //log information
        if (OriginD.size() % 1000 == 0)
            cout << ".";
        if (OriginD.size() % 10000 == 0)
            cout << OriginD.size() << " objects loaded" << endl;
    }

    cout << "Total number of objects: " << OriginD.size() << endl;
    fin.close();
    return;
}

void level::GlobalFilter(vector<int> &candidate) {
    candidate.clear();
    for (int i=0;i<OriginD.size();i++) candidate.push_back(i);
    //kskyband
    //k-onionlayer
    std::cout << "GlobalFilter done!" << std::endl;
}

void level::initIdx(fstream& log){
    vector<int> candidate;
    GlobalFilter(candidate);

    Allobj.clear();
    for (auto it=candidate.begin();it!=candidate.end();it++){
        Allobj.push_back(OriginD[*it]);
    }

    // Initial Grid
    vector<int> offset(dim,0);
    dominateG::EnumerateGrid(offset,0,div_num,dim,Allobj,Grid);
    std::cout << "Init Grid done!" << std::endl;

    idx.clear();
    kcell rootcell; rootcell.TobeRoot(candidate, dim);
    vector<kcell> Lzero;
    Lzero={rootcell};
    idx.push_back(Lzero);

    cout << "Init done!" << endl;
}

void level::Build(fstream& log) {
    vector<int> S1,Sk;
    initIdx(log);
    for (int k=1;k<=tau;k++){
        vector<kcell> this_level;  this_level.clear();
        for (auto cur_cell=idx[k-1].begin(); cur_cell!=idx[k-1].end(); cur_cell++){
            //LocalFilter: rskyband or Gridfilter
            rskyband(S1,Sk,*cur_cell);
            //GridFilter(S1,Sk,*cur_cell);
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

        cout << "LEVEL: " << k << endl;
        cout << "The region size of LEVEL " << k << ": " << idx[k].size() << endl;
        //cout << "The utk size of LEVEL " << k << ": " << utk_results.size() << endl;
        //cout << "Time Cost of LEVEL " << k << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;
        ///cout << "Average one_skyband size of LEVEL" << ": " << one_candidate / (float)L[level - 1].size() << endl;
        //cout << "Average spliting number of LEVEL" << ": " << spliting / (float)L[level - 1].size() << endl;
        //cout << "Average candidates of LEVEL" << ": " << ave_candidate / (float)L[level - 1].size() << endl;
        //cout << "Average verterices of region in LEVEL" << ": " << ave_hull_vertex / (float)L[level - 1].size() << endl;
        //cout << "subgraph_pruning in LEVEL" << ": " << subgraph_pruning / (float)L[level - 1].size() << endl;

    }
}

void level::rskyband(vector<int> &S1, vector<int> &Sk, kcell &cur_cell) {
    S1.clear();Sk.clear();
    for (auto i=cur_cell.Stau.begin();i!=cur_cell.Stau.end();i++){
        int cnt=0;
        for (auto j=cur_cell.Stau.begin();j!=cur_cell.Stau.end();j++){
            if (*i==*j) continue;
            if (dominateG::RegionDominate(cur_cell.r.V,Allobj[*i],Allobj[*j],dim)) cnt++;
        }
        if (cnt==0) S1.push_back(*i);
        if (cnt<(tau-cur_cell.curk)) Sk.push_back(*i);
    }
}

void level::GridFilter(vector<int> &S1, vector<int> &Sk, kcell &cur_cell) {
    unordered_map<int, set<int>> G; G.clear();
    for (auto v=cur_cell.r.V.begin();v!=cur_cell.r.V.end();v++){
        int cube_id=dominateG::FindCube(*v,dim);
        dominateG::MergeG(G,Grid[cube_id].G,cur_cell.Stau);
    }
    S1.clear();Sk.clear();
    for (auto p=G.begin();p!=G.end();p++){
        if (p->second.size()==0) S1.push_back(p->first);
        if (p->second.size()<(tau-cur_cell.curk)) Sk.push_back(p->first);
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
    qhull_adapter::ComputeVertex(cur_cell.r.H,cur_cell.r.V, cur_cell.r.innerPoint, dim);
}