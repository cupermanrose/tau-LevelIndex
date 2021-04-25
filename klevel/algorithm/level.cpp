//
// Created by cuper on 2021/3/30.
//

#include "level.h"
#include <chrono>
level::level(int a_dim, int a_tau){
    dim=a_dim;
    tau=a_tau;
    idx.clear();
    Allobj.clear();
    //Grid.clear();
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

        /*if (TEST) {
            if (OriginD.size() >= 50) break;
        }*/
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

void level::FreeMem(int k){
    if (k<0) return;
    for (auto it=idx[k].begin();it!=idx[k].end();it++){
        it->Stau.clear();
        unordered_set<int>().swap(it->Stau);
        it->topk.clear();
        unordered_set<int>().swap(it->topk);
        it->r.FreeMem();
    }
    idx[k].clear();
    vector<kcell>().swap(idx[k]);
}

void level::GlobalFilter(fstream& log, vector<int> &candidate) {
    candidate.clear();

    //kskyband
    vector<int> candidate_skyband, candidate_onionlayer;
    clock_t st=clock();
    kskyband(candidate_skyband, OriginD,tau);
    cout << "The number of options after kskyband: " << candidate_skyband.size() << std::endl;
    log << "The number of options after kskyband: " << candidate_skyband.size() << std::endl;
    cout << "kskyband Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC << endl;
    log << "kskyband Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC  << endl;
    //k-onionlayer
    vector<int> layer; layer.clear();
    st=clock();
    onionlayer(candidate_onionlayer, layer,  candidate_skyband,  OriginD, tau);
    cout << "onionlayer Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC << endl;
    log << "onionlayer Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC  << endl;

    candidate=candidate_onionlayer;

    Allobj.clear();
    int cnt=0;
    for (auto it=candidate.begin();it!=candidate.end();it++){
        Allobj.push_back(OriginD[*it]);
        levelId_2_dataId[cnt]=*it;
        ++cnt;
    }
    global_layer.clear();
    for(auto it=layer.begin();it!=layer.end();it++){
        global_layer.push_back(*it);
    }

    cout << "The number of options for building: " << Allobj.size() << std::endl;
    log << "The number of options for building: " << Allobj.size() << std::endl;
    std::cout << "GlobalFilter done!" << std::endl;
}

void level::initIdx(fstream& log){
    vector<int> candidate;
    GlobalFilter(log,candidate);

    // Initial Grid
    /*vector<int> offset(dim,0);
    dominateG::EnumerateGrid(offset,0,div_num,dim,Allobj,Grid);
    std::cout << "Init Grid done!" << std::endl;*/

    idx.clear();
    kcell rootcell; rootcell.TobeRoot(candidate, dim); rootcell.Get_HashValue();
    vector<kcell> Lzero;
    Lzero={rootcell};
    idx.push_back(Lzero);

    cout << "Init done!" << endl;
}

void level::Build(fstream& log, ofstream& idxout) {
    vector<int> S1,Sk;
    set<int> utk_set; utk_set.clear();
    int ave_S1=0,ave_Sk=0,ave_vertex=0;

    //for profiling
    clock_t tmp_profiling;
    double rskyband_time=0.0,verify_time=0.0,isFeasible_time=0.0,updateV_time=0.0;

    initIdx(log);
    clock_t level_zero_time=clock();
    for (int k=1;k<=tau;k++){

        clock_t level_k_time=clock();

        vector<kcell> this_level;  this_level.clear(); region_map.clear();
        for (auto cur_cell=idx[k-1].begin(); cur_cell!=idx[k-1].end(); cur_cell++){

            tmp_profiling=clock();
            LocalFilter(S1,Sk,*cur_cell,ave_S1,ave_Sk);
            rskyband_time+=(clock()-tmp_profiling);

            for (auto p=S1.begin();p!=S1.end();p++){

                kcell newcell;
                CreateNewCell(*p, S1, Sk, *cur_cell, newcell);

                bool verify;
                tmp_profiling=clock();
                //verify=VerifyDuplicate(*p,*cur_cell,Sk, this_level);
                verify=VerifyDuplicate(newcell,this_level);
                verify_time+=(clock()-tmp_profiling);
                if (!verify){ // just for profiling
                //if (!VerifyDuplicate(*p,*cur_cell,Sk, this_level)){ // merge s_tau
                    bool isFeasible;
                    tmp_profiling=clock();
                    isFeasible=lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim);
                    isFeasible_time+=(clock()-tmp_profiling);
                    if (isFeasible){ // just for profiling
                    //if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim)) { // compute innerPoint
                        utk_set.insert(newcell.objID);
                        this_level.emplace_back(newcell);
                        region_map.insert(make_pair(newcell.hash_value,this_level.size()-1));
                    }
                }
            }

        }

        //Compute V for each cell
        //discuss why we need recompute after all
        for (auto cur_cell=this_level.begin();cur_cell!=this_level.end();cur_cell++){
            UpdateH(*cur_cell);
            tmp_profiling=clock();
            UpdateV(*cur_cell,ave_vertex);
            updateV_time+=(clock()-tmp_profiling);
        }

        idx.emplace_back(this_level);

        WriteToDisk(k, idxout);
        print_info(k,level_zero_time,level_k_time,ave_S1,ave_Sk,ave_vertex, utk_set,log);
        profiling(k,level_zero_time,rskyband_time,verify_time,isFeasible_time,updateV_time,log);
        FreeMem(k-1);
    }
}

void level::LocalFilter(vector<int> &S1, vector<int> &Sk, kcell &cur_cell, int& ave_S1, int& ave_Sk) {
    rskyband(S1,Sk,cur_cell);
    //GridFilter(S1,Sk,cur_cell);
    //NoFilter(S1,Sk,cur_cell);
    ave_Sk+=Sk.size();ave_S1+=S1.size();
}

void level::NoFilter(vector<int> &S1, vector<int> &Sk, kcell &cur_cell) {
    S1.clear();Sk.clear();
    for (auto it=cur_cell.Stau.begin();it!=cur_cell.Stau.end();it++){
        S1.push_back(*it);
        Sk.push_back(*it);
    }
}

void level::rskyband(vector<int> &S1, vector<int> &Sk, kcell &cur_cell) {
    int nextk=cur_cell.curk+1;
    S1.clear();Sk.clear();
    if (cur_cell.r.V.size()==0) return;
    for (auto i=cur_cell.Stau.begin();i!=cur_cell.Stau.end();i++){
        int cnt=0;
        for (auto j=cur_cell.Stau.begin();j!=cur_cell.Stau.end();j++){
            if (*i==*j) continue;
            if (RegionDominate(cur_cell.r.V,Allobj[*i],Allobj[*j],dim)) cnt++;
            if (cnt>=(tau-cur_cell.curk)) break;
        }
        if ((cnt==0)&&(global_layer[*i]<=nextk)) S1.push_back(*i);
        if (cnt<(tau-cur_cell.curk)) Sk.push_back(*i);
    }
}

/*
void level::GridFilter(vector<int> &S1, vector<int> &Sk, kcell &cur_cell) {
    unordered_map<int, set<int>> G; G.clear();
    set<int> id; id.clear();
    for (auto v=cur_cell.r.V.begin();v!=cur_cell.r.V.end();v++){
        int cube_id=dominateG::FindCube(*v,dim);
        id.insert(cube_id);
    }
    for (auto cube_id=id.begin();cube_id!=id.end();cube_id++){
        auto cube=Grid.find(*cube_id);
        if (cube!=Grid.end()) dominateG::MergeG(G,cube->second.G,cur_cell.Stau);
    }

    S1.clear();Sk.clear();
    for (auto p=G.begin();p!=G.end();p++){
        if (p->second.size()==0) S1.push_back(p->first);
        if (p->second.size()<(tau-cur_cell.curk)) Sk.push_back(p->first);
    }
}
*/

bool level::VerifyDuplicate(kcell &newcell, vector<kcell> &this_level) {
    bool flag = false;

    auto r_id=region_map.find(newcell.hash_value);
    if (r_id!=region_map.end()){
        flag=true;
        for (auto it = newcell.Stau.begin(); it != newcell.Stau.end(); it++) {
            this_level[r_id->second].Stau.insert(*it);
        }
    }

    return flag;
}

/*
bool level::VerifyDuplicate(int p, kcell &cur_cell, vector<int>& Sk, vector<kcell> &this_level) {
    bool flag = false;

    for (auto r=this_level.begin();r!=this_level.end();r++){
        if (r->objID!=p) continue;
        if (r->topk.find(p)==r->topk.end()) continue;

        bool found=true;
        for (auto it=cur_cell.topk.begin();it!=cur_cell.topk.end();it++){
            if (r->topk.find(*it)==r->topk.end()){
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
*/

void level::CreateNewCell(int p, vector<int> &S1, vector<int> &Sk, kcell &cur_cell, kcell &newcell) {
    newcell.curk=cur_cell.curk+1;
    newcell.objID=p;
    newcell.topk=cur_cell.topk; newcell.topk.insert(p);
    newcell.Get_HashValue();
    newcell.Stau.clear();
    for (auto it=Sk.begin();it!=Sk.end();it++){
        if (*it!=p) newcell.Stau.insert(*it);
    }

    // just for verification whether p can top-1 in cur_cell
    newcell.r.V.clear();
    newcell.r.H.clear();
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

void level::UpdateV(kcell &cur_cell, int& ave_vertex) {
    qhull_adapter::ComputeVertex(cur_cell.r.H,cur_cell.r.V, cur_cell.r.innerPoint, dim);
    ave_vertex+=cur_cell.r.V.size();
}

void level::print_info(int k, clock_t& level_zero_time, clock_t & cur_time, int& ave_S1, int& ave_Sk, int& ave_vertex, set<int>& utk_set, fstream& log) {
    cout << "LEVEL: " << k << endl;
    cout << "The region size of LEVEL " << k << ": " << idx[k].size() << endl;
    cout << "The option number of LEVEL 1-" << k << ": " << utk_set.size() << endl;
    cout << "Average S1 of LEVEL" << ": " << ave_S1 / (float)idx[k-1].size() << endl;
    cout << "Average Sk of LEVEL" << ": " << ave_Sk / (float)idx[k-1].size() << endl;
    cout << "Average splitting of LEVEL" << ": " << (float)idx[k].size() / (float)idx[k-1].size() << endl;
    cout << "Average verterices of region in LEVEL" << ": " << (float) ave_vertex / (float) idx[k-1].size() << endl;
    cout << "Time Cost of LEVEL " << k << ": " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    cout << "Total Time Cost of LEVEL 1-" << k << ": " << (clock() - level_zero_time) / (float)CLOCKS_PER_SEC << endl;
    cout << endl;
    log << "LEVEL: " << k << endl;
    log << "The region size of LEVEL " << k << ": " << idx[k].size() << endl;
    log << "The option number of LEVEL 1-" << k << ": " << utk_set.size() << endl;
    log << "Average S1 of LEVEL" << ": " << ave_S1 / (float)idx[k-1].size() << endl;
    log << "Average Sk of LEVEL" << ": " << ave_Sk / (float)idx[k-1].size() << endl;
    log << "Average splitting of LEVEL" << ": " << (float)idx[k].size() / (float)idx[k-1].size() << endl;
    log << "Average verterices of region in LEVEL" << ": " << (float) ave_vertex / (float) idx[k-1].size() << endl;
    log << "Time Cost of LEVEL " << k << ": " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    log << "Total Time Cost of LEVEL 1-" << k << ": " << (clock() - level_zero_time) / (float)CLOCKS_PER_SEC << endl;
    log << endl;
    print_system_info(log);
    ave_Sk=ave_S1=0; ave_vertex=0.0;
}

void level::print_system_info(fstream &log) {
    using std::ios_base;
    using std::ifstream;
    using std::string;

    double vm_usage     = 0.0;
    double resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0 /1024.0; //MB
    resident_set = rss * page_size_kb /1024.0; //MB

    //cout << "Virtual Set Size(VSS) Cost: " << vm_usage << " MB" << endl;
    cout << "Resident Set Size(RSS) Cost: " << resident_set << " MB"  << endl;
    //log << "Virtual Set Size(VSS) Cost: " << vm_usage << " MB" << endl;
    log << "Resident Set Size(RSS) Cost: " << resident_set << " MB"  << endl;
}

void level::profiling(int k, clock_t &level_zero_time, double &rskyband_time, double &verify_time,
                      double &isFeasible_time, double &updateV_time, fstream &log) {
    cout << "rskyband time (%) of Level " << k << ": " << rskyband_time/(clock()-level_zero_time) * 100.0 <<endl;
    cout << "Verify Duplicate time (%) of Level " << k << ": " << verify_time/(clock()-level_zero_time) * 100.0 << endl;
    cout << "Is_Feasible time (%) of Level " << k << ": " << isFeasible_time/(clock()-level_zero_time) * 100.0 << endl;
    cout << "UpdateV time (%) of Level " << k << ": " << updateV_time/(clock()-level_zero_time) * 100.0 << endl;
    cout << endl;
    log << "rskyband time (%) of Level " << k << ": " << rskyband_time/(clock()-level_zero_time) * 100.0 << endl;
    log << "Verify Duplicate time (%) of Level " << k << ": " << verify_time/(clock()-level_zero_time) * 100.0 <<endl;
    log << "Is_Feasible time (%) of Level " << k << ": " << isFeasible_time/(clock()-level_zero_time) * 100.0 << endl;
    log << "UpdateV time (%) of Level " << k << ": " << updateV_time/(clock()-level_zero_time) * 100.0 <<  endl;
    log << endl;

}

void level::WriteToDisk(int k, ofstream &idxout) {
    int size=idx[k].size();
    idxout.write((char*) &size, sizeof(int));
    for (auto it=idx[k].begin();it!=idx[k].end();it++){
        it->WriteToDisk(idxout);
    }
}

void level::ReadFromDisk(int k, ifstream &idxin) {
    int size;
    idxin.read((char*) &size, sizeof(int));
    idx[k].clear();
    for (int i=0;i<size;i++){
        kcell tmp;
        tmp.ReadFromDisk(idxin);
        idx[k].emplace_back(tmp);
    }
}