
#include "level.h"
#include <chrono>
#include <algorithm>
#include <adapter/rtree_adapter.h>
#include <random>       // std::default_random_engine
#include <unistd.h>


bool apply_onion_from_file=false;
bool write_onion_to_file=false;
string anti_id_f;
string read_anti_dat;
level::level(int a_dim, int a_tau, int a_ik){
    dim=a_dim;
    tau=a_tau;
    ik=a_ik; // building levels
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

void level::LoadData(string datafile) {
    fstream fin(datafile, ios::in);
    if(!fin.is_open()){
        cout<<"fail to open "<<datafile<<endl;
        exit(-1);
    }
    OriginD.clear();
    vector<float> cl(dim);
    vector<float> cu(dim);
    while (true) {
        int id;

        fin >> id;
        if (fin.eof())
            break;

        vector<float> tmp; tmp.clear();
        for (int d = 0; d < dim; d++) fin >> cl[d];
        for (int d = 0; d < dim; d++) fin >> cu[d];
        for (int d = 0; d < dim; d++) tmp.push_back((cl[d]+cu[d])/2.0);

        OriginD.push_back(tmp);

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
    for (auto it=idx[k].begin();it!=idx[k].end();it++) it->FreeMem();
    idx[k].clear();
    vector<kcell>().swap(idx[k]);
}

void level::GlobalFilter(fstream& log, vector<int> &candidate) {
    candidate.clear();
    if(write_onion_to_file){
        build_onion(read_anti_dat, this->dim, this->tau);
        vector<vector<int>> onion;
        read_onion(anti_id_f, onion);
        int cnt=1;
        global_layer.clear();
        for(const auto &ly: onion){
            if(cnt>tau){
                break;
            }
            for(const auto &id:ly){
                candidate.push_back(id);
                global_layer.push_back(cnt);
            }
            ++cnt;
        }
//         for generate query file, do not delete
        vector<int> query(candidate.rbegin(), candidate.rbegin()+(10000<candidate.size()?10000:candidate.size()));
        for(auto &i: query){
            cout<<i+1<<endl;
        }
//        exit(0);
    }else if(apply_onion_from_file) {
        vector<vector<int>> onion;
        read_onion(anti_id_f, onion);
        int cnt=1;
        global_layer.clear();
        for(const auto &ly: onion){
            if(cnt>tau){
                break;
            }
            for(const auto &id:ly){
                candidate.push_back(id);
                global_layer.push_back(cnt);
            }
            ++cnt;
        }
        // for generate query file, do not delete
//        vector<int> query(candidate.rbegin(), candidate.rbegin()+(50<candidate.size()?50:candidate.size()));
//        for(auto &i: query){
//            cout<<i+1<<endl;
//        }
    }else {
        //kskyband
        vector<int> candidate_skyband, candidate_onionlayer;
        clock_t st = clock();
        kskyband(candidate_skyband, OriginD, tau);
        cout << "The number of options after kskyband: " << candidate_skyband.size() << std::endl;
        log << "The number of options after kskyband: " << candidate_skyband.size() << std::endl;
        cout << "kskyband Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC << endl;
        log << "kskyband Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC << endl;
        //k-onionlayer
        vector<int> layer;
        layer.clear();
        st = clock();
        onionlayer(candidate_onionlayer, layer, candidate_skyband, OriginD, tau);
        cout << "onionlayer Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC << endl;
        log << "onionlayer Cost: " << (clock() - st) / (float) CLOCKS_PER_SEC << endl;
        candidate = candidate_onionlayer;
        global_layer.clear();
        for(auto it=layer.begin();it!=layer.end();it++){
            global_layer.push_back(*it);
        }
    }
    Allobj.clear();
    int cnt=0;
    for (auto it=candidate.begin();it!=candidate.end();it++){
        Allobj.push_back(OriginD[*it]);
        levelId_2_dataId[cnt]=(*it)+1;
        ++cnt;
    }


    cout << "The number of options for building: " << Allobj.size() << std::endl;
    log << "The number of options for building: " << Allobj.size() << std::endl;
    std::cout << "GlobalFilter done!" << std::endl;
}

void level::initIdx(fstream& log){
    vector<int> candidate;
    GlobalFilter(log,candidate);
    idx.clear();
    kcell rootcell; rootcell.TobeRoot(candidate, dim); rootcell.Get_HashValue();
    vector<kcell> Lzero;
    Lzero={rootcell};
    idx.push_back(Lzero);

    cout << "Init done!" << endl;
}

void test_rt(Rtree& a_rtree, RtreeNode& node, set<int> &test, unordered_map<long int, RtreeNode *> &ramTree){
    if(node.isLeaf()){
        for (int i = 0; i <node.m_usedspace ; ++i) {
            assert(test.find(node.m_entry[i]->m_id+MAXPAGEID)==test.end());
            test.insert(node.m_entry[i]->m_id+MAXPAGEID);
        }
    }else{
        for (int i = 0; i <node.m_usedspace ; ++i) {
            assert(test.find(node.m_entry[i]->m_id)==test.end());
            test.insert(node.m_entry[i]->m_id);
            test_rt(a_rtree, *ramTree[node.m_entry[i]->m_id], test, ramTree);
        }
    }
}

void level::Build(fstream& log, ofstream& idxout) {
    vector<int> S1,Sk;
    set<int> utk_set; utk_set.clear();
    int ave_S1=0,ave_Sk=0,ave_vertex=0, cellsum=0,suc_split=0, HS_size=0;

    //for profiling
    clock_t tmp_profiling;
    double rskyband_time=0.0,verify_time=0.0,isFeasible_time=0.0,updateV_time=0.0;
    initIdx(log);
    clock_t level_zero_time=clock();
    for (int k=1;k<=ik;k++){
        clock_t level_k_time=clock();
        vector<kcell> this_level;
        region_map.clear();
        int valid_cell=0;
        for (auto cur_cell=idx[k-1].begin(); cur_cell!=idx[k-1].end(); cur_cell++){
            if ((cur_cell->curk!=0)&&(cur_cell->r.V.size()==0)) {
                continue;
            }
            valid_cell++;
            tmp_profiling=clock();
            LocalFilter(tau, S1,Sk,*cur_cell,ave_S1,ave_Sk);
            rskyband_time+=(clock()-tmp_profiling);

            for (auto p=S1.begin();p!=S1.end();p++){

                kcell newcell;
                CreateNewCell(*p, S1, Sk, *cur_cell, newcell);

                bool verify;
                tmp_profiling=clock();
                verify=VerifyDuplicate(newcell,this_level); // merge Stau
                verify_time+=(clock()-tmp_profiling);
                if (!verify){ // just for profiling
                    tmp_profiling=clock();
                    bool isFeasible=lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim); // compute innerPoint
                    isFeasible_time+=(clock()-tmp_profiling);
                    if (isFeasible){ // just for profiling
                        suc_split++;
                        utk_set.insert(newcell.objID);
                        this_level.emplace_back(newcell);
//                        UpdateH_S1(this_level.back(), S1);
                        region_map.insert(make_pair(newcell.hash_value,this_level.size()-1));
                    }else{
                        newcell.FreeMem();
                    }
                }else {
                    suc_split++;
                    newcell.FreeMem();
                }
            }
        }


        //Compute V for each cell
        for (auto cur_cell=this_level.begin();cur_cell!=this_level.end();cur_cell++){
            UpdateH(*cur_cell);
            HS_size+=cur_cell->r.H.size();
//            assert(lp_adapter::is_Feasible(cur_cell->r.H,cur_cell->r.innerPoint,dim)); // compute innerPoint
            tmp_profiling=clock();
            UpdateV(*cur_cell,ave_vertex);
            updateV_time+=(clock()-tmp_profiling);
        }
        cellsum+=this_level.size();
        idx.emplace_back(this_level);
        int ave_next = EdgeComputation(k-1);

        WriteToDisk(k-1, idxout);
        print_info(k, cellsum, valid_cell,level_zero_time,level_k_time,ave_S1,ave_Sk,ave_vertex, ave_next, suc_split, HS_size, utk_set, log);
        profiling(k,level_zero_time,rskyband_time,verify_time,isFeasible_time,updateV_time,log);
        FreeMem(k-1);
    }
    WriteToDisk(ik, idxout);

    cout << "The total size of index: " << cellsum << endl;
    log << "The total size of index: " << cellsum << endl;
}

void level::LocalFilter(int k, vector<int> &S1, vector<int> &Sk, kcell &cur_cell, int& ave_S1, int& ave_Sk) {
    if (cur_cell.curk==0){
        S1.clear();
        Sk.clear();
        for (int i=0;i<Allobj.size();i++){
            if (global_layer[i]==1) {
                S1.push_back(i);
            }
            Sk.push_back(i);
        }
    }
    else {
        rskyband(S1,Sk,cur_cell, k);
        //GridFilter(S1,Sk,cur_cell);
        //NoFilter(S1,Sk,cur_cell);
    }
    ave_Sk+=Sk.size();ave_S1+=S1.size();
}

void level::NoFilter(vector<int> &S1, vector<int> &Sk, kcell &cur_cell) {
    S1.clear();Sk.clear();
    for (auto it=cur_cell.Stau.begin();it!=cur_cell.Stau.end();it++){
        S1.push_back(*it);
        Sk.push_back(*it);
    }
}

void level::rskyband(vector<int> &S1, vector<int> &Sk, kcell &cur_cell, int k) {
    S1.clear();
    Sk.clear();
    if (cur_cell.r.V.size()==0) {
        return;
    }
    for (auto i=cur_cell.Stau.begin();i!=cur_cell.Stau.end();i++){
        int cnt=0;
        for (auto j=cur_cell.Stau.begin();j!=cur_cell.Stau.end();j++){
            if (*i==*j) {
                continue;
            }
            if (RegionDominate(cur_cell.r.V,Allobj[*i],Allobj[*j],dim)) {
                cnt++;
            }
            if (cnt>=(k-cur_cell.curk)) {
                break;
            }
        }
//        if ((cnt==0)&&(global_layer[*i]<=nextk)) {
        if (cnt==0) {
            S1.push_back(*i);
        }
        if (cnt<(k-cur_cell.curk)) {
            Sk.push_back(*i);
        }
    }
}

bool level::VerifyDuplicate(kcell &newcell, vector<kcell> &this_level) {
    bool flag = false;

    auto r_id=region_map.find(newcell.hash_value);
    if (r_id!=region_map.end()){
        flag=true;
        for (auto it = newcell.Stau.begin(); it != newcell.Stau.end(); it++) {
            bool foundInStau=false;
            for (auto p=this_level[r_id->second].Stau.begin();p!=this_level[r_id->second].Stau.end();p++){
                if (*p==*it) {
                    foundInStau=true;
                    break;
                }
            }
            if (!foundInStau) this_level[r_id->second].Stau.push_back(*it);
        }
    }

    return flag;
}

void level::CreateNewCell(int p, vector<int> &S1, vector<int> &Sk, kcell &cur_cell, kcell &newcell) {
    newcell.curk=cur_cell.curk+1;
    newcell.objID=p;
    newcell.topk.clear();
    newcell.topk=cur_cell.topk;
    newcell.topk.push_back(p);
    newcell.Get_HashValue();
    newcell.Stau.clear();
    for (auto it=Sk.begin();it!=Sk.end();it++){
        if (*it!=p) {
            newcell.Stau.push_back(*it);
        }
    }

    // just for verification whether p can top-1 in cur_cell
    newcell.r.V.clear();
    newcell.r.H.clear();
    newcell.r.H=cur_cell.r.H;

    for (auto it = S1.begin(); it != S1.end(); it++) {
        if (*it != p) {
//            assert(find(newcell.topk.begin(), newcell.topk.end(), *it)==newcell.topk.end());
            AddHS(p,*it,true,newcell.r.H);
        }
    }
    return;
}

void level::AddHS(int o1, int o2, bool side, vector<halfspace> &H) {
    halfspace hp;
    lp_adapter::ComputeHP(Allobj[o1],Allobj[o2],hp.w,dim);
    hp.side=side;
    bool flag=false;
    for (auto &i:hp.w) {
        if(i!=0){ // in case of 2 exact same options
            flag=true;
            break;
        }
    }
    if(flag){
        H.emplace_back(hp);
    }
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
//        assert(find(cur_cell.topk.begin(), cur_cell.topk.end(), *it)==cur_cell.topk.end());
        if (global_layer[*it]<=cur_cell.curk) {
            AddHS(p,*it,true,cur_cell.r.H);
        }
    }
}

void level::UpdateH_S1(kcell &cur_cell, vector<int> &S1) {
    // Replace halfspaces to compute the vertices of new region(merged)
    int p=cur_cell.objID;
    cur_cell.r.H.clear();
    for (auto topi:cur_cell.topk) {
        // the halfspace score(p)<score(*topi)
        if (topi != p) {
            AddHS(p, topi, false, cur_cell.r.H);
        }
    }

    for (auto otheri: S1) {
        // the halfspace score(p)>score(*otheri)
        if (global_layer[otheri] <= cur_cell.curk) {
            if(otheri!=p){
                AddHS(p, otheri, true, cur_cell.r.H);
            }
        }
    }
}

void level::UpdateV(kcell &cur_cell, int& ave_vertex) {
    if (dim>=3) qhull_adapter::ComputeVertex(cur_cell.r.H,cur_cell.r.V, cur_cell.r.innerPoint, dim);
    else qhull_adapter::ComputeVertex2D(cur_cell.r.H,cur_cell.r.V,cur_cell.r.innerPoint);
    ave_vertex+=cur_cell.r.V.size();
}

void level::print_info(int k, int cellsum, int valid_cell, clock_t& level_zero_time, clock_t & cur_time,
                       int& ave_S1, int& ave_Sk, int& ave_vertex, int& ave_next, int& suc_split, int& HS_size,
                       set<int>& utk_set, fstream& log) {
    cout << "LEVEL: " << k << endl;
    cout << "The region size of LEVEL " << k << ": " << idx[k].size() << endl;
    cout << "The option number of LEVEL 1-" << k << ": " << utk_set.size() << endl;
    cout << "Average S1 of LEVEL" << ": " << ave_S1 / (float)valid_cell << endl;
    cout << "Average Sk of LEVEL" << ": " << ave_Sk / (float)valid_cell << endl;
    cout << "Average splitting of LEVEL" << ": " << (float)suc_split / (float) valid_cell << endl;
    cout << "Average verterices of region in LEVEL" << ": " << (float) ave_vertex / (float) valid_cell << endl;
    cout << "Average next of region in LEVEL" << ": " << (float) ave_next / (float) valid_cell << endl;
    cout << "Average HS size of region in LEVEL" << ": " << (float) HS_size / (float) idx[k].size() << endl;
    cout << "Time Cost of LEVEL " << k << ": " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    cout << "Total region of LEVEL 1-" << k << ": " << cellsum << endl;
    cout << "Total Time Cost of LEVEL 1-" << k << ": " << (clock() - level_zero_time) / (float)CLOCKS_PER_SEC << endl;
    cout << endl;
    log << "LEVEL: " << k << endl;
    log << "The region size of LEVEL " << k << ": " << idx[k].size() << endl;
    log << "The option number of LEVEL 1-" << k << ": " << utk_set.size() << endl;
    log << "Average S1 of LEVEL" << ": " << ave_S1 / (float)valid_cell << endl;
    log << "Average Sk of LEVEL" << ": " << ave_Sk / (float)valid_cell << endl;
    log << "Average splitting of LEVEL" << ": " << (float)suc_split / (float)valid_cell << endl;
    log << "Average verterices of region in LEVEL" << ": " << (float) ave_vertex / (float) valid_cell << endl;
    log << "Average next of region in LEVEL" << ": " << (float) ave_next / (float) valid_cell << endl;
    log << "Average HS size of region in LEVEL" << ": " << (float) HS_size / (float) idx[k].size() << endl;
    log << "Time Cost of LEVEL " << k << ": " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    log << "Total region of LEVEL 1-" << k << ": " << cellsum << endl;
    log << "Total Time Cost of LEVEL 1-" << k << ": " << (clock() - level_zero_time) / (float)CLOCKS_PER_SEC << endl;
    log << endl;
    print_system_info(log);
    ave_Sk=ave_S1=ave_next=suc_split=HS_size=0; ave_vertex=0.0;
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
        // Space optimization
        //if (k<ik) it->WriteToDisk(idxout,false);
        //else it->WriteToDisk(idxout,true);
        it->WriteToDisk(idxout,true); // output Stau for all cells but not Halfspaces
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

int level::EdgeComputation(int k) {
    if (k==0){
        idx[k][0].Next.clear();
        for (int i=0;i<idx[k+1].size();i++) idx[k][0].Next.push_back(i);
        return idx[k+1].size();
    }
    unordered_map<size_t,vector<int>> k_map;
    for (int i=0;i<idx[k].size();i++){
        vector<int> tmp=idx[k][i].topk;
        std::sort(tmp.begin(),tmp.end()); // make it in order to output the unique hash_value for a set
        std::size_t seed = 0;
        for (auto it=tmp.begin();it!=tmp.end();it++){
            boost::hash_combine(seed, *it);
        }
        auto it=k_map.find(seed);
        if(it!=k_map.end()){
            it->second.emplace_back(i);
        }else{
            vector<int> id;
            id.emplace_back(i);
            k_map.emplace(seed, id);
        }
    }
    for (int i=0;i<idx[k+1].size();i++){
        int topk=idx[k+1][i].objID;
        vector<int> tmp;
        for(auto &id: idx[k+1][i].topk){
            if(id!=topk){
                tmp.emplace_back(id);
            }
        }
        std::sort(tmp.begin(),tmp.end()); // make it in order to output the unique hash_value for a set
        std:size_t seed=0;
        for(auto &id: tmp){
            boost::hash_combine(seed, id);
        }
        auto p_ids=k_map.find(seed);
        if (p_ids!=k_map.end()){
            for(auto &id:p_ids->second){
                idx[k][id].Next.push_back(i);
            }
        }
//        vector<int> tmp=idx[k+1][i].topk;
//        std::sort(tmp.begin(),tmp.end()); // make it in order to output the unique hash_value for a set
//        for (int j=0;j<k+1;j++){
//            std:size_t seed=0;
//            for (int t=0;t<k+1;t++){ // remove one option from idx[k+1][i] to find its parents
//                if (j==t) continue;
//                boost::hash_combine(seed, tmp[t]);
//            }
//            auto p_ids=k_map.find(seed);
//            if (p_ids!=k_map.end()){
//                for(auto &id:p_ids->second){
//                    idx[k][id].Next.push_back(i);
//                }
////                break;
//            }
//        }
    }
    int ave_next=0;
    for (auto it=idx[k].begin();it!=idx[k].end();it++) {
        ave_next+=it->Next.size();
//        if(it->Next.empty()){
//            ave_next+=it->Next.size();
//        }
    }
    return ave_next;
}

// Ins

void level::MergeCell(vector<kcell> &L_NoMerge, vector<kcell> &L_Merge) {
    L_Merge.clear();
    vector<kcell>().swap(L_Merge);
    region_map.clear();
    for (int i=0;i<L_NoMerge.size();i++){
        if (L_NoMerge[i].curk>ik) continue;
        L_NoMerge[i].Get_HashValue();
        if (!VerifyDuplicate(L_NoMerge[i],L_NoMerge)){
            region_map.insert(make_pair(L_NoMerge[i].hash_value,i));
        }
        else L_NoMerge[i].curk=ik+1;
    }
    for (int i=0;i<L_NoMerge.size();i++){
        if (L_NoMerge[i].curk>ik) continue;
        if (L_NoMerge[i].r.V.size()==0) continue;
        L_Merge.push_back(L_NoMerge[i]);
    }
}

void level::SplitCell(int p, int i, vector<kcell>& L) {
    // generate a new kcell within L[i].r
    bool flag = true;
    for (auto it=L[i].Stau.begin();it!=L[i].Stau.end();it++){
        if (RegionDominate(L[i].r.V, Allobj[p], Allobj[*it],dim)) {
            flag=false;
            break;
        }
    }
    if ((L[i].curk<ik)&&(global_layer[p]<=L[i].curk+1)&&(flag)){
        kcell newcell;
        newcell.curk=L[i].curk+1;
        newcell.objID=p;
        newcell.topk=L[i].topk;newcell.topk.push_back(p);
        newcell.Stau=L[i].Stau;
        newcell.r.H=L[i].r.H;
        for (auto it = newcell.topk.begin(); it != newcell.topk.end(); it++) {
            // the halfspace score(p)<score(*it)
            if (*it != p) AddHS(p,*it,false,newcell.r.H);
        }
        for (auto it = newcell.Stau.begin(); it != newcell.Stau.end(); it++) {
            // the halfspace score(p)>score(*it)
            AddHS(p,*it,true,newcell.r.H);
        }

        if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim)){
            int tmp;
            UpdateV(newcell,tmp);
            L.emplace_back(newcell);
        }
    }
    // update L[i].r
    if ((L[i].objID!=-1)&&(global_layer[p]<=L[i].curk)){
        if (RegionDominate(L[i].r.V, Allobj[p], Allobj[L[i].objID],dim)) {
            L[i].Stau.push_back(p);
        }
        else if (RegionDominate(L[i].r.V, Allobj[L[i].objID],Allobj[p], dim)){
            L[i].curk++;
            L[i].topk.push_back(p);
        }
        else if (flag) {
            // generate new kcell (R^-)
            kcell newcell;
            newcell.curk=L[i].curk+1;
            newcell.objID=L[i].objID;
            newcell.topk=L[i].topk;newcell.topk.push_back(p);
            newcell.Stau=L[i].Stau;
            newcell.r.H=L[i].r.H;
            AddHS(p,L[i].objID,true,newcell.r.H);
            if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim)){
                int tmp;
                UpdateV(newcell,tmp);
                L.emplace_back(newcell);
            }
            // update L[i] (R^+)
            L[i].Stau.push_back(p);
            AddHS(L[i].objID,p,true, L[i].r.H);
            if (!lp_adapter::is_Feasible(L[i].r.H,L[i].r.innerPoint,dim)){
                L[i].curk=tau+1; // delete this kcell
            }
            else{
                int tmp;
                UpdateV(L[i],tmp);
            }
        }
        else L[i].Stau.push_back(p);
    }
    else L[i].Stau.push_back(p);
}

void level::IncBuild(fstream& log, ofstream& idxout) {
    initIdx(log);
    kcell rootcell=idx[0][0]; rootcell.Stau.clear();
    vector<kcell> L,L_Merge; L={rootcell};

    int cnt=0;
    clock_t cur_time=clock();

    vector<int> ids(Allobj.size());
    iota(ids.begin(), ids.end(), 0);
//    shuffle(ids.begin(), ids.end(), default_random_engine(0));
    for(int id:ids){
//    for (int id=0;id<Allobj.size();id++){
        int size=L.size();
        for (int i=0;i<size;i++) {
            if ((L[i].r.V.size()==0)&(L[i].curk!=0)) L[i].curk=ik+1;
            if (L[i].curk>ik) {
                L[i].FreeMem();
                continue;
            }
            SplitCell(id, i, L);
        }
        cout << cnt++ << ": " << L.size() << endl;
    }

    MergeCell(L,L_Merge);
    L.clear();
    vector<kcell>().swap(L);
    for (int i=1;i<=ik;i++) idx.push_back(L);

    for (auto it=L_Merge.begin();it!=L_Merge.end();it++){
        if ((it->curk>ik)||(it->curk==0)) continue;
        idx[it->curk].push_back(*it);
    }
    //ComputeEdge();

    int cellsum=0;
    for (int i=0;i<=ik;i++) {
        cellsum+=idx[i].size();
        cout << "The region size of LEVEL " << i << ": " << idx[i].size() << endl;
        log << "The region size of LEVEL " << i << ": " << idx[i].size() << endl;
        WriteToDisk(i,idxout);
    }
    cout << "The total size of index: " << cellsum << endl;
    log << "The total size of index: " << cellsum << endl;

    cout << "Time Cost of IncBuild: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    log << "Time Cost of IncBuild: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;


    L_Merge.clear();
    vector<kcell>().swap(L_Merge);
}

void level::Build_nofilter(fstream& log, ofstream& idxout) {
    vector<int> S1,Sk;
    set<int> utk_set; utk_set.clear();
    int ave_S1=0,ave_Sk=0,ave_vertex=0, cellsum=0, suc_split=0, HS_size=0;

    //for profiling
    clock_t tmp_profiling;
    double rskyband_time=0.0,verify_time=0.0,isFeasible_time=0.0,updateV_time=0.0;

    initIdx(log);
    unordered_map<long int, RtreeNode *> ramTree;
    Rtree *rtree_rt= nullptr;
    build_rtree(rtree_rt, ramTree, Allobj);
    clock_t level_zero_time=clock();
    for (int k=1;k<=ik;k++){

        clock_t level_k_time=clock();
        vector<kcell> this_level;
        region_map.clear();
        int valid_cell=0;
        for (auto &cur_cell: idx[k-1]){
            if ((cur_cell.curk!=0)&&(cur_cell.r.V.size()==0)) continue;
            valid_cell++;
            tmp_profiling=clock();
            if (cur_cell.curk==0){
                S1.clear();
                Sk.clear();
                for (int i=0;i<Allobj.size();i++){
                    if (global_layer[i]==1) {
                        S1.push_back(i);
                    }
                    Sk.push_back(i);
                }
            }else{
                utk_rskyband(cur_cell.r.V, dim, *rtree_rt, S1, Allobj, ramTree, cur_cell.topk);
                utk_rskyband(cur_cell.r.V, dim, *rtree_rt, Sk, Allobj, ramTree, cur_cell.topk, k-cur_cell.curk);
//                unordered_set<int> s(cur_cell.Stau.begin(), cur_cell.Stau.end());
//                for (int i:cur_cell.topk) {
//                    s.erase(i);
//                }
//                Sk=vector<int>(s.begin(), s.end());
            }
            ave_S1+=S1.size();
            ave_Sk+=Sk.size();

//            NoFilter(S1,Sk,*cur_cell);
            rskyband_time+=(clock()-tmp_profiling);

            for (auto p=S1.begin();p!=S1.end();p++){
                kcell newcell;
                CreateNewCell(*p, S1, Sk, cur_cell, newcell);

                bool verify;
                tmp_profiling=clock();
                verify=VerifyDuplicate(newcell,this_level); // merge Stau
                verify_time+=(clock()-tmp_profiling);
                if (!verify){ // just for profiling
                    tmp_profiling=clock();
                    bool isFeasible=lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim); // compute innerPoint
                    isFeasible_time+=(clock()-tmp_profiling);
                    if (isFeasible){ // just for profiling
                        suc_split++;
                        utk_set.insert(newcell.objID);
                        this_level.emplace_back(newcell);
//                        UpdateH_S1(this_level.back(), S1);
//                        this_level.back().S1=S1;
                        region_map.insert(make_pair(newcell.hash_value,this_level.size()-1));
                    }else {
                        newcell.FreeMem();
                    }
                }else {
                    suc_split++;
                    newcell.FreeMem();
                }
            }

        }

        //Compute V for each cell
        for (auto cur_cell=this_level.begin();cur_cell!=this_level.end();cur_cell++){
            UpdateH(*cur_cell);
            HS_size+=cur_cell->r.H.size();
            tmp_profiling=clock();
            UpdateV(*cur_cell,ave_vertex);
            updateV_time+=(clock()-tmp_profiling);
        }
        cellsum+=this_level.size();
        idx.emplace_back(this_level);
        int ave_next = EdgeComputation(k-1);

        WriteToDisk(k-1, idxout);
        print_info(k, cellsum, valid_cell,level_zero_time,level_k_time,ave_S1,ave_Sk,ave_vertex, ave_next, suc_split, HS_size, utk_set,log);
        profiling(k,level_zero_time,rskyband_time,verify_time,isFeasible_time,updateV_time,log);
        FreeMem(k-1);
    }
    WriteToDisk(ik, idxout);
    cout << "The total size of index: " << cellsum << endl;
    log << "The total size of index: " << cellsum << endl;
}

// split next level cells from cur_cell, and all cells of NextCell are inside cur_cell
void level::SingleCellSplit(int qk, kcell &cur_cell, vector<kcell> &this_level) {
    if ((cur_cell.curk!=0)&&(cur_cell.r.V.size()==0)) return;
    vector<int> S1,Sk;
    int ave_S1=0,ave_Sk=0;
    LocalFilter(qk, S1,Sk,cur_cell,ave_S1,ave_Sk);
    for (auto p=S1.begin();p!=S1.end();p++){
        if (global_layer[*p]>cur_cell.curk+1) continue;
        kcell newcell;
        CreateNewCell(*p,S1,Sk,cur_cell,newcell);
        if (!VerifyDuplicate(newcell,this_level)){ // merge Stau
            if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim)) { // verify
                this_level.emplace_back(newcell);
                region_map.insert(make_pair(newcell.hash_value,this_level.size()-1));
            }
        }
    }
    return;
}

/*
// for large k queries
void level::SplitDFS(kcell& cell, vector<kcell> &L, ofstream& idxout, int& kcell_num) {
    if (cell.curk>=ik) return;
    vector<int> S1,Sk;
    int ave_S1=0,ave_Sk=0,ave_vertex=0;
    LocalFilter(tau, S1,Sk,cell,ave_S1,ave_Sk);
    for (auto p=S1.begin();p!=S1.end();p++){
        if (global_layer[*p]>cell.curk+1) continue;
        kcell newcell;
        //CreateNewCell(*p,S1,Sk,cell,newcell);
        newcell.curk=cell.curk+1;
        newcell.objID=*p;
        newcell.topk=cell.topk; newcell.topk.push_back(*p);
        newcell.Stau.clear();
        for (auto it=Sk.begin();it!=Sk.end();it++){
            if (*it!=*p) newcell.Stau.push_back(*it);
        }
        newcell.r.V.clear();
        newcell.r.H.clear();
        for (int i=0;i<newcell.topk.size();i++){
            for (int j=i+1;j<newcell.topk.size();j++){
                AddHS(newcell.topk[i],newcell.topk[j],true,newcell.r.H);
            }
        }
        for (auto it = S1.begin(); it != S1.end(); it++) {
            if (*it != *p) AddHS(*p,*it,true,newcell.r.H);
        }

        // verify
        if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,dim)) {
            kcell_num++;
            if (kcell_num%1000==0) cout << kcell_num <<endl;
            UpdateV(newcell, ave_vertex);
            SplitDFS(newcell,L, idxout, kcell_num);
            if (newcell.curk<ik) newcell.WriteToDisk(idxout,false);
            else newcell.WriteToDisk(idxout,true);
            //L.push_back(newcell);
        }
    }
    return;
}
*/



/*
void level::DFSBuild(fstream &log, ofstream &idxout) {
    vector<int> candidate; candidate.clear();
    GlobalFilter(log,candidate);
    kcell rootcell; rootcell.TobeRoot(candidate, dim);
    vector<kcell> L; L.clear();L.push_back(rootcell);


    clock_t cur_time=clock();
    int kcell_num=0;
    SplitDFS(rootcell,L, idxout,kcell_num);
    cout << "Index cell size: " << kcell_num << endl;
    cout << "Time Cost of DFSBuild: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    log << "Index cell size: " << kcell_num << endl;
    log << "Time Cost of DFSBuild: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    L.clear();
    vector<kcell>().swap(L);
}*/

/*
void level::SplitCell_Ins(int p, kcell& cur_cell, vector<kcell>& L) {
    bool flag = true;
    for (auto it=cur_cell.Stau.begin();it!=cur_cell.Stau.end();it++){ // check whether p can be top^(k+1)
        if (RegionDominate(cur_cell.r.V, Allobj[p], Allobj[*it],dim)) {
            flag=false;
            break;
        }
    }
    vector<int> S1,Stau;
    if (!flag){
        kcell newcell=cur_cell;
        newcell.Stau.push_back(p);
        L.emplace_back(newcell);
    }
    else {
        if ((cur_cell.curk < ik) && (global_layer[p] <= cur_cell.curk + 1)) { //c': generate a new kcell within L[i].r
            kcell newcell;
            S1=Stau=cur_cell.Stau;
            CreateNewCell(p, S1, Stau, cur_cell, newcell);
            if (lp_adapter::is_Feasible(newcell.r.H, newcell.r.innerPoint, dim)) L.emplace_back(newcell);
        }
        if (cur_cell.curk != 0) {
            // generate c^+ and c^-
            if (RegionDominate(cur_cell.r.V, Allobj[p], Allobj[cur_cell.objID], dim)) { // only c^+
                kcell newcell = cur_cell;
                newcell.Stau.push_back(p);
                L.emplace_back(newcell);
            }
            else if (RegionDominate(cur_cell.r.V, Allobj[cur_cell.objID], Allobj[p], dim) &&
                       (cur_cell.curk < ik)) { // only c^-
                kcell newcell = cur_cell;
                newcell.curk++;
                newcell.topk.push_back(p);
                L.emplace_back(newcell);
            }
            else {
                if (cur_cell.curk < ik) {
                    {// generate new kcell c^-
                        kcell newcell;
                        newcell.curk = cur_cell.curk + 1;
                        newcell.objID = cur_cell.objID;
                        newcell.topk = cur_cell.topk;
                        newcell.topk.push_back(p);
                        newcell.Stau = cur_cell.Stau;
                        newcell.r.H = cur_cell.r.H;
                        AddHS(p, cur_cell.objID, true, newcell.r.H);
                        if (lp_adapter::is_Feasible(newcell.r.H, newcell.r.innerPoint, dim)) L.emplace_back(newcell);
                    }
                    {// generate new kcell c^+
                        kcell newcell = cur_cell;
                        newcell.Stau.push_back(p);
                        AddHS(p, cur_cell.objID, false, newcell.r.H);
                        if (lp_adapter::is_Feasible(newcell.r.H, newcell.r.innerPoint, dim)) L.emplace_back(newcell);
                    }
                }
            }
        }
    }
}

void level::InsBuild(fstream& log, ofstream& idxout) {
    initIdx(log);
    kcell rootcell=idx[0][0]; rootcell.Stau.clear();
    vector<kcell> L, L_tmp; L={rootcell};

    int cnt=0, ave_vertex=0;
    clock_t cur_time=clock();
    for (int id=0;id<Allobj.size();id++){
        L_tmp.clear();
        vector<kcell>().swap(L_tmp);
        for (auto it=L.begin();it!=L.end();it++)  SplitCell_Ins(id,*it,L_tmp);
        if (id==Allobj.size()-1) MergeCell(L,L_tmp);
        L=L_tmp;
        for (auto it=L.begin();it!=L.end();it++) {
            //UpdateH(*it);
            if (it->r.H.empty()) it->r.V=rootcell.r.V;
            else UpdateV(*it,ave_vertex);
        }

        log << cnt << " : " << L.size() << endl;
        cout << cnt << " : " << L.size() << endl;
        cnt++;
        rootcell.Stau.push_back(id);
        L.push_back(rootcell);
    }


    for (auto it=L.begin();it!=L.end();it++){
        if ((it->curk>ik)||(it->curk==0)) continue;
        idx[it->curk].push_back(*it);
    }
    //ComputeEdge();
    int cellsum=0;
    for (int i=0;i<=ik;i++) {
        cellsum+=idx[i].size();
        cout << "The region size of LEVEL " << i << ": " << idx[i].size() << endl;
        log << "The region size of LEVEL " << i << ": " << idx[i].size() << endl;
        WriteToDisk(i,idxout);
    }
    cout << "The total size of index: " << cellsum << endl;
    log << "The total size of index: " << cellsum << endl;

    cout << "Time Cost of IncBuild: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;
    log << "Time Cost of IncBuild: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC << endl;

    L.clear();
    vector<kcell>().swap(L);
    L_tmp.clear();
    vector<kcell>().swap(L_tmp);
}
*/

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

