//
// Created by cuper on 2021/4/16.
//

#include "kspr.h"

void kspr::generate_query(level &idx, int q_num, vector<int> &q_list) {
    srand(0); // random seed
    q_list.clear();
    for (int i=0;i<q_num;i++){
        q_list.push_back(rand()%idx.Allobj.size());
    }
    /*cout<<"begin generate query of original id:";
    for(int &i: q_list){
        cout<<idx.levelId_2_dataId[i]<<"\n";
    }
    cout<<"end generate query of original id"<<endl;

    cout<<"begin generate query of k-level id:";
    for(int &i: q_list){
        cout<<i<<"\n";
    }
    cout<<"end generate query of k-level id"<<endl;*/
    return;
}

bool kspr::Find_qid(kcell &this_cell, int qid) {
    for (auto it=this_cell.topk.begin();it!=this_cell.topk.end();it++)
        if (*it==qid) return true;
    for (auto it=this_cell.Stau.begin();it!=this_cell.Stau.end();it++)
        if (*it==qid) return true;
    return false;
}

int kspr::single_query(level &idx, int k, int q_id, fstream &log) {
    int cnt=0;
    for (auto it=idx.idx[k].begin();it!=idx.idx[k].end();it++){
        bool find=false;
        for (auto p=it->topk.begin();p!=it->topk.end();p++)
                if (*p==q_id) find=true;
        if (find) cnt++;
    }
    return cnt;
}

int kspr::single_query_largek(level &idx, int k, int q_id, fstream &log) {
    vector<int> S1,Sk;
    int ave_S1=0,ave_Sk=0,ave_vertex=0;
    vector<vector<kcell>> tmp; tmp.clear();
    vector<kcell> init_level;  init_level.clear();
    for (auto it=idx.idx[idx.ik].begin();it!=idx.idx[idx.ik].end();it++){
        if (Find_qid(*it,q_id)) init_level.push_back(*it);
    }
    tmp.emplace_back(init_level);

    for (int i=0;i<k-idx.ik;i++){
        vector<kcell> this_level;  this_level.clear(); idx.region_map.clear();
        for (auto cur_cell=tmp[i].begin(); cur_cell!=tmp[i].end(); cur_cell++){
            if (!Find_qid(*cur_cell,q_id)) continue; // will not contribute to kspr query
            idx.LocalFilter(S1,Sk,*cur_cell,ave_S1,ave_Sk);

            for (auto p=S1.begin();p!=S1.end();p++){

                kcell newcell;
                idx.CreateNewCell(*p, S1, Sk, *cur_cell, newcell);

                if (!idx.VerifyDuplicate(newcell,this_level)){ // just for profiling
                    if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,idx.dim)){ // just for profiling
                        this_level.emplace_back(newcell);
                        idx.region_map.insert(make_pair(newcell.hash_value,this_level.size()-1));
                    }
                }
            }
        }
        //Compute V for each cell
        //discuss why we need recompute after all
        for (auto cur_cell=this_level.begin();cur_cell!=this_level.end();cur_cell++){
            idx.UpdateH(*cur_cell);
            idx.UpdateV(*cur_cell,ave_vertex);
        }

        tmp.emplace_back(this_level);
    }
    int cnt=0;
    for (auto it=tmp.back().begin();it!=tmp.back().end();it++)
        if(Find_qid(*it,q_id)) cnt++;
    return cnt;
}

void kspr::multiple_query(level &idx, int k, int q_num, fstream &log) {
    clock_t cur_time=clock();
    vector<int> q_list;
    generate_query(idx,q_num, q_list);
    for (int i=0;i<q_num;i++){
        int answer;
        if (k<=idx.ik) answer=single_query(idx, k,q_list[i],log);
        else answer=single_query_largek(idx,k,q_list[i],log);
        cout << "The answer of kspr query " << i << ": " << answer << endl;
        log << "The answer of kspr query " << i << ": " << answer << endl;
    }
    cout << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}
