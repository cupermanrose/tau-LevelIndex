//
// Created by cuper on 2021/4/16.
//

#include "kspr.h"
#include <fstream>
#include <vector>

void kspr::generate_query(level &idx, int q_num, vector<int> &q_list) {
    srand(0); // random seed
    q_list.clear();
    /*fstream kspr_query("/home/jiahaozhang/data/klevel/query/kspr/kspr_id_10000.txt",ios::in);
    for (int i=0;i<q_num;i++){
        int id;
        for (int i=0;i<100;i++) kspr_query >> id;
        q_list.push_back(id);
        //q_list.push_back(rand()%idx.Allobj.size());
    }
    kspr_query.close();*/
    for (int i=0;i<q_num;i++){
        q_list.push_back(rand()%idx.Allobj.size());
    }

    /*cout<<"begin generate query of original id:\n";
    for(int &i: q_list){
        cout<<idx.levelId_2_dataId[i]<<"\n";
    }
    cout<<"end generate query of original id"<<endl;

    cout<<"begin generate query of k-level id:\n";
    for(int &i: q_list){
        cout<<i<<"\n";
    }
    cout<<"end generate query of k-level id"<<endl;*/
    return;
}

bool kspr::Find_qid_topk(kcell &this_cell, int qid) {
    for (auto it=this_cell.topk.begin();it!=this_cell.topk.end();it++)
        if (*it==qid) return true;
    return false;
}

bool kspr::Find_qid_Stau(kcell &this_cell, int qid) {
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

void LocalFilter_kspr(level& idx, kcell& cell, vector<int>&S1, vector<int>&Sk, int qk){
    S1.clear();Sk.clear();
    for (auto i=cell.Stau.begin();i!=cell.Stau.end();i++){
        int cnt=0;
        for (auto j=cell.Stau.begin();j!=cell.Stau.end();j++){
            if (*i==*j) continue;
            if (RegionDominate(cell.r.V,idx.Allobj[*i],idx.Allobj[*j],idx.dim)) cnt++;
            if (cnt>=(qk-cell.curk)) break;
        }
        if (cnt==0) S1.push_back(*i);
        if (cnt<(qk-cell.curk)) Sk.push_back(*i);
    }
}
void kspr::SplitDFS(level& idx, kcell& cell, int pq, int qk, int& cnt) {
    if (!Find_qid_Stau(cell,pq)) return;
    if ((cell.curk>=qk)||(cell.r.V.size()==0)) return;
    vector<int> S1,Sk;
    LocalFilter_kspr(idx,cell,S1,Sk,qk);
    for (auto p=S1.begin();p!=S1.end();p++){
        kcell newcell;
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
                idx.AddHS(newcell.topk[i],newcell.topk[j],true,newcell.r.H);
            }
        }
        for (auto it = S1.begin(); it != S1.end(); it++) {
            if (*it != *p) idx.AddHS(*p,*it,true,newcell.r.H);
        }

        // verify
        if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,idx.dim)) {
            if (*p==pq) cnt++;
            else{
                int ave_vertex=0;
                idx.UpdateV(newcell, ave_vertex);
                SplitDFS(idx,newcell,pq, qk,cnt);
            }
        }
    }
    return;
}

int kspr::single_query_largek(level &idx, int k, int q_id, fstream &log) {
    int cnt=0;
    vector<int> S1,Sk;
    for (auto it=idx.idx[idx.ik].begin();it!=idx.idx[idx.ik].end();it++){
        if (Find_qid_topk(*it,q_id)) cnt++;
        else if (Find_qid_Stau(*it,q_id)) {
            cnt++;
            SplitDFS(idx,*it,q_id,k,cnt);
        }
    }

    /*for (int i=0;i<k-idx.ik;i++){
        vector<kcell> this_level;  this_level.clear(); idx.region_map.clear();
        //cout << i+idx.ik << ' ' << tmp[i].size() << endl;
        for (auto cur_cell=tmp[i].begin(); cur_cell!=tmp[i].end(); cur_cell++){
            if (Find_qid_topk(*cur_cell,q_id)) {cnt++;continue;}
            if (!Find_qid_Stau(*cur_cell,q_id)) continue; // will not contribute to kspr query
            idx.LocalFilter(k, S1,Sk,*cur_cell,ave_S1,ave_Sk);

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
    for (auto it=tmp.back().begin();it!=tmp.back().end();it++)
        if(Find_qid_topk(*it,q_id)) cnt++;*/
    return cnt;
}

void kspr::multiple_query(level &idx, int k, int q_num, fstream &log) {
    clock_t invert_time = clock();

    vector<vector<int>> inverted_idx;
    vector<int> ids;ids.clear();
    int cur_k;
    if (k<=idx.ik) cur_k=k;
    else cur_k=idx.ik;
    inverted_idx.clear();
    for (int i=0;i<idx.Allobj.size();i++) inverted_idx.push_back(ids);
    for (int i=0;i<idx.idx[cur_k].size();i++){
        for (auto it=idx.idx[cur_k][i].topk.begin();it!=idx.idx[cur_k][i].topk.end();it++){
            inverted_idx[*it].push_back(i);
        }
    }

    log << "InvertedIdx from k-level building time: " << (clock() - invert_time) / (float) CLOCKS_PER_SEC << endl;
    cout << "InvertedIdx from k-level building time: " << (clock() - invert_time) / (float) CLOCKS_PER_SEC << endl;

    vector<int> q_list;
    generate_query(idx,q_num, q_list);
    clock_t cur_time=clock();
    for (int i=0;i<q_num;i++){
        int answer=0;
        for (auto it=idx.levelId_2_dataId.begin();it!=idx.levelId_2_dataId.end();it++){
            if (q_list[i]==it->second){
                if (k<=idx.ik) answer=inverted_idx[it->first].size();
                else answer=single_query_largek(idx,k,it->first,log);
                break;
            }
        }
        cout << "The answer of kspr query " << i <<"("<< q_list[i]<<")"<< ": " << answer << endl;
        log << "The answer of kspr query " << i <<"("<< q_list[i]<<")"<< ": " << answer << endl;
    }
    cout << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}
