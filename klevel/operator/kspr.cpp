
#include "kspr.h"
#include <fstream>
#include <vector>

void kspr::generate_query(level &idx, int q_num, vector<int> &q_list) {
    srand(0); // random seed
    q_list.clear();
    fstream kspr_query("~/data/klevel/query/kspr/kspr_id_10000.txt",ios::in);
    for (int i=0;i<q_num;i++){
        int id;
        for (int i=0;i<100;i++) kspr_query >> id;
        q_list.push_back(id);
        //q_list.push_back(rand()%idx.Allobj.size());
    }
    kspr_query.close();
    /*for (int i=0;i<q_num;i++){
        q_list.push_back(rand()%idx.Allobj.size());
    }
//    multimap<double, int> heap;
//    vector<float> w(idx.dim, 1.0/idx.dim);
//    for (int j = 0; j < idx.Allobj.size(); ++j) {
//        float s=GetScore(w, idx.Allobj[j], idx.dim);
//        heap.emplace(s, j);
//        if(heap.size()>q_num){
//            heap.erase(heap.begin());
//        }
//    }
//    for(auto &it:heap){
//        q_list.emplace_back(it.second);
//    }
    }*/

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

void kspr::single_query(level &idx, int k, int q_id, int& visit_sum, int& result_sum, fstream &log) {
    int visit=0, result=0, ave_vertex=0;
    vector<vector<kcell>> queue;
    for (int i=0;i<=k;i++){
        queue.push_back({});
    }
    queue[0].push_back(idx.idx[0][0]); // only contains rootcell
    set<int> results; results.clear();
    for (int i=0;i<k;i++){
        set<int> hash_set; hash_set.clear();
        idx.region_map.clear();
        for (int j=0;j<queue[i].size();j++){
            kcell& cur_cell=queue[i][j];
            if (Find_qid_topk(cur_cell,q_id)){
                result++;
                continue;
            }
            if (Find_qid_Stau(cur_cell,q_id)){
                if (cur_cell.curk<idx.ik){
                    for (auto it=cur_cell.Next.begin();it!=cur_cell.Next.end();it++){
                        if (hash_set.find(*it)==hash_set.end()) {
                            queue[i+1].push_back(idx.idx[i+1][*it]);
                            hash_set.insert(*it);
                        }
                    }
                }
                else{ // for large k
                    idx.SingleCellSplit(k,cur_cell,queue[i+1]);
                }
            }
        }
        for (int j=0;j<queue[i+1].size();j++){
            idx.UpdateH(queue[i+1][j]);
            idx.UpdateV(queue[i+1][j],ave_vertex);
        }
        visit+=queue[i+1].size();
    }
    visit_sum+=visit;
    result_sum+=result;
    cout << "Visiting cells of kspr query: " <<  visit << endl;
    log << "Visiting cells of kspr query: " <<  visit << endl;
    cout << "Result cells of kspr query: " << result << endl;
    log << "Result cells of kspr query: " << result << endl;
    return;
}

void kspr::multiple_query(level &idx, int k, int q_num, fstream &log) { // dag travel
    vector<int> q_list;
    generate_query(idx,q_num, q_list);
    clock_t cur_time=clock();
    int visit_sum=0;
    int result_sum=0;
    for (int i=0;i<q_num;i++){
        for (auto it=idx.levelId_2_dataId.begin();it!=idx.levelId_2_dataId.end();it++){
            if (q_list[i]==it->second){
                clock_t qb=clock();
                cout << "kspr query " << i <<"("<< q_list[i]<<")"<< ": " << endl;
                log << "kspr query " << i <<"("<< q_list[i]<<")"<< ": " << endl;
                single_query(idx,k,q_list[i],visit_sum, result_sum, log);
                clock_t qe=clock();
                cout << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
                log << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
                break;
            }
        }
    }
    cout << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    cout << "Average visit cell: " << (float) visit_sum / (float) q_num << endl;
    log << "Average visit cell: " << (float) visit_sum / (float) q_num << endl;
    cout << "Average result cell: " << (float) result_sum / (float) q_num << endl;
    log << "Average result cell: " << (float) result_sum / (float) q_num << endl;
    return;
}

/*
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
*/
