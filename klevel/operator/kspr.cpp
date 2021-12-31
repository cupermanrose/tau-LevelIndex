
#include "kspr.h"
#include <fstream>
#include <vector>

void kspr::generate_query(level &idx, int q_num, vector<int> &q_list) {
    srand(100); // random seed
    for (int i=0;i<q_num;i++){
        q_list.push_back(rand()%idx.Allobj.size());
    }
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
//        for (auto it=idx.levelId_2_dataId.begin();it!=idx.levelId_2_dataId.end();it++){
//            if (q_list[i]==it->second){
                clock_t qb=clock();
                cout << "kspr query " << i <<"("<< q_list[i]<<")"<< ": " << endl;
                log << "kspr query " << i <<"("<< q_list[i]<<")"<< ": " << endl;
                single_query(idx,k,q_list[i],visit_sum, result_sum, log);
                clock_t qe=clock();
                cout << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
                log << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
//                break;
//            }
//        }
    }
    cout << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    cout << "Average visit cell: " << (float) visit_sum / (float) q_num << endl;
    log << "Average visit cell: " << (float) visit_sum / (float) q_num << endl;
    cout << "Average result cell: " << (float) result_sum / (float) q_num << endl;
    log << "Average result cell: " << (float) result_sum / (float) q_num << endl;
    return;
}

