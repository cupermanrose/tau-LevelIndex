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
    return;
}

int kspr::single_query(level &idx, int k, int q_id, fstream &log) {
    int cnt=0;
    for (auto it=idx.idx[k].begin();it!=idx.idx[k].end();it++){
        if (it->topk.find(q_id)!=it->topk.end()) cnt++;
    }
    return cnt;
}

void kspr::multiple_query(level &idx, int k, int q_num, fstream &log) {
    clock_t cur_time=clock();
    vector<int> q_list;
    generate_query(idx,q_num, q_list);
    for (int i=0;i<q_num;i++){
        int answer=single_query(idx, k,q_list[i],log);
        cout << "The answer of kspr query " << i << ": " << answer << endl;
        log << "The answer of kspr query " << i << ": " << answer << endl;
    }
    cout << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}
