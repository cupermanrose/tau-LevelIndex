
#include "topk.h"
#include "oru.h"
#include <vector>


void dtopk::multiple_query(level& idx, int k, int q_num, fstream& log){
    vector<vector<float>> ws; // user preference weight vectors
    generate_query(idx, q_num, ws);
    clock_t cur_time=clock();
    int visit_sum=0;
    int result_sum=0;
    for (int i=0;i<q_num;i++){
        clock_t qb=clock();
        cout << "topk query " << i <<"("<<ws[i]<<")"<< ": " << endl;
        log << "topk query " << i <<"("<< ws[i]<<")"<< ": " << endl;
        single_query(idx,k, ws[i],visit_sum, result_sum, log);
        clock_t qe=clock();
        cout << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
        log << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
    }
    cout << "Average topk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}

void dtopk::generate_query(level& idx, int q_num, vector<vector<float>>& ws_ret){
    oru::generate_query(idx, q_num, ws_ret);
}

inline bool satisfy_constr(vector<float>& point, vector<float>& constr, bool sign){
    float dot=0;
    for (int i = 0; i <point.size(); ++i) {
        dot+=point[i]*constr[i];
    }
    return !((dot>constr[point.size()-1]) ^ sign);
}

bool point_in_convex(vector<float>& point, vector<halfspace> &convex){
    // check:
    // if sign:
    //     convex[:, 0:d-1] * point > convex[:, d]
    // else:
    //     convex[:, 0:d-1] * point < convex[:, d]
    assert(!convex.empty());
    assert(point.size()+1==convex[0].w.size());
    for (halfspace& constr: convex) {
        if(!satisfy_constr(point, constr.w, constr.side)){
            return false;
        }
    }
    return true;
}


void dtopk::single_query(level& idx, int k, vector<float> &w, int& visit_sum, int& result_sum, fstream& log){
    kcell *cur_cell=&idx.idx[0][0];
    int i=0;
    bool flag=true;
    visit_sum=0;
    while(flag){
        flag=false;
        if(cur_cell->curk>=idx.ik){
            vector<kcell> children;
            idx.SingleCellSplit(k,*cur_cell, children);
            for(kcell &c: children){
                visit_sum++;
                idx.UpdateH(c);
                if(point_in_convex(w, c.r.H)){
                    cur_cell=&c;
                    flag=true;
                    break;
                }
            }
        }else{
            for (int & it : cur_cell->Next){
                visit_sum++;
                idx.UpdateH(idx.idx[i+1][it]);
                if(point_in_convex(w, idx.idx[i+1][it].r.H)){
                    cur_cell=&idx.idx[i+1][it];
                    flag=true;
                    break;
                }
            }
        }
        ++i;
    }
    if(i<k){
        cout<<"visit cell number: "<<visit_sum<<endl;
        cout<<"query exists due to calculation accuracy when k="<<i<<endl;
        cout<<"the topi options are "<<cur_cell->topk<<endl;
    }

}
