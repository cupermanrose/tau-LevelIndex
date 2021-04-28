//
// Created by cuper on 2021/4/16.
//

#include "utk.h"

void utk::generate_query(level &idx, int q_num, float utk_side_length, vector<vector<float>>& q_list) {
    srand(0); // random seed
    q_list.clear();
    for (int i=0;i<q_num;i++){
        vector<float> v((idx.dim-1)*2,0.0);
        float res=rand()/RAND_MAX-utk_side_length*(idx.dim-1);
        if (res<0) res=1.0-utk_side_length*(idx.dim-1);

        for (int d=0;d<idx.dim-1;d++){
            v[d*2] = res*rand()/RAND_MAX;
            v[d*2+1] = v[d*2]+utk_side_length;
            res=res-v[d*2];
        }
        q_list.push_back(v);
    }
    return;
}

bool utk::isIn(vector<float>& v, vector<float>& Qregion, int dim){
    float EPS_isIn=0.000001;
    for (int i = 0; i < dim - 1; i++) {
        if ((v[i] + EPS_isIn >= Qregion[2 * i]) && (v[i] - EPS_isIn <= Qregion[2 * i + 1])) continue;
        return false;
    }
    return true;
}

bool utk::isIn(vector<float> &v, vector<halfspace> &H, int dim) {
    for (auto it = H.begin(); it != H.end(); it++) {
        float sum = 0.0;
        for (int i = 0; i < dim - 1; i++) {
            sum = sum + v[i] * it->w[i];
        }
        if (it->side == false) {
            if (sum <= it->w[dim - 1]) continue;
            return false;
        }
        else {
            if (sum >= it->w[dim - 1]) continue;
            return false;
        }
    }
    return true;
}

bool utk::Intersect(vector<float> &Qregion, region& r, int dim) {
    if (r.V.empty()) return false;
    for (auto it=r.V.begin();it!=r.V.end();it++){
        if (isIn(*it,Qregion,dim)) return true;
    }

    int bitset = 1 << (dim-1);
    for (int i = 0; i < bitset; i++) {
        vector<float> Qv; Qv.clear();
        int tmp_i = i;
        for (int d = 0; d < dim - 1; d++) {
            if (tmp_i % 2 == 0) Qv.push_back(Qregion[d * 2]);
            else Qv.push_back(Qregion[d * 2 + 1]);
            tmp_i = tmp_i >> 1;
        }
        if (isIn(Qv, r.H, dim)) return true;
    }
    return false;
}

int utk::single_query(level &idx, int k, vector<float> &Qregion, fstream &log) {
    unordered_set<int> results; results.clear();
    for (auto it=idx.idx[k].begin();it!=idx.idx[k].end();it++){
        if (Intersect(Qregion,it->r, idx.dim)){
            for (auto p=it->topk.begin();p!=it->topk.end();p++){
                results.insert(*p);
            }
        }
    }
    return results.size();
}

int utk::single_query_largek(level &idx, int k, vector<float> &Qregion, fstream &log) {
    vector<int> S1,Sk;
    int ave_S1=0,ave_Sk=0,ave_vertex=0;
    vector<vector<kcell>> tmp; tmp.clear();
    vector<kcell> init_level;  init_level.clear();
    for (auto it=idx.idx[idx.ik].begin();it!=idx.idx[idx.ik].end();it++){
        if (Intersect(Qregion,it->r, idx.dim)) init_level.push_back(*it);
    }
    tmp.emplace_back(init_level);

    for (int i=0;i<k-idx.ik;i++){
        vector<kcell> this_level;  this_level.clear(); idx.region_map.clear();
        for (auto cur_cell=tmp[i].begin(); cur_cell!=tmp[i].end(); cur_cell++){
            if (!Intersect(Qregion,cur_cell->r, idx.dim)) continue; // will not contribute to utk query
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
    unordered_set<int> results; results.clear();
    for (auto it=tmp.back().begin();it!=tmp.back().end();it++) {
        if (Intersect(Qregion,it->r, idx.dim)){
            for (auto p=it->topk.begin();p!=it->topk.end();p++){
                results.insert(*p);
            }
        }
    }
    return results.size();
}

void utk::multiple_query(level &idx, int k, int q_num, float utk_side_length, fstream &log) {
    clock_t cur_time=clock();
    vector<vector<float>> q_list;

    generate_query(idx,q_num, utk_side_length, q_list);
    /*string queryfile="/home/jiahaozhang/data/klevel/query/utk/utk_query_4d_l0.001.txt";
    fstream queryout(queryfile, ios::out);
    for (int i=0;i<q_num;i++){
        for (int j=0;j<idx.dim-1;j++)
            queryout << q_list[i][2*j] << ' ' << q_list[i][2*j+1] << ' ';
        queryout << endl;
    }
    queryout.close();
    return;*/

    for (int i=0;i<q_num;i++){
        vector<int> utk_results;
        int answer;
        if (k<=idx.ik) answer=single_query(idx, k,q_list[i],log);
        else answer=single_query_largek(idx,k,q_list[i],log);
        cout << "The answer of utk query " << i << ": " << answer << endl;
        log << "The answer of utk query " << i << ": " << answer << endl;
    }
    cout << "Average utk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average utk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}

