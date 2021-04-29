//
// Created by cuper on 2021/4/28.
//

#include "oru.h"

void oru::generate_query(level &idx, int q_num, vector<vector<float>>& q_list) {
    srand(0); // random seed
    q_list.clear();
    for (int i=0;i<q_num;i++){
        vector<float> v(idx.dim-1,0.0);
        float res=1.0;
        for (int d=0;d<idx.dim-1;d++){
            v[d] = res*rand()/RAND_MAX;
            res=res-v[d];
        }
        q_list.push_back(v);
    }
    return;
}

bool oru::isIn(vector<float> &v, vector<halfspace> &H, int dim) {
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

float oru::GetDistance(vector<float>& q, region& r, int dim){
    if (isIn(q,r.H,dim)) return 0.0;
    float dis=1.0;
    for (auto it=r.V.begin();it!=r.V.end();it++){
        float tmp=0.0, resq=1.0,resv=1.0;
        for (int d=0;d<dim-1;d++){
            tmp=tmp+(q[d]-(*it)[d])*(q[d]-(*it)[d]);
            resq=resq-q[d];
            resv=resv-(*it)[d];
        }
        tmp=sqrt(tmp+resq*resv);
        if (tmp<dis) dis=tmp;
    }
    return dis;
}

float oru::kcell_filter(vector<kcell> &L, vector<bool>& filter, int ret_size, vector<float>& q, int dim) {
    vector<pair<float, int>> dis2q; dis2q.clear();
    for (int i=0;i<L.size();i++){
        float dis=GetDistance(q,L[i].r,dim);
        dis2q.push_back(make_pair(dis,i));
    }
    sort(dis2q.begin(),dis2q.end());
    unordered_set<int> ret;ret.clear();
    for (auto it=dis2q.begin();it!=dis2q.end();it++){
        int id=it->second;
        filter[id]=true;
        for (auto p=L[id].topk.begin();p!=L[id].topk.end();p++){
            ret.insert(*p);
        }
        if (ret.size()>=ret_size) return it->first;
    }
    return 0.0;
}

float oru::single_query(level &idx, int k, int ret_size, vector<float>& q, fstream &log) {
    vector<bool> filter(idx.idx[k].size(),false);
    return kcell_filter(idx.idx[k],filter,ret_size,q,idx.dim);
}

float oru::single_query_largek(level &idx, int k, int ret_size, vector<float>& q, fstream &log) {
    vector<int> S1,Sk;
    int ave_S1=0,ave_Sk=0,ave_vertex=0;
    // init level
    vector<vector<kcell>> tmp; tmp.clear();
    vector<kcell> init_level;  init_level.clear();
    vector<bool> filter(idx.idx[idx.ik].size(),false);
    float ans=kcell_filter(idx.idx[idx.ik],filter,ret_size,q,idx.dim);
    for (int i=0;i<idx.idx[idx.ik].size();i++){
        if (filter[i]) init_level.push_back(idx.idx[idx.ik][i]);
    }
    tmp.emplace_back(init_level);
    filter.clear();
    for (int i=0;i<tmp[0].size();i++) filter.push_back(true);


    for (int i=0;i<k-idx.ik;i++){
        vector<kcell> this_level;  this_level.clear(); idx.region_map.clear();
        cout << i+idx.ik << ' ' << tmp[i].size() << endl;
        int id=-1;
        for (auto cur_cell=tmp[i].begin(); cur_cell!=tmp[i].end(); cur_cell++){
            if (!filter[++id]) continue;

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

        filter.clear();
        for (int j=0;j<tmp[i+1].size();j++) filter.push_back(false);
        ans=kcell_filter(tmp[i+1],filter,ret_size,q,idx.dim);
    }
    return ans;;
}

void oru::multiple_query(level &idx, int k, int ret_size, int q_num, fstream &log) {
    clock_t cur_time=clock();
    vector<vector<float>> q_list;
    generate_query(idx,q_num, q_list);
    for (int i=0;i<q_num;i++){
        float answer;
        if (k<=idx.ik) answer=single_query(idx, k, ret_size, q_list[i],log);
        else answer=single_query_largek(idx,k,ret_size, q_list[i],log);
        cout << "The answer of oru query " << i << ": " << answer << endl;
        log << "The answer of oru query " << i << ": " << answer << endl;
    }
    cout << "Average oru query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average oru query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}
