//
// Created by cuper on 2021/4/28.
//

#include "oru.h"
#include "qp_adapter.h"
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
    double eps=0.0001;
    for (auto it = H.begin(); it != H.end(); it++) {
        float sum = 0.0;
        for (int i = 0; i < dim - 1; i++) {
            sum = sum + v[i] * it->w[i];
        }
        if (it->side == false) {
            if (sum+eps >= it->w[dim - 1]) return false;
        }
        else {
            if (sum-eps <= it->w[dim - 1]) return false;
        }
    }
    return true;
}

float oru::GetDistance(vector<float>& q, region& r, int dim){
    vector<vector<float>> halfspaces;
    for(auto &i: r.H){
        halfspaces.emplace_back(i.w);
        if(i.side){
            for(auto &j: halfspaces.back()){
                j=-j;
            }
        }
    }
    return getDistance(q, halfspaces);
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

float oru::single_query(level &idx, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree,
                        int k, int ret_size, vector<float>& q, fstream &log) {
    //vector<bool> filter(idx.idx[k].size(),false);
    //return kcell_filter(idx.idx[k],filter,ret_size,q,idx.dim);

    float init_dis=0.0001;
    while (true){
        if (init_dis>1.0) return 1.0; // k-level option is less than ret_size
        vector<float> ql,qu;
        ql.clear();qu.clear();
        for (int d=0;d<idx.dim-1;d++){
            ql.push_back(q[d]-init_dis);
            qu.push_back(q[d]+init_dis);
        }
        vector<int> kcellID;kcellID.clear();
        RangeQueryFromRtree(rt,ramTree,ql,qu,kcellID);
        unordered_set<int> results; results.clear();
        for (auto it=kcellID.begin();it!=kcellID.end();it++){
            for (auto p=idx.idx[k][*it].topk.begin();p!=idx.idx[k][*it].topk.end();p++){
                results.insert(*p);
            }
        }

        if (results.size()>=ret_size){
            vector<pair<float, int>> dis2q; dis2q.clear();
            for (int i=0;i<kcellID.size();i++){
                int id=kcellID[i];
                idx.UpdateH(idx.idx[k][id]);
                //float dis=GetDistance(q,idx.idx[k][id].r, idx.dim); // replace!!!
                vector<vector<float>> HS; HS.clear();
                for (auto it=idx.idx[k][id].r.H.begin();it!=idx.idx[k][id].r.H.end();it++){
                    HS.push_back(it->w);
                    if(it->side){
                        for(auto &j: HS.back()){
                            j=-j;
                        }
                    }
                }
                float dis=getDistance(q,HS);
                dis2q.push_back(make_pair(dis,id));
            }
            sort(dis2q.begin(),dis2q.end());
            unordered_set<int> ret;ret.clear();
            for (auto it=dis2q.begin();it!=dis2q.end();it++){
                int id=it->second;
                for (auto p=idx.idx[k][id].topk.begin();p!=idx.idx[k][id].topk.end();p++){
                    ret.insert(*p);
                }
                if (ret.size()>=ret_size) return it->first;
            }
            return 0.0;
        }
        else init_dis= init_dis*2;
    }
}


// --------------------------------- new code since 2021/7/8 -------------------------------------
float oru::single_query(level &idx, int k, int ret_size, vector<float>& q, fstream &log) {
    float oru_ret_dis=INFINITY;
    vector<vector<kcell>> &cells = idx.idx;
    if(cells.empty() || cells[0].empty()){
        return oru_ret_dis;
    }
    kcell &root=cells[0][0];
    multimap<double, kcell> heap;
    heap.emplace(INFINITY, root);
    unordered_set<int> ret_option;
    unordered_set<kcell*> added;
    while(!heap.empty() && ret_option.size()<ret_size){
        kcell nearest_cell=heap.begin()->second;
        oru_ret_dis=heap.begin()->first;
        heap.erase(heap.begin());
        if(nearest_cell.curk==k){
            for(auto &i:nearest_cell.topk){
                ret_option.insert(i);
            }
        }else{
            for(auto &i:nearest_cell.Next){
                kcell *child_cell=&cells[nearest_cell.curk+1][i];// cell
                if(added.find(child_cell)!=added.end()){
                    continue;
                }else{
                    added.insert(child_cell);
                }
                idx.UpdateH(*child_cell);
                vector<vector<float>> HS;
                for (auto & it : child_cell->r.H){
                    HS.push_back(it.w);
                    if(it.side){// this is importance !!!!!!!
                        for(auto &j: HS.back()){
                            j=-j;
                        }
                    }
                }
                double dis=getDistance(q,HS);
                heap.emplace(dis, *child_cell);
            }
            if(nearest_cell.curk<k && nearest_cell.Next.empty()){
                vector<kcell> NextCell;
                idx.SingleCellSplit(k, nearest_cell,NextCell);
                for (auto &child_cell:NextCell) {
                    idx.UpdateH(child_cell);
                    vector<vector<float>> HS;
                    for (auto & it : child_cell.r.H){
                        HS.push_back(it.w);
                        if(it.side){// this is importance !!!!!!!
                            for(auto &j: HS.back()){
                                j=-j;
                            }
                        }
                    }
                    double dis=getDistance(q,HS);
                    heap.emplace(dis, child_cell);
                }
            }
        }
    }
    return oru_ret_dis;
}

void oru::multiple_query(level &idx, int k, int ret_size, int q_num, fstream &log) {
    clock_t cur_time=clock();
    vector<vector<float>> q_list;
    generate_query(idx,q_num, q_list);
    for (int i=0;i<q_num;i++){
        float answer=single_query(idx, k, ret_size, q_list[i],log);
        cout << "The answer of oru query " << i << ": " << answer << endl;
        log << "The answer of oru query " << i << ": " << answer << endl;
    }
    cout << "Average oru query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average oru query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}