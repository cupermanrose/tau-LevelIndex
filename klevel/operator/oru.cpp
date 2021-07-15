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
    //if (isIn(q,r.H,dim)) return 0.0;
//    float dis=1.0;
//    for (auto it=r.V.begin();it!=r.V.end();it++){
//        float tmp=0.0, resq=1.0,resv=1.0;
//        for (int d=0;d<dim;d++){
//            tmp=tmp+(q[d]-(*it)[d])*(q[d]-(*it)[d]);
//            resq=resq-q[d];
//            resv=resv-(*it)[d];
//        }
//        tmp=sqrt(tmp+(resq-resv)*(resq-resv));
//        if (tmp<dis) dis=tmp;
//    }
//    return dis;
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

/*
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
*/

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