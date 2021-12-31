
#include "oru.h"
#include "qp_adapter.h"
void oru::generate_query(level &idx, int q_num, vector<vector<float>>& q_list) {
    srand(0); // random seed
    q_list.clear();
    for (int i=0;i<q_num;i++){ // generate user preference uniform under the constraint \sum v_i=1,
        // to understand this code, you may need some knowledge of definite in
        vector<float> v(idx.dim-1,0.0);
        float res=1.0;
        for (int d=0;d<idx.dim-1;d++){
            v[d] = res*(1.0-pow((double)random()/RAND_MAX,  1.0/(idx.dim-d-1)));
            res-=v[d];
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

inline size_t oru_non_order_butk_hash(const vector<int> &input, int exactK){
    vector<int> tmp=input;
    sort(tmp.begin(),tmp.end()); // make it in order to output the unique hash_value for a set
    size_t seed = 0;
    for (int & it : tmp){
        if(it!=exactK){
            boost::hash_combine(seed, it);
        }
    }
    boost::hash_combine(seed, exactK);
    return seed;
}

ostream & operator<<(ostream &out, vector<int> &v){
    if (v.empty()) {
        return out;
    }
    out << v[0];
    for (auto i = 1; i < v.size(); ++i) {
        out << ", " << v[i];
    }
    return out;
}
double largeKTime=0.0;
// --------------------------------- new code since 2021/7/8 -------------------------------------
float oru::single_query(level &idx, int k, int ret_size, vector<float>& q, fstream &log,
    int &push_cnt, int& pop_cnt) {
    float oru_ret_dis=INFINITY;
    vector<vector<kcell>> &cells = idx.idx;
    unordered_set<int> ret_option;
    if(cells.empty() || cells[0].empty()){
        return oru_ret_dis;
    }
    if(k>=ret_size){
        multimap<double, int> topk;
        for (int i=0; i<idx.Allobj.size(); ++i) {
            double s=GetScore(q, idx.Allobj[i], idx.dim);
            topk.emplace(s, i);
            if(topk.size()>=k){
                topk.erase(topk.begin());
            }
        }
        oru_ret_dis=topk.rbegin()->first;
        // the comment below is for future use
//        for (auto &s_id: topk) {
//            ret_option.emplace(s_id.second);
//        }
        return oru_ret_dis;
    }
    kcell &root=cells[0][0];
    multimap<double, kcell*> heap;
    heap.emplace(INFINITY, &root);
//    unordered_set<kcell*> added;
    unordered_set<size_t> nr_added;
    vector<kcell*> toBeFree;
    push_cnt=0;
    pop_cnt=0;
    while(!heap.empty() && ret_option.size()<ret_size){
        kcell *nearest_cell=heap.begin()->second;
        oru_ret_dis=heap.begin()->first;
        heap.erase(heap.begin());
        ++pop_cnt;
        for(auto &i:nearest_cell->topk){
            ret_option.insert(i);
        }
        if(nearest_cell->curk<k){
            clock_t b=clock();
            if(nearest_cell->Next.empty() || nearest_cell->curk>=idx.ik){
//                vector<kcell> NextCell;
//                idx.SingleCellSplit(k, nearest_cell,NextCell);
                vector<int> S1,Sk;
                int ave_S1=0,ave_Sk=0;
                idx.LocalFilter(k, S1,Sk,*nearest_cell,ave_S1,ave_Sk);
                for (auto p=S1.begin();p!=S1.end();p++){
                    if (idx.global_layer[*p]>nearest_cell->curk+1) continue;
                    auto *newcell=new kcell();
                    idx.CreateNewCell(*p,S1,Sk, *nearest_cell, *newcell);
                    std::size_t hashv=oru_non_order_butk_hash(newcell->topk, newcell->objID);
                    if(nr_added.find(hashv)!=nr_added.end()){
                        delete (newcell);
                        continue;
                    }
                    idx.UpdateH(*newcell);
                    vector<vector<float>> HS;
                    for (auto & it : newcell->r.H){
                        HS.push_back(it.w);
                        if(it.side){// this is importance !!!!!!!
                            for(auto &j: HS.back()){
                                j=-j;
                            }
                        }
                    }
                    double dis=getDistance(q,HS);
                    int UNUSED;
                    if(dis!=INFINITY){
                        idx.UpdateV(*newcell,UNUSED);
                        nr_added.insert(hashv);
                        heap.emplace(dis, newcell);
                        ++push_cnt;
                        toBeFree.emplace_back(newcell);
                    }else{
                        delete (newcell);
                    }

                }

            }else{
                for(auto &i:nearest_cell->Next){
                    assert(nearest_cell->curk<=idx.ik);
                    kcell *child_cell=&cells[nearest_cell->curk+1][i];// cell
                    std::size_t hashv=oru_non_order_butk_hash(child_cell->topk, child_cell->objID);
                    if(nr_added.find(hashv)!=nr_added.end()){
                        continue;
                    }else{
                        nr_added.insert(hashv);
                    }
                    idx.UpdateH(*child_cell);
                    int UNUSED;
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
                    heap.emplace(dis, child_cell);
                    ++push_cnt;
                }

            }

            clock_t e=clock();
            if(nearest_cell->curk>=idx.ik){
                largeKTime+=(e - b) / (float)CLOCKS_PER_SEC;
            }
        }
    }
    clock_t b=clock();
    for(auto &i:toBeFree){
        delete(i);
    }
    clock_t e=clock();
    largeKTime+=(e - b) / (float)CLOCKS_PER_SEC;

    return oru_ret_dis;
}

void oru::multiple_query(level &idx, int k, int ret_size, int q_num, fstream &log) {
    clock_t cur_time=clock();
    vector<vector<float>> q_list;
    generate_query(idx,q_num, q_list);
    int tt_pop_cnt=0, tt_push_cnt=0;
    largeKTime=0.0;
    for (int i=0;i<q_num;i++){
        clock_t qb=clock();
        int pop_cnt=0, push_cnt=0;
        float answer=single_query(idx, k, ret_size, q_list[i],log, push_cnt, pop_cnt);
        clock_t qe=clock();
        tt_push_cnt+=push_cnt;
        tt_pop_cnt+=pop_cnt;
        cout << "The answer of oru query " << i << ": " << answer << endl;
        log << "The answer of oru query " << i << ": " << answer << endl;
        cout << "# of visit cells " << i << ": " << push_cnt << endl;
        log << "# of visit cells " << i << ": " << push_cnt << endl;
        cout << "# of result cells " << i << ": " << pop_cnt << endl;
        log << "# of result cells " << i << ": " << pop_cnt << endl;
        cout << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
        log << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
    }
    if(tt_push_cnt!=0 && tt_pop_cnt!=0){
        cout << "Average # of oru visit cells: " << (double)tt_push_cnt/q_num << endl;
        log << "Average # of oru visit cells: " << (double)tt_push_cnt/q_num << endl;
        cout << "Average # of oru result cells: " << (double)tt_pop_cnt/q_num << endl;
        log << "Average # of oru result cells: " << (double)tt_pop_cnt/q_num << endl;
    }
    cout << "Average oru query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average oru query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    cout << "Large k oru query time: " << largeKTime / (float) q_num << endl;
    log << "Large k oru query time: " << largeKTime / (float) q_num << endl;
    return;
}