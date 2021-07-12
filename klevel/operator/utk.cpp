//
// Created by cuper on 2021/4/16.
//

#include "utk.h"

void utk::generate_query(level &idx, int q_num, float utk_side_length, vector<vector<float>>& q_list) {
    srand(0); // random seed
    q_list.clear();
    for (int i=0;i<q_num;i++){
        vector<float> v((idx.dim-1)*2,0.0);
        /*float res=rand()/RAND_MAX-utk_side_length*(idx.dim-1);
        if (res<0) res=1.0-utk_side_length*(idx.dim-1);*/
        int res=rand()%10000;
        cout << res << endl;
        for (int d=0;d<idx.dim-1;d++){
            int tmp=rand() % res;
            v[d*2] = (float) tmp / 10000.0;
            v[d*2+1] = v[d*2]+utk_side_length;
            res=res-tmp;
        }
        q_list.push_back(v);
    }
    return;
}

bool utk::isIn(vector<float>& v, vector<float>& Qregion, int dim){
    for (int i = 0; i < dim - 1; i++) {
        if ((v[i] >= Qregion[2 * i]) && (v[i] <= Qregion[2 * i + 1])) continue;
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

    /*AddQregion(Qregion,r,dim);
    if (lp_adapter::is_Feasible(r.H,r.innerPoint,dim)) return true;
    else return false;*/

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
    //TODO This might be wrong
    return false;
}

void utk::AddQregion(vector<float> &Qregion, region &r, int dim) {
    for (int i=0;i<dim-1;i++){
        halfspace tmp_HP; tmp_HP.w.clear();
        for (int j=0;j<dim;j++) tmp_HP.w.push_back(0.0);
        tmp_HP.w[i]=1.0;

        tmp_HP.w[dim-1]=Qregion[2*i]; tmp_HP.side=true;
        r.H.push_back(tmp_HP);

        tmp_HP.w[dim-1]=Qregion[2*i+1]; tmp_HP.side=false;
        r.H.push_back(tmp_HP);
    }
    return;
}

void utk::single_query(level &idx, int k, vector<float> &Qregion, fstream &log) {
    int visit=0;
    list<kcell> queue; queue={idx.idx[0][0]}; // only contains rootcell
    set<pair<int,int>> hash_set; hash_set.clear();
    set<int> results; results.clear();
    vector<kcell> NextCell;
    while (!queue.empty()){
        auto cur_cell=queue.front(); queue.pop_front(); visit++;
        if (Intersect(Qregion,cur_cell.r, idx.dim)){
            for (auto it=cur_cell.topk.begin();it!=cur_cell.topk.end();it++) results.insert(*it);
            if (cur_cell.curk<k){
                if ((cur_cell.curk<idx.ik)&&(cur_cell.Next.size()!=0)) { // non-best-performance
                    for (auto it=cur_cell.Next.begin();it!=cur_cell.Next.end();it++){
                        if (hash_set.find(make_pair(cur_cell.curk+1, *it))==hash_set.end()) {
                            queue.push_back(idx.idx[cur_cell.curk+1][*it]);
                            hash_set.insert(make_pair(cur_cell.curk+1,*it));
                            idx.UpdateH(queue.back()); // generate halfspace representation
                        }
                    }
                }
                else{ // for large k
                    idx.SingleCellSplit(k,cur_cell,NextCell);
                    for (auto it=NextCell.begin();it!=NextCell.end();it++){
                        queue.push_back(*it);
                    }
                }
            }
        }
    }
    NextCell.clear();
    vector<kcell>().swap(NextCell);
    cout << "Visiting cells of utk query: " <<  visit << endl;
    log << "Visiting cells of utk query: " <<  visit << endl;
    //cout << "Result cells of utk query: " <<  visit << endl;
    //log << "Result cells of utk query: " <<  visit << endl;
    cout << "Results of utk query: " << results.size() << endl;
    log << "Results of utk query: " << results.size() << endl;
    return;
}

void utk::multiple_query(level &idx, int k, int q_num, float utk_side_length, fstream &log) { // dag travel
    vector<vector<float>> q_list;
    generate_query(idx,q_num, utk_side_length, q_list);
    clock_t cur_time=clock();
    for (int i=0;i<q_num;i++){
        cout << "utk query " << i << ": " << endl;
        log << "utk query " << i << ": " << endl;
        single_query(idx,k,q_list[i],log);
    }
    cout << "Average utk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average utk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}

/*
void LocalFilter_utk(level& idx, kcell& cell, vector<int>&S1, vector<int>&Sk, int qk){
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
*/

/*
int utk::single_query(level &idx, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree,
                      int k, vector<float> &Qregion, fstream &log) {
    vector<float> ql,qu;
    ql.clear();qu.clear();
    for (int i=0;i<idx.dim-1;i++){
        ql.push_back(Qregion[i*2]);
        qu.push_back(Qregion[i*2+1]);
    }

    vector<int> kcellID;
    RangeQueryFromRtree(rt,ramTree,ql,qu,kcellID);
    unordered_set<int> results; results.clear();

    for (auto it=kcellID.begin();it!=kcellID.end();it++){
        if (Intersect(Qregion, idx.idx[k][*it].r,idx.dim)){
            for (auto p=idx.idx[k][*it].topk.begin();p!=idx.idx[k][*it].topk.end();p++){
                results.insert(*p);
            }
        }
    }
    return results.size();

}
*/

/*
void utk::SplitDFS(level& idx, kcell& cell, vector<float> &Qregion, int qk, unordered_set<int>& results) {
    if ((cell.curk>=qk)||(cell.r.V.size()==0)) return;
    for (auto it=cell.topk.begin();it!=cell.topk.end();it++) results.insert(*it);
    vector<int> S1,Sk;
    LocalFilter_utk(idx,cell,S1,Sk,qk);
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
            if (!Intersect(Qregion,newcell.r, idx.dim)) continue; // will not contribute to utk query
            else{
                int ave_vertex=0;
                idx.UpdateV(newcell, ave_vertex);
                SplitDFS(idx,newcell,Qregion, qk,results);
            }
        }
    }
    return;
}
*/

/*int utk::single_query_largek(level &idx, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree,
                             int k, vector<float> &Qregion, fstream &log) { // COMPUTE FROM LEVEL-0
    vector<float> ql,qu;
    ql.clear();qu.clear();
    for (int i=0;i<idx.dim-1;i++){
        ql.push_back(Qregion[i*2]);
        qu.push_back(Qregion[i*2+1]);
    }
    vector<int> kcellID;
    RangeQueryFromRtree(rt,ramTree,ql,qu,kcellID);

    unordered_set<int> results; results.clear();
    for (auto it=kcellID.begin();it!=kcellID.end();it++){
        SplitDFS(idx,idx.idx[idx.ik][*it],Qregion,k,results);
    }
    return results.size();


//    vector<int> S1,Sk;
//    int ave_S1=0,ave_Sk=0,ave_vertex=0;
//    vector<vector<kcell>> tmp; tmp.clear();
//    vector<kcell> init_level;  init_level.clear();
//     unordered_set<int> results; results.clear();
//    tmp.emplace_back(init_level);
//    for (int i=0;i<k-idx.ik;i++){
//        vector<kcell> this_level;  this_level.clear(); idx.region_map.clear();
//        for (auto cur_cell=tmp[i].begin(); cur_cell!=tmp[i].end(); cur_cell++){
//            if (!Intersect(Qregion,cur_cell->r, idx.dim)) continue; // will not contribute to utk query
//            for (auto p=cur_cell->topk.begin();p!=cur_cell->topk.end();p++){
//                results.insert(*p);
//            }
//            idx.LocalFilter(k, S1,Sk,*cur_cell,ave_S1,ave_Sk);
//            for (auto p=S1.begin();p!=S1.end();p++){
//                kcell newcell;
//                idx.CreateNewCell(*p, S1, Sk, *cur_cell, newcell);
//                if (!idx.VerifyDuplicate(newcell,this_level)){
//                    if (lp_adapter::is_Feasible(newcell.r.H,newcell.r.innerPoint,idx.dim)){
//                        this_level.emplace_back(newcell);
//                        idx.region_map.insert(make_pair(newcell.hash_value,this_level.size()-1));
//                    }
//                }
//            }
//        }
//        //Compute V for each cell
//        for (auto cur_cell=this_level.begin();cur_cell!=this_level.end();cur_cell++){
//            idx.UpdateH(*cur_cell);
//            idx.UpdateV(*cur_cell,ave_vertex);
//        }
//        tmp.emplace_back(this_level);
//    }
//
//    return results.size();
}*/

/*void utk::multiple_query(level &idx, int k, int q_num, float utk_side_length, fstream &log) {
    clock_t rtree_time = clock();
    Rtree *rt = nullptr;
    unordered_map<long int, RtreeNode *> ramTree;
    if (k<=idx.ik) BuildRtree(idx.idx[k], rt, ramTree, idx.dim);
    else BuildRtree(idx.idx[idx.ik], rt, ramTree, idx.dim);
    log << "R-tree from k-level building time: " << (clock() - rtree_time) / (float) CLOCKS_PER_SEC << endl;
    cout << "R-tree from k-level building time: " << (clock() - rtree_time) / (float) CLOCKS_PER_SEC << endl;

    clock_t cur_time=clock();
    vector<vector<float>> q_list;

    generate_query(idx,q_num, utk_side_length, q_list);
//   string queryfile="/home/kemingli/klevel/query/utk/utk_query_4d_l0.001.txt";
//    fstream queryout(queryfile, ios::in);
//    for (int i=0;i<q_num;i++){
//        float tmp;
//        for (int j=0;j<(idx.dim-1)*2;j++) {
//            queryout>>tmp;
//            cout<<tmp<<" ";
//        }
//        cout << endl;
//    }
//    return;
//    for (int i=0;i<q_num;i++){
//        for (int j=0;j<idx.dim-1;j++)
//            queryout << q_list[i][2*j] << ' ' << q_list[i][2*j+1] << ' ';
//        queryout << endl;
//    }
//    queryout.close();
//    return;

    for (int i=0;i<q_num;i++){
        vector<int> utk_results;
        int answer;
        if (k<=idx.ik) answer=single_query(idx, rt, ramTree, k,q_list[i],log);
        else answer=single_query_largek(idx, rt, ramTree, k,q_list[i],log);
        cout << "The answer of utk query " << i << ": " << answer << endl;
        log << "The answer of utk query " << i << ": " << answer << endl;
    }
    cout << "Average utk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << "Average utk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    return;
}*/

