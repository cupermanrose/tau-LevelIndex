//
// Created by cuper on 2021/3/30.
//

#include "k_level_lib.h"
#include "rtree_adapter.h"
void BuildIndex(level& idx, string datafile, fstream& log, string idxfile) {
    ofstream idxout(idxfile);
    idx.LoadData(datafile);
    idx.Build(log,idxout);
    log.close();
    idxout.close();
}

void IncBuildIndex(level& idx, string datafile, fstream& log, string idxfile) {
    ofstream idxout(idxfile);
    idx.LoadData(datafile);
    idx.IncBuild(log,idxout);
    log.close();
    idxout.close();
}

void DFSBuildIndex(level& idx, string datafile, fstream& log, string idxfile) {
    ofstream idxout(idxfile);
    idx.LoadData(datafile);
    idx.DFSBuild(log,idxout);
    log.close();
    idxout.close();
}

void LoadIndex(level& idx, string datafile, fstream& log, string idxfile) {
    ifstream idxin(idxfile);
    idx.LoadData(datafile);
    vector<int> candidate;
    idx.GlobalFilter(log,candidate); //TODO remove this and load from index file

    vector<kcell> empty_level; empty_level.clear();
    idx.idx.emplace_back(empty_level);
    for (int k=1;k<=idx.tau;k++){
        if (k>idx.ik) break;
        idx.idx.push_back(empty_level);
        idx.ReadFromDisk(k,idxin);
        cout << "The size of level " << k << ": "<<idx.idx[k].size() << endl;
        log << "The size of level " << k << ": "<<idx.idx[k].size() << endl;
    }
    idxin.close();

    // Build R-tree
    /*
     BuildRtree(vector<kcell> L, rtree* rt){
        for (int i=0;i<idx.idx[idx.ik].size();i++){
            float* cl = new float[idx.dim];
            float* cu = new float[idx.dim];
            for (int d = 0;  d<idx.dim ; d++) {
                for (each v in idx.idx[i]) {
                    cl[d] = min(cl[d], v[d]);
                    cu[d] = max(cu[d], v[d]);
                }
            }
            Hypercube hc(dim, cl, cu);
            p[i] = new RtreeNodeEntry(i, hc); // i is kcell id
        }
     }
     RangeQuery(rtree* rt, float* ql, float* qu, vector<int> ids){
     }*/
    cout << "Load Index Done!" << endl;
    log << "Load Index Done!" << endl;
    idx.print_system_info(log);
}

inline vector<vector<float>> vertex2BOX(const vector<kcell> &L){
    // the convex hull defined by vertexes to lower bound upper bound box
    // dim-1 vertex to dim box
    if(L.empty() || L.begin()->r.V.empty()){
        return vector<vector<float>>();
    }
    int dim = L.begin()->r.V.size()+1;
    vector<vector<float>> ret;
    for(auto &iter: L){
        vector<float> box(dim*2);
        for (int j = 0; j < dim; ++j) {
            box[j]=1;
        }
        float bias=1;
        for(auto &vertex: iter.r.V){

            for(auto &attr: vertex){
                bias-=attr;
            }

            for (int i = 0; i < vertex.size(); ++i) {
                box[i]=min(box[i], vertex[i]);
            }
            box[dim-1]=min(box[dim-1], bias);

            for (int i = 0; i < vertex.size(); ++i) {
                box[i+dim]=max(box[i+dim], vertex[i]);
            }
            box[dim-1+dim]=max(box[dim-1+dim], bias);
        }
        ret.emplace_back(box);
    }
    return ret;
}

inline void BuildRtree(const vector<kcell> &L, Rtree* rt, unordered_map<long int, RtreeNode*>& ramTree){
    auto boxes=vertex2BOX(L);
    box2rtree(rt, ramTree, boxes);
}

template<typename PTS>
inline bool boxIntersection(PTS &boxl, PTS &boxu,
                            vector<float> &targetl, vector<float> &targetu){
    // to satisfy:
    //     x_i>=boxl_i && x_i>=targetl_i && x_i<=boxu_i && x_i<=targetu_i
    int dim=targetl.size();
    for (int i = 0; i < dim; ++i) {
        if(max(boxl[i], targetl[i]) > min(boxu[i], targetu[i])){
            return false;
        }
    }
    return true;
}

void rtree_boxInter(vector<int> &ret, const Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree,
                    vector<float> &ql, vector<float> &qu){
    multiset<int> heap;
    multiset<int>::iterator heapIter;
    RtreeNode *node;
    int pageID;
    heap.emplace(rtree_rt->m_memory.m_rootPageID);
    while (!heap.empty()) {
        heapIter = heap.begin();
        pageID = *heapIter;
        heap.erase(heapIter);
        node = ramTree[pageID];
        if (node->isLeaf()) {
            for (int i = 0; i < node->m_usedspace; i++) {
                if (boxIntersection(node->m_entry[i]->m_hc.getLower(), node->m_entry[i]->m_hc.getUpper(),
                                    ql, qu)) {
                    ret.push_back(node->m_entry[i]->m_id);
                }
            }
        } else {
            for (int i = 0; i < node->m_usedspace; i++) {
                if (boxIntersection(node->m_entry[i]->m_hc.getLower(), node->m_entry[i]->m_hc.getUpper(),
                                    ql, qu)) {
                    heap.emplace(node->m_entry[i]->m_id);
                }
            }
        }
    }
}

inline float sum(vector<float> &v){
    float ret=0;
    for(auto &f:v){
        ret+=f;
    }
    return ret;
}

void RangeQuery(const vector<kcell> &L, vector<float> &ql, vector<float> &qu, vector<int> &ret_ids){
    Rtree* rt=nullptr;
    unordered_map<long int, RtreeNode*> ramTree;
    BuildRtree(L, rt, ramTree);
    vector<float> target_l(ql);
    vector<float> target_u(qu);

    target_l.push_back(1-sum(ql));
    target_u.push_back(1-sum(qu));

    rtree_boxInter(ret_ids, rt, ramTree, target_l, target_u);
}

void RangeQueryFromRtree(Rtree* rt, unordered_map<long int, RtreeNode*> ramTree,
        vector<float> &ql, vector<float> &qu, vector<int> &ret_ids){
    vector<float> target_l(ql);
    vector<float> target_u(qu);

    target_l.push_back(1-sum(ql));
    target_u.push_back(1-sum(qu));

    rtree_boxInter(ret_ids, rt, ramTree, target_l, target_u);
}

