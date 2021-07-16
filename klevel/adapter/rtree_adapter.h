#include "filemem.h"
#include "hypercube.h"
#include "point.h"
#include "rentry.h"
#include "rnode.h"
#include "tgs.h"
#include "rtree.h"

void release_rtree(Rtree *rt);

void rtreeRAM(const Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);

int countRecords(Rtree& a_rtree, int pageID);

void aggregateRecords(Rtree& a_rtree);

template<typename VVF>
void build_rtree(Rtree* &rtree_rt, unordered_map<long int, RtreeNode*>& ramTree,VVF &data){
    if(data.empty()){
        return;
    }
    int dim=data[0].size();
    RtreeNodeEntry** p = new RtreeNodeEntry*[data.size()];
    for (int id=0;id<data.size();++id)
    {
        float* cl = new float[dim];
        float* cu = new float[dim];
        for (int i = 0; i <dim ; ++i) {
            cl[i]=data[id][i]-SIDELEN;
            cu[i]=data[id][i]+SIDELEN;
        }
        Hypercube hc(dim, cl, cu);
        p[id] = new RtreeNodeEntry(id, hc);
    }
    // build rtree
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    //FileMemory mem(PAGESIZE, "./result/index.txt", RtreeNodeEntry::fromMem, true);
    string indexf_name="index"+to_string(time(nullptr)%1000)+"txt"; // TODO specific by user
    FileMemory *mem=new FileMemory(PAGESIZE, indexf_name.c_str(), RtreeNodeEntry::fromMem, true);
    rtree_rt = TGS::bulkload(*mem, dim, maxChild, maxChild, (int)maxChild*0.3, (int)maxChild*0.3, p, data.size(), false);
//            cout << "[Rtree build done]" << endl;

    // in-memory rtree
//            cout << "cache R-tree into memory" << endl;
    rtreeRAM(*rtree_rt, ramTree);

    // aggregate rtree
    aggregateRecords(*rtree_rt);
}


template<typename VVF>
void box2rtree(Rtree* &rtree_rt, unordered_map<long int, RtreeNode*>& ramTree,VVF &data){
    if(data.empty()){
        return;
    }
    int dim=data[0].size()/2;
    RtreeNodeEntry** p = new RtreeNodeEntry*[data.size()];
    for (int id=0;id<data.size();++id)
    {
        Hypercube hc(dim, &(data[id][0]), &(data[id][dim]));
        p[id] = new RtreeNodeEntry(id, hc);
    }
    // build rtree
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    //FileMemory mem(PAGESIZE, "./result/index.txt", RtreeNodeEntry::fromMem, true);
    //TODO not using a file, if running multiple rtree, this would cause problem
    string indexf_name="index"+to_string(time(nullptr)%1000)+"txt"; // TODO specific by user
    FileMemory *mem=new FileMemory(PAGESIZE, indexf_name.c_str(), RtreeNodeEntry::fromMem, true);
    rtree_rt = TGS::bulkload(*mem, dim, maxChild, maxChild, (int)maxChild*0.3, (int)maxChild*0.3, p, data.size(), false);
//            cout << "[Rtree build done]" << endl;

    // in-memory rtree
//            cout << "cache R-tree into memory" << endl;
    rtreeRAM(*rtree_rt, ramTree);

    // aggregate rtree
    aggregateRecords(*rtree_rt);
}


// Apply R-tree in query processing
template<typename PTS>
bool boxIntersection(PTS &boxl, PTS &boxu,
                            vector<float> &targetl, vector<float> &targetu){
    // to satisfy:
    //     x_i>=boxl_i && x_i>=targetl_i && x_i<=boxu_i && x_i<=targetu_i
    int dim=targetl.size();
    for (int i = 0; i < dim-1; ++i) {
        if(max(boxl[i], targetl[i]) > min(boxu[i], targetu[i])){
            return false;
        }
    }
    return true;
}

void rtree_boxInter(vector<int> &ret, const Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree,
                    vector<float> &ql, vector<float> &qu);
