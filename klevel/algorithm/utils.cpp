#include "utils.h"
#include "filemem.h"
#include "hypercube.h"
#include "point.h"
#include "rentry.h"
#include "rnode.h"
#include "tgs.h"
#include "rtree_adapter.h"
#include <chrono>
using namespace std;

template<typename V1, typename V2>
double dot(const V1& v1, const V2& v2, size_t size){
    double ret=0;
    for (int i = 0; i <size ; ++i) {
        ret+=v1[i]*v2[i];
    }
    return ret;
}

template<typename V, typename U>
bool v1_dominate_v2(const V &v1, const U &v2, size_t size) {
    /*
     * /tpara V array, pointer
     */
    for (auto i = 0; i < size; ++i) {
        if (v1[i] < v2[i]) {
            return false;
        }
    }
    return true;
}

template<typename V>
bool v1_dominate_v2(const V &v1, const V &v2) {
    /*
     * /tpara V vector, vector
     */
    return v1_dominate_v2(v1, v2, v1.size());
}

template<typename V>
bool dominatedByK(const int dimen, const V &pt, const vector<int> &kskyband, const vector<vector<float>>&PG, int k)
{
    // see if pt is dominated by k options of kskyband
    if (kskyband.empty())
        return false;

    int count = 0;
    for (long pid : kskyband)
    {
        bool dominated = true;
        for (int i = 0; i < dimen; i++)
        {
            if (PG[pid][i] < pt[i])
            {
                dominated = false;
                break;
            }
        }
        if (dominated) {
            count++;
            if(count>=k){
                return true;
            }
        }
    }
    return false;
}

template<typename V, typename U>
bool v1_rdominate_v2(const V &v1, const U &v2, size_t size,
                     const std::vector<std::vector<float>> &region) {
    /*
     * /tpara V array, pointer
     */
//    for (auto &v:region) {
//        if (v1[i] < v2[i]) {
//            return false;
//        }
//    }
//TODO
    return true;
}

template<typename V>
bool v1_rdominate_v2(const V &v1, const V &v2,
                     const std::vector<std::vector<float>> &region) {
    /*
     * /tpara V vector, vector
     */
    //TODO
    return v1_dominate_v2(v1, v2, v1.size());
}

template<typename V>
bool rdominatedByK(const int dimen, const V &pt, const vector<int> &kskyband, const vector<vector<float>>&PG, int k,
                   const std::vector<std::vector<float>> &region)
{
    //TODO

    // see if pt is dominated by k options of kskyband
    if (kskyband.empty())
        return false;

    int count = 0;
    for (long pid : kskyband)
    {
        bool dominated = true;
        for (int i = 0; i < dimen; i++)
        {
            if (PG[pid][i] < pt[i])
            {
                dominated = false;
                break;
            }
        }
        if (dominated) {
            count++;
            if(count>=k){
                return true;
            }
        }
    }
    return false;
}


template<typename V1, typename V2>
float dist(const V1& v1, const V2& v2, int dim){
    float ret=0;
    float tmp;
    for (int i = 0; i < dim; ++i) {
        tmp=(v1[i]-v2[i]);
        ret+=tmp*tmp;
    }
    return sqrt(ret);
}

//class ch{
//    unordered_set<int> rest;
//    unordered_map<int, int> pdtid_layer;
//    unordered_map<int, int> dominated_cnt;// only use for build k-convex-hull
//
//    vector<vector<int>> chs;
//
//    vector<int> EMPTY;
//    float** pointSet;
//    int d;
//public:
//    unordered_set<int> all;
//    vector<int> rskyband;
//    unordered_map<int, vector<int>> A_p;
//    unordered_map<int, unordered_set<int>> do_map;
//    unordered_map<int, unordered_set<int>> dominated_map;
//
//    ch(vector<int> &idxes, float** &point_set, int dim){
//        this->rskyband=idxes;
//        this->rest=unordered_set<int>(idxes.begin(), idxes.end());
//        this->all=unordered_set<int>(idxes.begin(), idxes.end());
//        this->pointSet=point_set;
//        this->d=dim;
//        build_do_re(idxes, point_set, dim);
//    }
//
//    void fast_non_dominate_sort(
//            const unordered_map<int, unordered_set<int>> &do_map,
//            unordered_map<int, int>& dominated_cnt,
//            const vector<int> &last_layer){
//        for (int opt:last_layer) {
//            auto iter=do_map.find(opt);
//            if(iter!=do_map.end()){
//                for(int dominated:iter->second){
//                    auto cnt_iter=dominated_cnt.find(dominated);
//                    if(cnt_iter!=dominated_cnt.end()){
//                        cnt_iter->second-=1;
//                    }
//                }
//            }
//        }
//    }
//
//    void build_do_re(vector<int> &idxes, float** &point_set, int dim){
//        for (int i:idxes){
//            dominated_cnt[i]=0;
//            do_map[i]=unordered_set<int>();
//            dominated_map[i]=unordered_set<int>();
//        }
//        for (int ii = 0;ii<idxes.size();++ii) {
//            int i=idxes[ii];
//            for (int ji = ii+1; ji <idxes.size() ; ++ji) {
//                int j=idxes[ji];
//                if(v1_dominate_v2(point_set[i], point_set[j], dim)){
//                    do_map[i].insert(j);
//                    dominated_map[j].insert(i);
//                    ++dominated_cnt[j];
//                }else if(v1_dominate_v2(point_set[j], point_set[i], dim)){
//                    do_map[j].insert(i);
//                    dominated_map[i].insert(j);
//                    ++dominated_cnt[i];
////                }else{// non-v1_dominate_v2
//                }
//            }
//        }
//    }
//
//    const vector<int>& get_next_layer(){
//        vector<vector<double>> square_vertexes(d+1, vector<double>(d));
//        if(!chs.empty()){
//            fast_non_dominate_sort(do_map, dominated_cnt, chs[chs.size()-1]);
//        }
//        vector<int> rest_v;
//        for(int i:rest){
//            auto iter=dominated_cnt.find(i);
//            if(iter!=dominated_cnt.end() && iter->second<=0){
//                rest_v.push_back(i);
//            }
//        }
//        cout<<"no. of points to build convex hull: "<<rest_v.size()<<endl;
//        vector<int> ch;
//        if(rest_v.size()>=d+1){
//            Qhull q;
//            qhull_user qu;
//            auto begin=chrono::steady_clock::now();
//            build_qhull(rest_v, pointSet, square_vertexes, &q);
//            auto end=chrono::steady_clock::now();
//            chrono::duration<double> elapsed_seconds= end-begin;
//            cout<<"finish build convex hull: "<<elapsed_seconds.count()<<endl;
//            ch=qu.get_CH_pointID(q, rest_v);
//            qu.get_neiVT_of_VT(q, rest_v, A_p);
//        }else{
//            for(int i:rest_v){
//                vector<int> tmp;
//                for(int j:rest_v){
//                    if(i!=j){
//                        tmp.push_back(j);
//                    }
//                }
//                A_p[i]=tmp;
//                ch.push_back(i);
//            }
//        }
//        chs.push_back(ch);
//        for (int idx:ch) {
//            pdtid_layer[idx] =  chs.size();
//            rest.erase(idx);
//        }
//
//        return chs.back();
//    }
//
//    int get_option_layer(int option){
//        auto iter=pdtid_layer.find(option);
//        if(iter!=pdtid_layer.end()){
//            return iter->second;
//        }else{
//            return -1; // not in current i layers
//        }
//    }
//
//    const vector<int>& get_neighbor_vertex(int option){
////        assert(option>=0 && option <=objCnt);
//        auto lazy_get=A_p.find(option);
//        if(lazy_get!=A_p.end()){
//            return lazy_get->second;
//        }else{
//            return EMPTY;
//        }
//    }
//
//    const vector<int>& get_layer(int which_layer){
//        // layer is starting from 1
//        while(chs.size()<which_layer && !rest.empty()){
//            this->get_next_layer();
//        }
//        if(chs.size()<which_layer || which_layer<=0){
//            return EMPTY;
//        }
//        return this->chs[which_layer-1];
//    }
//
//    ~ch(){
//    }
//};

void kskyband(
        vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k,
        bool rtree,
        const Rtree *rtree_rt){
    if(rtree){
        if(rtree_rt){
            unordered_map<long int, RtreeNode*> ramTree;
            rtreeRAM(*rtree_rt, ramTree);
            kskyband_rtree(ret,data,k,rtree_rt, ramTree);
        }else{
            Rtree* rtree = nullptr;
            unordered_map<long int, RtreeNode*> ramTree;
            build_rtree(rtree, ramTree, data);
            kskyband_rtree(ret,data,k,rtree, ramTree);
            // TODO release mem of rtree
        }
    }else{
        kskyband_nortree(ret,data,k);
    }

}

void kskyband_nortree(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k) {
    /*
     * the k-skyband contains thoes records that are dominated by fewer than k others
     */
    vector<int> do_cnt(data.size(), 0);
    for (auto i = 0; i < data.size(); ++i) {
        for (auto j = i + 1; j < data.size(); ++j) {
            if (do_cnt[i] >= k) {
                break;
            }
            if (v1_dominate_v2(data[i], data[j])) {
                ++do_cnt[j];
            } else if (v1_dominate_v2(data[j], data[i])) {
                ++do_cnt[i];
            }
        }
        if (do_cnt[i] < k) {
            ret.push_back(i);
        }
    }
}

void kskyband_rtree(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k,
        const Rtree *rtree_rt,
        unordered_map<long int, RtreeNode*>& ramTree){
    if(data.empty()){
        return;
    }
    int dimen=data[0].size();
    RtreeNode* node;
    multimap<float, int> heap;
    multimap<float, int>::iterator heapIter;
    vector<float> ORIGNIN(dimen, 1.0);
    float mindist;
    for (int i = 0; i < dimen; i++)
        ORIGNIN[i] = 1;

    int pageID;
    float dist_tmp;

    heap.emplace(INFINITY, rtree_rt->m_memory.m_rootPageID);

    while (!heap.empty())
    {
        heapIter = heap.begin();
        dist_tmp = heapIter->first;
        pageID = heapIter->second;
        heap.erase(heapIter);

        if (pageID > MAXPAGEID)
        {
            if (!dominatedByK(dimen, data[pageID - MAXPAGEID], ret, data, k))
            {
                ret.push_back(pageID - MAXPAGEID);
            }
        }
        else
        {
            //node = a_rtree.m_memory.loadPage(pageID);
            node = ramTree[pageID];
            if (node->isLeaf())
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    if (!dominatedByK(dimen, data[node->m_entry[i]->m_id], ret, data, k))
                    {
                        mindist = dist(data[node->m_entry[i]->m_id], ORIGNIN, dimen);
                        heap.emplace(mindist, node->m_entry[i]->m_id + MAXPAGEID);
                    }
                }
            }
            else
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    if (!dominatedByK(dimen, node->m_entry[i]->m_hc.getUpper(), ret, data, k))
                    {
                        mindist =dist(node->m_entry[i]->m_hc.getUpper(), ORIGNIN, dimen);
                        heap.emplace(mindist, node->m_entry[i]->m_id);
                    }
                }
            }
        }
    }
}
