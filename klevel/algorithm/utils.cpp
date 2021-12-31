#include "utils.h"
#include "filemem.h"
#include "hypercube.h"
#include "point.h"
#include "rentry.h"
#include "rnode.h"
#include "tgs.h"
#include "rtree_adapter.h"
#include <chrono>
#include "qhull_adapter.h"

using namespace std;

template<typename V1, typename V2>
double dot(const V1 &v1, const V2 &v2, size_t size) {
    double ret = 0;
    for (int i = 0; i < size; ++i) {
        ret += v1[i] * v2[i];
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
bool dominatedByK(const int dimen, const V &pt, const vector<int> &kskyband, const vector<vector<float>> &PG, int k) {
    // see if pt is dominated by k options of kskyband
    if (kskyband.empty())
        return false;

    int count = 0;
    for (long pid : kskyband) {
        bool dominated = true;
        for (int i = 0; i < dimen; i++) {
            if (PG[pid][i] < pt[i]) {
                dominated = false;
                break;
            }
        }
        if (dominated) {
            count++;
            if (count >= k) {
                return true;
            }
        }
    }
    return false;
}

template<typename V1, typename V2>
float dist(const V1 &v1, const V2 &v2, int dim) {
    float ret = 0;
    float tmp;
    for (int i = 0; i < dim; ++i) {
        tmp = (v1[i] - v2[i]);
        ret += tmp * tmp;
    }
    return sqrt(ret);
}

void test_rt2(Rtree& a_rtree, RtreeNode& node, set<int> &test, unordered_map<long int, RtreeNode *> &ramTree){
    if(node.isLeaf()){
        for (int i = 0; i <node.m_usedspace ; ++i) {
            assert(test.find(node.m_entry[i]->m_id+MAXPAGEID)==test.end());
            test.insert(node.m_entry[i]->m_id+MAXPAGEID);
        }
    }else{
        for (int i = 0; i <node.m_usedspace ; ++i) {
            assert(test.find(node.m_entry[i]->m_id)==test.end());
            test.insert(node.m_entry[i]->m_id);
            test_rt2(a_rtree, *ramTree[node.m_entry[i]->m_id], test, ramTree);
        }
    }
}


void kskyband(
        vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k,
        bool rtree,
        const Rtree *rtree_rt) {
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    if (rtree) {
        cout<<"using k-skyband with rtree"<<endl;
        if (rtree_rt) {
            unordered_map<long int, RtreeNode *> ramTree;
            rtreeRAM(*rtree_rt, ramTree);
            kskyband_rtree(ret, data, k, rtree_rt, ramTree);
        } else {
            Rtree *rtree = nullptr;
            unordered_map<long int, RtreeNode *> ramTree;
            build_rtree(rtree, ramTree, data);
            set<int> test;
            RtreeNode *node=ramTree[rtree->m_memory.m_rootPageID];
            test_rt2(*rtree, *node, test, ramTree);
            cout<<"finish init rtree"<<endl;
            kskyband_rtree(ret, data, k, rtree, ramTree);
            // TODO release mem of rtree
        }
    } else {
        cout<<"using k-skyband without rtree"<<endl;
        kskyband_nortree(ret, data, k);
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< "using "<<elapsed_seconds.count() <<" sec to apply k-skyband"<<endl;
}

void kskyband_nortree(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k) {
    /*
     * the k-skyband contains thoes records that are dominated by fewer than k others
     */
    vector<int> do_cnt(data.size(), 0);
    vector<int> k_cnt(k, 0);
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
            k_cnt[do_cnt[i]]+=1;
        }
    }
    cout<<"for k-skyband, # of options:"<<endl;
    for (int l = 0; l <k ; ++l) {
        cout<<"@k="<<l+1<<":num="<<k_cnt[l]<<endl;
    }
}

void kskyband_rtree(
        std::vector<int> &ret,
        const std::vector<std::vector<float>> &data,
        const int k,
        const Rtree *rtree_rt,
        unordered_map<long int, RtreeNode *> &ramTree) {
    if (data.empty()) {
        return;
    }
    int dimen = data[0].size();
    RtreeNode *node;
    multimap<float, int> heap;
    multimap<float, int>::iterator heapIter;
    vector<float> ORIGNIN(dimen, 1.0);
    float mindist;
    for (int i = 0; i < dimen; i++)
        ORIGNIN[i] = 1;

    int pageID;
    float dist_tmp;

    heap.emplace(INFINITY, rtree_rt->m_memory.m_rootPageID);
    while (!heap.empty()) {
        heapIter = heap.begin();
        dist_tmp = heapIter->first;
        pageID = heapIter->second;
        heap.erase(heapIter);

        if (pageID >= MAXPAGEID) {
            if (!dominatedByK(dimen, data[pageID - MAXPAGEID], ret, data, k)) {
                ret.push_back(pageID - MAXPAGEID);
            }
        } else {
            //node = a_rtree.m_memory.loadPage(pageID);
            node = ramTree[pageID];
            if (node->isLeaf()) {
                for (int i = 0; i < node->m_usedspace; i++) {
                    if (!dominatedByK(dimen, data[node->m_entry[i]->m_id], ret, data, k)) {
                        mindist = dist(data[node->m_entry[i]->m_id], ORIGNIN, dimen);
                        heap.emplace(mindist, node->m_entry[i]->m_id + MAXPAGEID);
                    }
                }
            } else {
                for (int i = 0; i < node->m_usedspace; i++) {
                    if (!dominatedByK(dimen, node->m_entry[i]->m_hc.getUpper(), ret, data, k)) {
                        mindist = dist(node->m_entry[i]->m_hc.getUpper(), ORIGNIN, dimen);
                        heap.emplace(mindist, node->m_entry[i]->m_id);
                    }
                }
            }
        }
    }
}


template<typename VVF>
void build_qhull(const vector<int> &opt_idxes, VVF &PG, vector<vector<double>> &square_vertexes, Qhull *q_ptr) {
    int dim = square_vertexes[0].size();
    square_vertexes.clear();
    square_vertexes.emplace_back(dim);
    string s = to_string(dim) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for (int opt_idx:opt_idxes) {
        for (int i = 0; i < dim; ++i) {
            s += to_string(PG[opt_idx][i]) + " ";
        }
    }
    for (vector<double> &square_vertex : square_vertexes) {
        for (float j : square_vertex) {
            s += to_string(j) + " ";
        }
    }
    istringstream is(s);
    RboxPoints rbox;
    rbox.appendPoints(is);
    q_ptr->runQhull(rbox, "QJ");
}

void onionlayer(vector<int> &ret, vector<int>& ret_layer ,vector<int> &candidate, vector<vector<float>> &data, int k) {
    if (data.empty()) {
        return;
    }
    int dim = data[0].size();
    ch c(candidate, data, dim);
    cout<<"for onion layer, # of options:"<<endl;
    for (int i = 1; i <= k; ++i) {
        for (int id: c.get_layer(i)) {
            ret.push_back(id);
            ret_layer.push_back(i);
        }
        cout << "#k=" << i << ";num=" << c.get_layer(i).size() << endl;
    }
    cout << "#total=" << ret.size() << endl;
}

float GetScore(vector<float> &w, vector<float> &p, int dim) { // w[dim-1]=1.0-sigma(w[0] to w[dim-2])
    float score = 0.0;
    float res = 1.0;
    for (int i = 0; i < dim - 1; i++) {
        score = score + w[i] * p[i];
        res = res - w[i];
    }
    score = score + res * p[dim - 1];
    return score;
}


ch::ch(vector<int> &idxes, vector<vector<float>> &point_set, int dim) : pointSet(point_set) {
    this->rskyband = idxes;
    this->rest = unordered_set<int>(idxes.begin(), idxes.end());
    this->all = unordered_set<int>(idxes.begin(), idxes.end());
//        this->pointSet=point_set;
    this->d = dim;
    build_do_re(idxes, point_set, dim);
}

void ch::fast_non_dominate_sort(
        const unordered_map<int, unordered_set<int>> &do_map,
        unordered_map<int, int> &dominated_cnt,
        const vector<int> &last_layer) {
    for (int opt:last_layer) {
        auto iter = do_map.find(opt);
        if (iter != do_map.end()) {
            for (int dominated:iter->second) {
                auto cnt_iter = dominated_cnt.find(dominated);
                if (cnt_iter != dominated_cnt.end()) {
                    cnt_iter->second -= 1;
                }
            }
        }
    }
}

void ch::build_do_re(vector<int> &idxes, vector<vector<float>> &point_set, int dim) {
    for (int i:idxes) {
        dominated_cnt[i] = 0;
        do_map[i] = unordered_set<int>();
        dominated_map[i] = unordered_set<int>();
    }
    for (int ii = 0; ii < idxes.size(); ++ii) {
        int i = idxes[ii];
        for (int ji = ii + 1; ji < idxes.size(); ++ji) {
            int j = idxes[ji];
            if (v1_dominate_v2(point_set[i], point_set[j], dim)) {
                do_map[i].insert(j);
                dominated_map[j].insert(i);
                ++dominated_cnt[j];
            } else if (v1_dominate_v2(point_set[j], point_set[i], dim)) {
                do_map[j].insert(i);
                dominated_map[i].insert(j);
                ++dominated_cnt[i];
//                }else{// non-dominate
            }
        }
    }
}

const vector<int> &ch::get_next_layer() {
    vector<vector<double>> square_vertexes(d + 1, vector<double>(d));
    if (!chs.empty()) {
        fast_non_dominate_sort(do_map, dominated_cnt, chs[chs.size() - 1]);
    }
    vector<int> rest_v;
    for (int i:rest) {
        auto iter = dominated_cnt.find(i);
        if (iter != dominated_cnt.end() && iter->second <= 0) {
            rest_v.push_back(i);
        }
    }
//        cout<<"no. of points to build convex hull: "<<rest_v.size()<<endl;
    vector<int> ch;
    if (rest_v.size() >= d + 1) {
        Qhull q;
        qhull_user qu;
        auto begin = chrono::steady_clock::now();
        build_qhull(rest_v, pointSet, square_vertexes, &q);
        auto end = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds = end - begin;
//            cout<<"finish build convex hull: "<<elapsed_seconds.count()<<endl;
        ch = qu.get_CH_pointID(q, rest_v);
        qu.get_neiVT_of_VT(q, rest_v, A_p);
    } else {
        for (int i:rest_v) {
            vector<int> tmp;
            for (int j:rest_v) {
                if (i != j) {
                    tmp.push_back(j);
                }
            }
            A_p[i] = tmp;
            ch.push_back(i);
        }
    }
    chs.push_back(ch);
    for (int idx:ch) {
        pdtid_layer[idx] = chs.size();
        rest.erase(idx);
    }

    return chs.back();
}

int ch::get_option_layer(int option) {
    auto iter = pdtid_layer.find(option);
    if (iter != pdtid_layer.end()) {
        return iter->second;
    } else {
        return -1; // not in current i layers
    }
}

const vector<int> &ch::get_neighbor_vertex(int option) {
//        assert(option>=0 && option <=objCnt);
    auto lazy_get = A_p.find(option);
    if (lazy_get != A_p.end()) {
        return lazy_get->second;
    } else {
        return EMPTY;
    }
}

const vector<int> &ch::get_layer(int which_layer) {
    // layer is starting from 1
    while (chs.size() < which_layer && !rest.empty()) {
        this->get_next_layer();
    }
    if (chs.size() < which_layer || which_layer <= 0) {
        return EMPTY;
    }
    return this->chs[which_layer - 1];
}

ch::~ch() {
}

void build_onion(const string &s, int dim, int tau){
//    int dim=7;
//    int tau=20; // NBA: tau=30
    string input=s+".dat";
    fstream fin(input, ios::in);
    vector<vector<float>> data;
    vector<float> cl(dim);
    vector<float> cu(dim);
    while (true) {
        int id;
        fin >> id;
        if (fin.eof())
            break;
        vector<float> tmp;
        for (int d = 0; d < dim; d++) fin >> cl[d];
        for (int d = 0; d < dim; d++) fin >> cu[d];
        for (int d = 0; d < dim; d++) tmp.push_back((cl[d]+cu[d])/2.0);

        data.push_back(tmp);
        if (data.size() % 1000 == 0)
            cout << ".";
        if (data.size() % 10000 == 0)
            cout << data.size() << " objects loaded" << endl;
    }

    cout << "Total number of objects: " << data.size() << endl;
    fin.close();
//    onionlayer(candidate_onionlayer, layer,  candidate_skyband,  OriginD, tau);
//    vector<int> in_idx(data.size());
//    iota(in_idx.begin(), in_idx.end(), 0);
    vector<int> in_idx;
    kskyband(in_idx, data,tau, false);

    if (data.empty()) {
        return;
    }
    ch c(in_idx, data, dim);
    fstream log;
    string logfile=s+".ch";
    log.open(logfile, ios::out);

    for (int i = 1; i <= tau; ++i) {
        log<<"#"<<i<<":";
        for (int id: c.get_layer(i)) {
            log<<id<<" ";
        }
        log<<endl;
    }
    log.close();
}

void read_onion(const string &filename, vector<vector<int>> &ret){
    FILE *in=fopen(filename.c_str(), "r");
    int cur;
    char UNUSED;
    int flag=1;
    vector<int> one_onion;
    while(flag){
        flag=fscanf(in, "%d", &cur);
        if(feof(in)){
            break;
        }
        if(flag){
            one_onion.push_back(cur);
        }else{
            flag=fscanf(in, "%c", &UNUSED);
            if(flag){
                if(UNUSED=='#'){
                    flag=fscanf(in, "%d", &cur);
                }
                flag=fscanf(in, "%c", &UNUSED);
                if(cur!=1){
                    ret.push_back(one_onion);
                    one_onion=vector<int>();
                }
            }// else eof
        }
    }
    fclose(in);
    ret.push_back(one_onion);
    cout<<ret.size()<<endl;
    for(auto &i:ret){
        cout<<i.size()<<endl;
    }
}

float utk_orderScore(vector<float>& pivot, vector<float>& entry){
    float ret=0;
    for (int i = 0; i < pivot.size(); i++){
        ret += pivot[i] * entry[i];
    }
    return ret;
}


bool r_dominate(vector<vector<float>>& vs, vector<float> &v1, vector<float>& v2){
    for (auto &v:vs) {
        vector<float> vc=v;
        float s=0;
        for (float val:vc) {
            s+=val;
        }
        vc.push_back(1.0-s);
        float s1=utk_orderScore(vc, v1);
        float s2=utk_orderScore(vc, v2);
        if(s1<s2){
            return false;
        }
    }
    return true;
}

bool utk_countRegionDominator(int dimen, vector<float> &pt, vector<int>& rskyband, vector<vector<float>> &PG, vector<vector<float>>& vs, const int k)
{
    // pt, dim
    // PG, dim
    // vs, region vertexes, dim-1
    if(rskyband.size()<k){
        return true;
    }
    vector<float> record(dimen,0);
    int count = 0;
    for (int i : rskyband){
        if(r_dominate(vs, PG[i], pt)) {
            count++;
            if (count >= k) {
                return false;
            }
        }
    }
    return count<k;
}

void utk_rskyband(vector<vector<float>>& region_v, const int dimen, Rtree& a_rtree, vector<int>& rskyband,
        vector<vector<float>>&PG, unordered_map<long int, RtreeNode *> &ramTree, vector<int>&topk, int k){
    // region_v: dim-1
    // dimen:  dim
    rskyband.clear();
    vector<float> pivot(dimen, 0);
    for (auto &v: region_v) {
        for (int i = 0; i <v.size() ; ++i) {
            pivot[i]+=v[i];
        }
    }
    for (float & l : pivot) {
        l/=region_v.size();
    }
    float one=1.0;
    for (float & l: pivot) {
        one-=l;
    }
    pivot[dimen-1]=one;
    RtreeNode* node;
    priority_queue<pair<float, int>> heap;
    int NegPageid;
    float maxscore;
    int pageID;
    float tmpScore;
    unordered_set<long int> dominators;
    heap.emplace(INFINITY, a_rtree.m_memory.m_rootPageID);
    vector<float> pt(dimen, 0);
    while (!heap.empty()){
        tmpScore = heap.top().first;
        pageID = heap.top().second;
        heap.pop();
        if (pageID >= MAXPAGEID){
            if (utk_countRegionDominator(dimen, PG[pageID - MAXPAGEID], rskyband, PG, region_v, k)){
                if(find(topk.begin(), topk.end(), pageID-MAXPAGEID)==topk.end()){
                    rskyband.push_back(pageID - MAXPAGEID);
                }
            }
        }else{
            node = ramTree[pageID];
            if (node->isLeaf()){
                for (int i = 0; i < node->m_usedspace; i++){
                    if (utk_countRegionDominator(dimen, PG[node->m_entry[i]->m_id], rskyband, PG, region_v, k)){
                        maxscore = utk_orderScore(pivot, PG[node->m_entry[i]->m_id]);
                        heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
                    }
                }
            }
            else{
                for (int i = 0; i < node->m_usedspace; i++){
                    for (int j = 0; j < dimen; j++){
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                    }
                    if (utk_countRegionDominator(dimen, pt, rskyband, PG, region_v, k)){
                        maxscore = utk_orderScore(pivot, pt);
                        heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
                    }
                }
            }
        }
    }
}

