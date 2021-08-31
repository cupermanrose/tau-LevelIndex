
#ifndef K_LEVEL_UTILS_H
#define K_LEVEL_UTILS_H
#include <vector>
#include <string>
#include "rtree.h"
using std::unordered_set;
using std::unordered_map;
using std::vector;

void kskyband(
        vector<int> &ret,
        const vector<vector<float>> &data,
        const int k,
        bool rtree=true,
        const Rtree *rtree_rt= nullptr);

void kskyband_rtree(
        vector<int> &ret,
        const vector<vector<float>> &data,
        const int k,
        const Rtree *rtree_rt,
        unordered_map<long int, RtreeNode*> &ramTree);

void kskyband_nortree(
        vector<int> &ret,
        const vector<vector<float>> &data,
        const int k);

void onionlayer(vector<int> &ret, vector<int>& ret_layer, vector<int>& candidate, vector<vector<float>>& data, int k);

float GetScore(vector<float> &w, vector<float> &p, int dim); // w[dim-1]=1.0-sigma(w[0] to w[dim-2])


// Return true: oj dominates oi in this region;
template<typename POINT>
bool RegionDominate(vector<POINT> &V, vector<float> &oi, vector<float> &oj, int dim) {
    float s1, s2;
    bool flag=false;
    for (auto it=V.begin();it!=V.end();it++){
        s1=GetScore(*it,oi,dim);
        s2=GetScore(*it,oj,dim);
        if (s1>s2){
            return false;
        }else if(s1<s2){
            flag=true;
        }
    }
    return flag;
}

class ch{
    unordered_set<int> rest;
    unordered_map<int, int> pdtid_layer;
    unordered_map<int, int> dominated_cnt;// only use for build k-convex-hull

    vector<vector<int>> chs;

    vector<int> EMPTY;
    vector<vector<float>>& pointSet;
    int d;
public:
    unordered_set<int> all;
    vector<int> rskyband;
    unordered_map<int, vector<int>> A_p;
    unordered_map<int, unordered_set<int>> do_map;
    unordered_map<int, unordered_set<int>> dominated_map;
    ch(vector<int> &idxes, vector<vector<float>> &point_set, int dim);

    void fast_non_dominate_sort(
            const unordered_map<int, unordered_set<int>> &do_map,
            unordered_map<int, int>& dominated_cnt,
            const vector<int> &last_layer);

    void build_do_re(vector<int> &idxes, vector<vector<float>> &point_set, int dim);

    const vector<int>& get_next_layer();

    int get_option_layer(int option);

    const vector<int>& get_neighbor_vertex(int option);

    const vector<int>& get_layer(int which_layer);

    ~ch();
};

void build_onion(const std::string &s, int dim, int tau);
void read_onion(const std::string &filename, vector<vector<int>> &ret);


/*
void pivotRegion(vector<float>& R, vector<float>& pivot)
{
    int dim = R.size() / 2;
    pivot.clear();
    for (int i = 0; i < dim; i++)
        pivot.push_back((R[2 * i] + R[2 * i + 1]) / 2);
}

bool isRdominated(const int  dim, vector<float>& R, float focal[], vector<float>& entry, bool& fDe)
{
    // generate HP;
    vector<float> tmpHS;
    float entry_d = entry[dim - 1];
    float focal_d = focal[dim - 1];
    for (int d = 0; d < dim - 1; d++)
    {
        tmpHS.push_back((focal[d] - focal_d) - (entry[d] - entry_d));
    }
    tmpHS.push_back(entry_d - focal_d);


    // verify each vertic
    int posCount = 0;
    int negCount = 0;
    int totVertics = pow(2, dim - 1);
    for (int i = 0; i < totVertics; i++)
    {
        stringstream ss;
        string tmpS = bitset<MAXDIMEN>(i).to_string();
        tmpS = tmpS.substr(tmpS.size() - (dim - 1), dim - 1);

        float sum = 0;
        for (int si = 0; si < tmpS.size(); si++)
        {
            if (tmpS[si] == '0')
                sum += R[si * 2] * tmpHS[si];
            else if (tmpS[si] == '1')
                sum += R[si * 2 + 1] * tmpHS[si];
            else
                cout << "bug here!!!" << endl;
        }

        if (sum > tmpHS[dim - 1])
            negCount++;
        else
            posCount++;
    }

    // dominationship
    if (negCount == totVertics)
    {
        fDe = true;
        return true;
    }
    else if (posCount == totVertics)
    {
        fDe = false;
        return true;
    }
    else
        return false;
}

bool countRegionDominator(int dimen, float pt[], vector<long int>& rskyband, float* PG[], vector<float>& R, const int k, unordered_set<long int>& dominators)
{
    vector<float> record(dimen,0);
    dominators.clear();
    bool fDe;
    int count = 0;

    for (int i = 0; i < rskyband.size(); i++)
    {
        for (int di = 0; di < dimen; di++)
            record[di] = (PG[rskyband[i]][di] + PG[rskyband[i]][di + dimen]) / 2;
        fDe = false;
        if (isRdominated(dimen, R, pt, record, fDe))
        {
            if (fDe == false)
            {
                count++;
                dominators.insert(rskyband[i]);
            }
        }
    }
    if (count < k)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Naive_rskyband(vector<float>& region, const int dimen, Rtree& a_rtree, vector<long int>& rskyband, float* PG[], const int k)
{
    vector<float> pivot;
    pivotRegion(region, pivot);

    RtreeNode* node;
    priority_queue<pair<float, int>> heap;
    int NegPageid;

    float pt[10];
    float maxscore;
    int pageID;
    float tmpScore;
    unordered_set<long int> dominators;

    heap.push(make_pair(INFINITY, a_rtree.m_memory.m_rootPageID));

    while (!heap.empty())
    {
        tmpScore = heap.top().first;
        pageID = heap.top().second;
        heap.pop();

        if (pageID > MAXPAGEID)
        {
            for (int j = 0; j < dimen; j++)
                pt[j] = (PG[pageID - MAXPAGEID][j] + PG[pageID - MAXPAGEID][j + dimen]) / 2;
            if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
            {
                rskyband.push_back(pageID - MAXPAGEID);
                //daGraph[pageID - MAXPAGEID] = dominators;
            }
        }
        else
        {
            node = a_rtree.m_memory.loadPage(pageID);
            if (node->isLeaf())
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    for (int j = 0; j < dimen; j++)
                        pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;

                    if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
                    {
                        maxscore = orderScore(pivot, pt, dimen);
                        heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
                    }
                }
            }
            else
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    for (int j = 0; j < dimen; j++)
                    {
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                    }
                    if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
                    {
                        maxscore = orderScore(pivot, pt, dimen);
                        heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
                    }
                }
            }
        }
    }
}
*/

#endif //K_LEVEL_UTILS_H


