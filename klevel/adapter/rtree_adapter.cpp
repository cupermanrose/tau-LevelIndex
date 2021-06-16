#include "rtree_adapter.h"

void release_rtree(Rtree *rt){}

void rtreeRAM(const Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree)
{
    ramTree.clear();
    queue<long int> H;
    RtreeNode* node;
    H.push(rt.m_memory.m_rootPageID);
    long int pageID;
    while (!H.empty())
    {
        pageID = H.front();
        H.pop();
        node = rt.m_memory.loadPage(pageID);
        ramTree[pageID] = node;
        if (node->isLeaf() == false)
        {
            for (int i = 0; i < node->m_usedspace; i++)
            {
                H.push(node->m_entry[i]->m_id);
            }
        }
    }
}

int countRecords(Rtree& a_rtree, int pageID)
{
    int sumRecords = 0;
    RtreeNode* node = a_rtree.m_memory.loadPage(pageID);
    if (node->isLeaf())
    {
        sumRecords = node->m_usedspace;
    }
    else
    {
        for (int i = 0; i < node->m_usedspace; i++)
        {
            sumRecords += countRecords(a_rtree, node->m_entry[i]->m_id);
        }
    }
    delete node;
    return sumRecords;
}

void aggregateRecords(Rtree& a_rtree)
{
    queue<long int> H;
    RtreeNode* node;
    H.push(a_rtree.m_memory.m_rootPageID);
    long int pageID;

    while (!H.empty())
    {
        pageID = H.front();
        H.pop();
        node = a_rtree.m_memory.loadPage(pageID);

        if (node->isLeaf() == false)
        {
            for (int i = 0; i < node->m_usedspace; i++)
            {
                node->m_entry[i]->num_records = countRecords(a_rtree, node->m_entry[i]->m_id);
                H.push(node->m_entry[i]->m_id);
            }
        }
        a_rtree.m_memory.writePage(pageID, node);
    }
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