//
// Created by cuper on 2021/4/1.
//

#include "dominateG.h"

int dominateG::offsetToID(vector<int> &offset, int dim) {
    int sum=0; int res=1;
    for (int i=1;i<=dim-1;i++){
        sum=sum+res*offset[i-1];
        res=res*div_num;
    }
    return sum;
}

void dominateG::GetCube(vector<point> &cube, vector<int> &offset, int pos, int dim) {
    float w_offset=1.0/div_num; // side length of the cube
    vector<float> w(dim,0);
    w[dim-1]=1.0;

    for (int i=0;i<dim-1;i++){
        int pos_i = pos >> i % 2;// pos_i = 0 means offset[i]; pos_i = 1 means offset[i]+1
        w[dim-1]=w[dim-1]-(offset[i]+pos_i)*w_offset;
        w[i]=(offset[i]+pos_i)*w_offset;
    }
    if (w[dim-1]<-EPS) return; // invalid vertex, some grid are not cubes
    point p; p.w=w;
    cube.emplace_back(p);
}

float dominateG::GetScore(vector<float> &w, vector<float> &p, int dim) { // w[dim-1]=1.0-sigma(w[0] to w[dim-2])
    float score=0.0;
    for (int i=0;i<dim;i++) score=score+w[i]*p[i];
    return score;
}

// return true if oi<oj in every point of V
bool dominateG::RegionDominate(vector<point> &V, vector<float> &oi, vector<float> &oj, int dim) {
    for (auto it=V.begin();it!=V.end();it++){
        if (GetScore(it->w,oi,dim)>GetScore(it->w,oj,dim)) return false;
    }
    return true;
}

void dominateG::GetGraph(vector<int> &offset, vector<vector<float>> &Allobj, int dim) {
    Gid=offsetToID(offset,dim);
    vector<point> cube; // (dim-1)-dimensional cube
    for (int i=0;i<1<<(dim-1);i++)
        GetCube(cube,offset, i, dim);

    G.clear();
    for (int i=0;i<Allobj.size();i++){
        set<int> tmp; tmp.clear();
        for (int j=0;j<Allobj.size();j++){
            if (i!=j){
                if (RegionDominate(cube, Allobj[i],Allobj[j],dim)) tmp.insert(j);
            }
        }
        G.insert(make_pair(i,tmp));
    }
}

void dominateG::EnumerateGrid(vector<int> &cur_offset, int cur_dim, int remains, int dim, vector<vector<float>>& Allobj, vector<dominateG>& Grid) {
    if (cur_dim==(dim-1)){
        if (remains<=0) return;
        dominateG G;
        G.GetGraph(cur_offset,Allobj,dim);
        Grid.emplace_back(G);
        return;
    }
    cur_dim++;
    for (int i=0;i<remains;i++){
        cur_offset[cur_dim-1]=i;
        EnumerateGrid(cur_offset,cur_dim,remains-i,dim,Allobj,Grid);
    }
}

int dominateG::FindCube(point &v,int dim) {
    float w_offset=1.0/div_num;
    vector<int> offset(dim,0);
    for (int d=0;d<dim-1;d++) offset[d]=v.w[d]/w_offset;
    return offsetToID(offset,dim);
}

void dominateG::MergeG(unordered_map<int, set<int>> &G, unordered_map<int, set<int>> &G_cube, set<int> &S) {

    if (G.size()==0){
        for (auto i=S.begin();i!=S.end();i++){
            set<int> edge; edge.clear();
            auto it=G_cube.find(*i);
            for (auto j=S.end();j!=S.end();j++){
                if (*i==*j) continue;
                if (it->second.find(*j)!=it->second.end()) edge.insert(*j);
            }
            G.insert(make_pair(*i,edge));
        }
    }
    else{
        for (auto i=G.begin();i!=G.end();i++){
            set<int> edge;edge.clear();
            auto j=G_cube.find(i->first);
            for (auto it=i->second.begin();it!=i->second.end();it++){
                if (j->second.find(*it)!=j->second.end()) edge.insert(*it);
            }
            i->second=edge;
        }
    }

    return;
}