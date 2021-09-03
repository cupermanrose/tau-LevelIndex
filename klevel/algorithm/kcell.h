
#ifndef K_LEVEL_KCELL_H
#define K_LEVEL_KCELL_H

#include "region.h"
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <boost/functional/hash.hpp>

using namespace std;

class kcell {
public:
    int curk, objID;
    size_t hash_value;
    region r;
    vector<int> topk; // the top-kth in this region
    vector<int> S1;
    vector<int> Stau; // top-tau candidates set
    vector<int> Next;

public:
    kcell();
    ~kcell();

    void WriteToDisk(ofstream& Outfile, bool OutStau);
    void ReadFromDisk(ifstream& Infile);
    void Get_HashValue();
    void TobeRoot(vector<int>& candidates, int dim);
    void FreeMem();
};


#endif //K_LEVEL_KCELL_H
