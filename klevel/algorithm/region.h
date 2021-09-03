

#ifndef K_LEVEL_REGION_H
#define K_LEVEL_REGION_H

#include <vector>
#include <fstream>
#include "lp_adapter.h"
#include "qhull_adapter.h"

using namespace std;

class region {
public:

    vector<halfspace> H; // halfspace representation
    vector<vector<float>> V; // vertice representation
    vector<float> innerPoint;

public:
    region();
    ~region();

    void WriteToDisk(ofstream& Outfile, bool OutHS);
    void ReadFromDisk(ifstream& Infile);
    void ToBeRoot(int dim);
    void ComputeHP(vector<float>& w, vector<float>& o1, vector<float>& o2, int dim);
    void FreeMem();
};


#endif //K_LEVEL_REGION_H
