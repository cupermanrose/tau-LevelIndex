#include "RboxPoints.h"
#include "QhullError.h"
#include "QhullQh.h"
#include "QhullFacet.h"
#include "QhullFacetList.h"
#include "QhullFacetSet.h"
#include "QhullLinkedList.h"
#include "QhullPoint.h"
#include "QhullUser.h"
#include "QhullVertex.h"
#include "QhullVertexSet.h"
#include "Qhull.h"
#include <iostream>
#include <iomanip> // setw
#include <unordered_map>
#include <unordered_set>
//TODO merge to qhull_adapter
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::cout;
using std::endl;
using std::string;
using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetSet;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexSet;
using orgQhull::RboxPoints;
using orgQhull::Coordinates;

#include "k_level_lib.h"
#include "kspr.h"
#include "utk.h"

#define building 1

enum query_type{kspr, utk};

int main(int argc, char* argv[]) {
    int dim=4;
    int tau=40;
    int k=10;
    int q_num=20;
    float utk_side_length=0.01;
    char* datafile="data/U400K4.dat";// TODO 改成相对路径
    char* logfile="result/debug_400K_40.log";// TODO 改成相对路径
    fstream log(logfile, ios::out);
    string idxfile="result/U400K4_"+to_string(tau)+".idx";// TODO 改成相对路径
    query_type query=utk;
    level idx(dim,tau);
    if (building){
        BuildIndex(idx, datafile, log, idxfile);
    }
    else{
        LoadIndex(idx, datafile, log, idxfile);
        switch (query) {
            case kspr:
                kspr::multiple_query(idx, k, q_num, log);
                break;
            case utk:
                utk::multiple_query(idx, k, q_num, utk_side_length, log);
                break;
        }
    }

    cout << "DONE" << endl;
    return 0;
}


