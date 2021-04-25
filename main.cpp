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
    int dim=6;
    int tau=20;
    char* datafile="/home/kemingli/klevel/data/inde/U400K4.dat";// TODO 改成相对路径
    string logfile="/home/kemingli/klevel/result/idx/U400K"+to_string(dim)+"_tau"+to_string(tau)+".log";// TODO 改成相对路径
    string idxfile="/home/kemingli/klevel/result/idx/U400K"+to_string(dim)+"_"+to_string(tau)+".idx";// TODO 改成相对路径
    cout<<"dim"<<dim<<endl;
    cout<<"tau"<<tau<<endl;
    cout<<string(datafile)<<endl;
    fstream log(logfile, ios::out);
//    return 0;
//    string idxfile="/home/jiahaozhang/data/klevel/results/ANTI400K_d"+to_string(dim)+"_tau"+to_string(tau)+".idx";// TODO 改成相对路径

    level idx(dim,tau);
    if (building){
        BuildIndex(idx, datafile, log, idxfile);
    }
    else{
        LoadIndex(idx, datafile, log, idxfile);

        query_type query=kspr;
        int k=10;
        int q_num=20;

        switch (query) {
            case kspr:
                kspr::multiple_query(idx, k, q_num, log);
                break;
            case utk:
                float utk_side_length=0.01;
                utk::multiple_query(idx, k, q_num, utk_side_length, log);
                break;
        }
    }

    cout << "DONE" << endl;
    return 0;
}


