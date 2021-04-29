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
#define BUILDIDX "buildidx"
#define LOADIDX "loadidx"
#define KSPR "kspr"

enum func_type{buildidx, loadidx};
enum query_type{kspr, utk, oru};

void Config(int dim, int tau, string root_directory, string filename,
            string& datafile, string& logfile, string& idxfile, fstream& log,
            string func_str, func_type& func, string query_str, query_type& query){

    datafile=root_directory+"data/"+filename+".dat";
    logfile=root_directory+"log/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+".log";// TODO in case of user doesn't have this folder
    idxfile=root_directory+"index/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+".idx";
    log.open(logfile, ios::out);
    cout<<logfile<<endl;

    if (func_str==BUILDIDX) {
        func=buildidx;
    }else if (func_str==LOADIDX) {
        func=loadidx;
    }else {
        cout << "Unknown function!" << endl;
    }

    if (query_str=="kspr") query=kspr;
    else if (query_str=="utk") query=utk;
    else if (query_str=="oru") query=oru;
    else cout << "Unknown query!" << endl;
}

void ParameterInput(int argc, char* argv[], int& dim, int& tau,
                    string& root_directory, string& filename, string& func_str, int& q_num, int& k, string& query_str){
    dim=4;
    tau=1;
    root_directory="/home/kemingli/klevel/";
    filename="inde/U400K4";
    func_str=LOADIDX;
    q_num=50;
    k=1;
    query_str="utk";
}

int main(int argc, char* argv[]) {
    int dim,tau,q_num, k;
    string datafile, logfile, idxfile, root_directory, filename, func_str,query_str;
    fstream log;
    query_type query;
    func_type func;
    ParameterInput(argc, argv, dim,tau,root_directory,filename,func_str, q_num, k, query_str);
    Config(dim,tau,root_directory,filename,datafile,logfile,idxfile,log, func_str, func, query_str, query);

    level idx(dim,tau);

    switch (func) {
        case buildidx:
            BuildIndex(idx, datafile, log, idxfile);
            break;
        case loadidx:
            LoadIndex(idx, datafile, log, idxfile);
            switch (query) {
                case kspr:
                    kspr::multiple_query(idx, k, q_num, log);
                    break;
                case utk:
                    float utk_side_length=0.01;
                    utk::multiple_query(idx, k, q_num, utk_side_length, log);
                    break;
                /*case oru:
                    float ret_size=10;
                    break;*/
            }
    }

    cout << "DONE" << endl;
    log.close();
    return 0;
}


