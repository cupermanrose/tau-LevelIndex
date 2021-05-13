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
#include "oru.h"

enum func_type{buildidx, loadidx};
enum query_type{kspr, utk, oru};
enum build_type{IncBuild,DFSBuild,BFSBuild};

void Config(int dim, int tau, int ik, string root_directory, string filename,
            string& datafile, string& logfile, string& idxfile, fstream& log,
            string func_str, func_type& func, string query_str, query_type& query, build_type& build, string build_str){

    if (func_str=="buildidx") func=buildidx;
    else if (func_str=="loadidx") func=loadidx;
    else cout << "Unknown function!" << endl;

    if (query_str=="kspr") query=kspr;
    else if (query_str=="utk") query=utk;
    else if (query_str=="oru") query=oru;
    else cout << "Unknown query!" << endl;

    if (build_str=="IncBuild") build=IncBuild;
    else if (build_str=="DFSBuild") build=DFSBuild;
    else if (build_str=="BFSBuild") build=BFSBuild;
    else cout << "Unknown building method!" << endl;

    datafile=root_directory+"data/"+filename+".dat";
    if (func==buildidx){
        logfile=root_directory+"log/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+".log";
    }
    else if (func==loadidx) {
        switch (query){
            case kspr:
                logfile=root_directory+"log/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_query"+query_str+".ans";
                break;
            case utk:
                logfile=root_directory+"log/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_query"+query_str+".ans";
                break;
            case oru:
                logfile=root_directory+"log/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_query"+query_str+".ans";
                break;
        }
    }
    idxfile=root_directory+"index/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_ik"+to_string(ik)+".idx";
    log.open(logfile, ios::out);
}

void ParameterInput(int argc, char* argv[], int& dim, int& tau, int& ik,
                    string& root_directory, string& filename, string& func_str, string& build_str, int& q_num, int& k, string& query_str){
    dim=4;
    tau=10; // NBA: tau=30
    ik=10;
    root_directory="/home/jiahaozhang/data/klevel/";
    filename="inde/U100K4";
    //filename="real/HOUSE6D";
    func_str="buildidx";
    build_str="IncBuild";

    q_num=5;
    k=1;
    query_str="oru";

}

int main(int argc, char* argv[]) {
    int dim,tau,ik,q_num,k;
    float utk_side_length;
    string datafile, logfile, idxfile, root_directory, filename, func_str,query_str, build_str;
    fstream log;
    query_type query;
    func_type func;
    build_type build;
    ParameterInput(argc, argv, dim,tau,ik,root_directory,filename,func_str, build_str, q_num, k, query_str);
    Config(dim,tau,ik,root_directory,filename,datafile,logfile,idxfile,log, func_str, func, query_str, query, build, build_str);

    level idx(dim,tau,ik);

    switch (func) {
        case buildidx:
            switch(build){
                case IncBuild:
                    IncBuildIndex(idx,datafile,log,idxfile);
                    break;
                case DFSBuild:
                    DFSBuildIndex(idx,datafile,log,idxfile);
                    break;
                case BFSBuild:
                    BuildIndex(idx, datafile, log, idxfile);
                    break;
            }
            break;
        case loadidx:
            LoadIndex(idx, datafile, log, idxfile);
            switch (query) {
                case kspr:
                    kspr::multiple_query(idx, 1, q_num, log);
                    kspr::multiple_query(idx, 5, q_num, log);
                    kspr::multiple_query(idx, 10, q_num, log);
                    kspr::multiple_query(idx, 15, q_num, log);
                    kspr::multiple_query(idx, 20, q_num, log);
                    kspr::multiple_query(idx, 30, q_num, log);
                    //kspr::multiple_query(idx, 40, q_num, log);
                    //kspr::multiple_query(idx, 50, q_num, log);
                    break;
                case utk:
                    utk_side_length=0.001;
                    cout << utk_side_length << endl;
                    log << utk_side_length << endl;
                    utk::multiple_query(idx, 1, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 5, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 10, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 15, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 20, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 30, q_num, utk_side_length, log);
                    utk_side_length=0.005;
                    cout << utk_side_length << endl;
                    log << utk_side_length << endl;
                    utk::multiple_query(idx, 1, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 5, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 10, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 15, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 20, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 30, q_num, utk_side_length, log);
                    utk_side_length=0.01;
                    cout << utk_side_length << endl;
                    log << utk_side_length << endl;
                    utk::multiple_query(idx, 1, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 5, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 10, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 15, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 20, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 30, q_num, utk_side_length, log);
                    utk_side_length=0.05;
                    cout << utk_side_length << endl;
                    log << utk_side_length << endl;
                    utk::multiple_query(idx, 1, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 5, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 10, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 15, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 20, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 30, q_num, utk_side_length, log);
                    utk_side_length=0.1;
                    cout << utk_side_length << endl;
                    log << utk_side_length << endl;
                    utk::multiple_query(idx, 1, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 5, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 10, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 15, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 20, q_num, utk_side_length, log);
                    utk::multiple_query(idx, 30, q_num, utk_side_length, log);
                    break;
                case oru:
                    k=1;
                    cout << k << endl;
                    oru::multiple_query(idx,k, 10,q_num,log);
                    oru::multiple_query(idx,k, 30,q_num,log);
                    oru::multiple_query(idx,k, 50,q_num,log);
                    oru::multiple_query(idx,k, 70,q_num,log);
                    oru::multiple_query(idx,k, 90,q_num,log);
                    k=5;
                    cout << k << endl;
                    oru::multiple_query(idx,k, 10,q_num,log);
                    oru::multiple_query(idx,k, 30,q_num,log);
                    oru::multiple_query(idx,k, 50,q_num,log);
                    oru::multiple_query(idx,k, 70,q_num,log);
                    oru::multiple_query(idx,k, 90,q_num,log);
                    k=10;
                    cout << k << endl;
                    oru::multiple_query(idx,k, 10,q_num,log);
                    oru::multiple_query(idx,k, 30,q_num,log);
                    oru::multiple_query(idx,k, 50,q_num,log);
                    oru::multiple_query(idx,k, 70,q_num,log);
                    oru::multiple_query(idx,k, 90,q_num,log);
                    k=15;
                    cout << k << endl;
                    oru::multiple_query(idx,k, 10,q_num,log);
                    oru::multiple_query(idx,k, 30,q_num,log);
                    oru::multiple_query(idx,k, 50,q_num,log);
                    oru::multiple_query(idx,k, 70,q_num,log);
                    oru::multiple_query(idx,k, 90,q_num,log);
                    k=20;
                    cout << k << endl;
                    oru::multiple_query(idx,k, 10,q_num,log);
                    oru::multiple_query(idx,k, 30,q_num,log);
                    oru::multiple_query(idx,k, 50,q_num,log);
                    oru::multiple_query(idx,k, 70,q_num,log);
                    oru::multiple_query(idx,k, 90,q_num,log);
                    break;
            }
            break;
    }

    cout << "DONE" << endl;
    log.close();
    return 0;
}


