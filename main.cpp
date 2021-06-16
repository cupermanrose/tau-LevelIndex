

#include "k_level_lib.h"
#include "kspr.h"
#include "utk.h"
#include "oru.h"
#define BUILDIDX "buildidx"
#define LOADIDX "loadidx"
#define KSPR "kspr"
#define UTK "utk"
#define ORU "oru"


extern bool apply_onion_from_file;
extern bool write_onion_to_file;
extern string anti_id_f;
extern string read_anti_dat;


enum func_type{buildidx, loadidx};
enum query_type{kspr, utk, oru};
enum build_type{IncBuild,DFSBuild,BFSBuild};

void Config(int dim, int tau, int ik, string root_directory, string filename,
            string& datafile, string& logfile, string& idxfile, fstream& log,
            string func_str, func_type& func, string query_str, query_type& query, build_type& build, string build_str){

    if (func_str==BUILDIDX) func=buildidx;
    else if (func_str==LOADIDX) func=loadidx;
    else cout << "Unknown function!" << endl;

    if (query_str==KSPR) query=kspr;
    else if (query_str==UTK) query=utk;
    else if (query_str==ORU) query=oru;
    else cout << "Unknown query!" << endl;

    if (build_str=="IncBuild") build=IncBuild;
    else if (build_str=="DFSBuild") build=DFSBuild;
    else if (build_str=="BFSBuild") build=BFSBuild;
    else cout << "Unknown building method!" << endl;

    datafile=root_directory+"data/"+filename+".dat";
    if (func==buildidx){
        logfile=root_directory+"log/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_"+build_str+".log";
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
    idxfile=root_directory+"index/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_ik"+to_string(ik)+"_noHS_"+build_str+".idx";
    log.open(logfile, ios::out);
}

void ParameterInput(int argc, char* argv[], int& dim, int& tau, int& ik,
                    string& root_directory, string& filename, string& func_str, string& build_str, int& q_num, int& k, string& query_str){
    dim=4;
    tau=10; // NBA: tau=30
    ik=10;
    root_directory="/home/jiahaozhang/data/klevel/";
    filename="inde/U400K4";
    //func_str="buildidx";
    func_str="loadidx";
    build_str="BFSBuild";
    anti_id_f=root_directory+"data/"+filename+"tau100.ch";
    apply_onion_from_file=false;
    write_onion_to_file=false;
    q_num=10;
    k=10;
    query_str="utk";

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
                case kspr:{
                    kspr::multiple_query(idx, 10, q_num, log);
                    break;
                }
                case utk: {
                    clock_t rtree_time = clock();
                    Rtree *rt = nullptr;
                    unordered_map<long int, RtreeNode *> ramTree;
                    BuildRtree(idx.idx[idx.ik], rt, ramTree);
                    log << "R-tree from k-level building time: " << (clock() - rtree_time) / (float) CLOCKS_PER_SEC
                        << endl;
                    cout << "R-tree from k-level building time: " << (clock() - rtree_time) / (float) CLOCKS_PER_SEC
                         << endl;
                    utk::multiple_query(idx, rt, ramTree, 10, q_num, 0.01, log);
                    break;
                }
                case oru: {
                    clock_t rtree_time = clock();
                    Rtree *rt = nullptr;
                    unordered_map<long int, RtreeNode *> ramTree;
                    BuildRtree(idx.idx[idx.ik], rt, ramTree);
                    log << "R-tree from k-level building time: " << (clock() - rtree_time) / (float) CLOCKS_PER_SEC
                        << endl;
                    cout << "R-tree from k-level building time: " << (clock() - rtree_time) / (float) CLOCKS_PER_SEC
                         << endl;
                    oru::multiple_query(idx, rt, ramTree, 10, 50, q_num, log);
                    break;
                }
            }
            break;
    }
    cout << "DONE" << endl;
    log.close();
    return 0;
}


