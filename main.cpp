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
enum build_type{INS,PT, PTF};

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

    if (build_str=="INS") build=INS;
    else if (build_str=="PT") build=PT;
    else if (build_str=="PTF") build=PTF;
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
    idxfile=root_directory+"index/"+filename+"_dim"+to_string(dim)+"_tau"+to_string(tau)+"_ik"+to_string(ik)+"_"+build_str+".idx";
    cout<<idxfile<<endl;
    log.open(logfile, ios::out);
}

void ParameterInput(int argc, char* argv[], int& dim, int& tau, int& ik,
                    string& root_directory, string& filename, string& func_str, string& build_str, int& q_num, int& k, string& query_str){
    dim=5; tau=10; ik=10;
//    root_directory="~/data/klevel/";
    root_directory="~/klevel/";
    filename="inde/U400K5";
    func_str="buildidx"; // buildidx loadidx
//    func_str="loadidx";
    build_str="PT"; // INS PT PTF
    query_str="oru"; // kspr utk oru
    q_num=20; // # of query
    k=10; // query k

    // auxiliary parameter
    anti_id_f=root_directory+"data/"+filename+".ch";
    read_anti_dat=root_directory+"data/"+filename;
    apply_onion_from_file= false;
    write_onion_to_file=false;
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
                case INS:
                    INSBuild(idx, datafile, log, idxfile);
                    break;
                case PT:
                    PTBuild(idx, datafile, log, idxfile);
                    break;
                case PTF:
                    PTFBuild(idx, datafile, log, idxfile);
                    break;
            }
            break;
        case loadidx:
            LoadIndex(idx, datafile, log, idxfile);
            switch (query) {
                case kspr:{
                    kspr::multiple_query(idx, k, q_num, log);
                    break;
                }
                case utk: {
                    utk::multiple_query(idx, k, q_num, 0.01, log);
                    break;
                }
                case oru: {
                    oru::multiple_query(idx, k, 50, q_num, log);
                    break;
                }
            }
            break;
    }
    cout << "DONE" << endl;
    log.close();
    return 0;
}


