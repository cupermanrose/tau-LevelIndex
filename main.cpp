#include "k_level_lib.h"
#include "kspr.h"
#include "utk.h"
#include "oru.h"
#include "topk.h"
#define BUILDIDX "buildidx"
#define LOADIDX "loadidx"
#define KSPR "kspr"
#define UTK "utk"
#define ORU "oru"
#define TOPK "topk"

extern bool apply_onion_from_file;
extern bool write_onion_to_file;
extern string anti_id_f;
extern string read_anti_dat;


enum func_type{buildidx, loadidx};
enum query_type{kspr, utk, oru, topk};
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
    else if (query_str==TOPK) query=topk;
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
    dim=4; tau=10; ik=10;
//    root_directory="~/data/klevel/";
    root_directory="~/klevel/";
    filename="real/HOTEL4D";
//    func_str="buildidx"; // buildidx loadidx
    func_str="loadidx";
    build_str="PTF"; // INS PT PTF
    query_str="topk"; // kspr utk oru topk
    q_num=30; // # of query
    k=10; // query k

    // auxiliary parameter
    anti_id_f=root_directory+"data/"+filename+".ch";
    read_anti_dat=root_directory+"data/"+filename;
    apply_onion_from_file= false;
    write_onion_to_file=false;
}

std::string strip(std::string in)
{
    in.erase(std::remove_if(in.begin(), in.end(), [] (std::string::value_type ch)
                            { return isspace(ch); }
    ), in.end());
    return in;
}

void read_configFile(int argc, char* argv[], int& dim, int& tau, int& ik,
                     string& root_directory, string& filename, string& func_str, string& build_str, int& q_num, int& k, string& query_str){
    // default setting
    dim=4; tau=10; ik=10;
//    root_directory="~/data/klevel/";
    root_directory="~/klevel/";
    filename="inde/U400K4";
//    func_str="buildidx"; // buildidx loadidx
    func_str="loadidx";
    build_str="PTF"; // INS PT PTF
    query_str="topk"; // kspr utk oru topk
    q_num=30; // # of query
    k=1; // query k

    const char *config=argv[1];
    char str[1000]; // only support line len within 1000 char
    ifstream fconfig(config);
    while(fconfig.getline(str, 1000)){
        string s(str);
        auto pos_c=s.find("#");
        if(pos_c!=s.npos){
            s=s.substr(0, pos_c);
        }

        s=strip(s);

        auto pos_eq=s.find("=");
        if(pos_eq==s.npos){
            continue;
        }
        string name=s.substr(0, pos_eq);
        string property=s.substr(pos_eq+1, s.npos);

        if(strcmp(name.c_str(), "dim")==0){
            dim=stoi(property);
        }else if(strcmp(name.c_str(), "tau")==0){
            tau=stoi(property);
        }else if(strcmp(name.c_str(), "ik")==0){
            ik=stoi(property);
        }else if(strcmp(name.c_str(), "root_directory")==0){
            root_directory=property;
            if(root_directory[0]=='~'){
                const unsigned long maxLen = 64;
                char currentUser[maxLen];
                getlogin_r(currentUser, maxLen);
                root_directory= string("/home/") + string(currentUser)+string("/")+root_directory.substr(1, root_directory.npos);
            }
        }else if(strcmp(name.c_str(), "filename")==0){
            filename=property;
        }else if(strcmp(name.c_str(), "func_str")==0){
            func_str=property;
        }else if(strcmp(name.c_str(), "query_str")==0){
            query_str=property;
        }else if(strcmp(name.c_str(), "q_num")==0){
            q_num=stoi(property);
        }else if(strcmp(name.c_str(), "k")==0){
            k=stoi(property);
        }else if(strcmp(name.c_str(), "build_str")==0) {
            build_str=property;
        }else{
            cout<<"unrecognized property name "<< name <<endl;
            exit(-1);
        }
    }
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
    read_configFile(argc, argv, dim,tau,ik,root_directory,filename,func_str, build_str, q_num, k, query_str);
//    ParameterInput(argc, argv, dim,tau,ik,root_directory,filename,func_str, build_str, q_num, k, query_str);
    Config(dim,tau,ik,root_directory,filename,datafile,logfile,idxfile,log, func_str, func, query_str, query, build, build_str);
    cout<<datafile<<endl;
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
                case topk: {
                    dtopk::multiple_query(idx, k, q_num, log);
                }
            }
            break;
    }
    cout << "DONE" << endl;
    log.close();
    return 0;
}


