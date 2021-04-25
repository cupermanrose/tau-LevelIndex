#include <iostream>

#include "k_level_lib.h"
#include "kspr.h"
#include "utk.h"

#define building 1

enum query_type{kspr, utk};

int main(int argc, char* argv[]) {
    int dim=6;
    int tau=30;
    char* datafile="/home/kemingli/klevel/data/real/NBA8D.dat";// TODO 改成相对路径
    string logfile="/home/kemingli/klevel/result/idx/NBA8D"+to_string(dim)+"_tau"+to_string(tau)+".log";// TODO 改成相对路径
    string idxfile="/home/kemingli/klevel/result/idx/NBA8D"+to_string(dim)+"_"+to_string(tau)+".idx";// TODO 改成相对路径
    cout<<"dim"<<dim<<endl;
    cout<<"tau"<<tau<<endl;
    cout<<string(datafile)<<endl;
    return 0;
    fstream log(logfile, ios::out);
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


