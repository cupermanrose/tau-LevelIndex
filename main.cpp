// 1. using config file
//void fun1();
// 2. directly pass parameter and using API
//void fun2();

#include "k_level_lib.h"

int main(int argc, char* argv[]) {
    int dim=4;
    int tau=10;
    char* datafile="/home/jiahaozhang/data/klevel/U100K4.dat";
    char* logfile="/home/jiahaozhang/data/klevel/U100K4.log";
    k_level_lib::BuildIdx(dim, tau, datafile, logfile);
    return 0;
}
