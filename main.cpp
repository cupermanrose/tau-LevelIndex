// 1. using config file
//void fun1();
// 2. directly pass parameter and using API
//void fun2();

#include "k_level_lib.h"

int main(int argc, char* argv[]) {
    int dim=4;
    int tau=10;
    char* datafile="/home/jiahaozhang/data/klevel/U100K4.dat";// TODO 改成相对路径
    char* logfile="/home/jiahaozhang/projects/klevel/result/debug.log";// TODO 改成相对路径
    k_level_lib::BuildIdx(dim, tau, datafile, logfile);
    return 0;
}
