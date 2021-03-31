// 1. using config file
//void fun1();
// 2. directly pass parameter and using API
//void fun2();

#include "k_level_lib.h"

int main(int argc, char* argv[]) {
    int dim,tau;
    string datafile;
    k_level_lib::BuildIdx(dim,tau,datafile);
    return 0;
}
