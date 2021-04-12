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

int main(int argc, char* argv[]) {
    int dim=4;
    int tau=10;
    char* datafile="/home/jiahaozhang/data/klevel/U100K4.dat";// TODO 改成相对路径
    char* logfile="/home/jiahaozhang/projects/klevel/debug.log";// TODO 改成相对路径
    k_level_lib::BuildIdx(dim, tau, datafile, logfile);
    return 0;
}


