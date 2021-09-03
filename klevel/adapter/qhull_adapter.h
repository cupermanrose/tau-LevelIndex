
#ifndef K_LEVEL_QHULLADAPTER_H
#define K_LEVEL_QHULLADAPTER_H

#include "lp_adapter.h"
//#include "libqhull_r/libqhull_r.h"
//#include "libqhull_r/qhull_ra.h"
#include "Qhull.h"
#include "QhullUser.h"
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
#include <unordered_map>
#include <unordered_set>

class qhull_adapter{
public:
    qhull_adapter();
    ~qhull_adapter();
    static void ComputeVertex(vector<halfspace>& H, vector<vector<float>>& V, vector<float>& innerPoint, int& dim);
    static void ComputeVertex2D(vector<halfspace>& H, vector<vector<float>>& V, vector<float>& innerPoint);
};


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

class qhull_user{

#define POINT_ID int
#define REGION std::vector<std::vector<double>>
public:
    qhull_user();

    static void points_at_half_inter(std::vector<std::vector<float>> &ret, const realT *pointCoordinates, int p_num, vector<float> &innerPoint);

    void get_neiVT_of_VT(Qhull &q, const std::vector<int>&pd_ids, std::unordered_map<int, std::vector<int>> &ret);

    std::unordered_map<int, std::vector<int>> get_neiVT_of_VT(Qhull &q, const std::vector<int>&pd_ids);

    std::vector<std::vector<int>> get_points_of_facets(Qhull &q, const std::vector<int> &pd_ids);

    std::vector<int> get_CH_pointID(Qhull &q, const std::vector<int> &pd_ids);

    std::unordered_map<POINT_ID, REGION> get_neiFacetsNorm_of_point(Qhull &q, const std::vector<int> &pd_ids);

    void get_neiFacetsNorm_of_point(Qhull &q, const std::vector<int> &pd_ids, std::unordered_map<POINT_ID, REGION> &ret);

    std::vector<std::vector<double>> get_norm_of_all_facets(Qhull &q);

    std::vector<std::vector<double>> get_cone_norms(Qhull &q, std::vector<std::vector<double>> &points);

    std::vector<std::vector<int>> get_neiFacets_of_points(Qhull &q, const std::vector<int> &pd_ids);
};


#endif //K_LEVEL_QHULLADAPTER_H