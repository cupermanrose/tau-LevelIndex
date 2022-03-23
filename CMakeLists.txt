cmake_minimum_required(VERSION 3.10)
project(k_level)

set(CMAKE_CXX_STANDARD 14)
### add qhull
include_directories(/usr/local/include/libqhullcpp)
### qhullstatic reentrant library
add_library(libqhullstatic_r STATIC IMPORTED)
set_property(TARGET libqhullstatic_r PROPERTY IMPORTED_LOCATION /usr/local/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qhullstatic_r${CMAKE_STATIC_LIBRARY_SUFFIX})
add_dependencies(libqhullstatic_r qhull)
### qhullcpp library
add_library(libqhullcpp STATIC IMPORTED)
set_property(TARGET libqhullcpp PROPERTY IMPORTED_LOCATION /usr/local/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qhullcpp${CMAKE_STATIC_LIBRARY_SUFFIX})
set_property(TARGET libqhullcpp PROPERTY INTERFACE_LINK_LIBRARIES libqhullstatic_r)
add_dependencies(libqhullcpp qhull)

### add osqp
find_package(osqp REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed -O3 -Wall -ldl")



include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/bfp)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/bfp/bfp_LUSOL)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/bfp/bfp_LUSOL/LUSOL)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/colamd)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/lp_solve)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/lpsolve55)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib/shared)
include_directories(${CMAKE_SOURCE_DIR}/lib/lp_lib)
include_directories(${CMAKE_SOURCE_DIR}/klevel)
include_directories(${CMAKE_SOURCE_DIR}/klevel/algorithm)
include_directories(${CMAKE_SOURCE_DIR}/klevel/adapter)
include_directories(${CMAKE_SOURCE_DIR}/klevel/operator)
include_directories(${CMAKE_SOURCE_DIR}/lib/rtree)
include_directories(${CMAKE_SOURCE_DIR}/lib/qhull_lib)




add_executable(${PROJECT_NAME}
        main.cpp
        klevel/k_level_lib.h
        klevel/k_level_lib.cpp

        klevel/algorithm/region.h
        klevel/algorithm/region.cpp
        klevel/algorithm/kcell.h
        klevel/algorithm/kcell.cpp
        klevel/algorithm/level.h
        klevel/algorithm/level.cpp
        klevel/algorithm/utils.h
        klevel/algorithm/utils.cpp
        #[[klevel/algorithm/dominateG.h
        klevel/algorithm/dominateG.cpp]]

        klevel/adapter/lp_adapter.h
        klevel/adapter/lp_adapter.cpp
        klevel/adapter/qhull_adapter.h
        klevel/adapter/qhull_adapter.cpp
        klevel/adapter/rtree_adapter.h
        klevel/adapter/rtree_adapter.cpp
        klevel/adapter/qp_adapter.h
        klevel/adapter/qp_adapter.cpp

        klevel/operator/kspr.h
        klevel/operator/kspr.cpp
        klevel/operator/utk.h
        klevel/operator/utk.cpp
        klevel/operator/oru.h
        klevel/operator/oru.cpp
        klevel/operator/topk.h
        klevel/operator/topk.cpp


        lib/lp_lib/bfp/bfp_LUSOL/lp_LUSOL.h
        lib/lp_lib/bfp/bfp_LUSOL/LUSOL/lusol.h
        lib/lp_lib/bfp/lp_BFP.h
        lib/lp_lib/colamd/colamd.h
        lib/lp_lib/ini.h
        lib/lp_lib/shared/commonlib.h
        lib/lp_lib/shared/mmio.h
        lib/lp_lib/shared/myblas.h
        lib/lp_lib/lp_crash.h
        lib/lp_lib/lp_Hash.h
        lib/lp_lib/lp_lib.h
        lib/lp_lib/lp_matrix.h
        lib/lp_lib/lp_MDO.h
        lib/lp_lib/lp_mipbb.h
        lib/lp_lib/lp_MPS.h
        lib/lp_lib/lp_presolve.h
        lib/lp_lib/lp_price.h
        lib/lp_lib/lp_pricePSE.h
        lib/lp_lib/lp_report.h
        lib/lp_lib/lp_scale.h
        lib/lp_lib/lp_simplex.h
        lib/lp_lib/lp_SOS.h
        lib/lp_lib/lp_utils.h
        lib/lp_lib/lp_wlp.h
        lib/lp_lib/bfp/bfp_LUSOL/lp_LUSOL.c
        lib/lp_lib/bfp/bfp_LUSOL/LUSOL/lusol.c
        lib/lp_lib/colamd/colamd.c
        lib/lp_lib/ini.c
        lib/lp_lib/shared/commonlib.c
        lib/lp_lib/shared/mmio.c
        lib/lp_lib/shared/myblas.c
        lib/lp_lib/lp_crash.c
        lib/lp_lib/lp_Hash.c
        lib/lp_lib/lp_lib.c
        lib/lp_lib/lp_matrix.c
        lib/lp_lib/lp_MDO.c
        lib/lp_lib/lp_mipbb.c
        lib/lp_lib/lp_MPS.c
        lib/lp_lib/lp_params.c
        lib/lp_lib/lp_presolve.c
        lib/lp_lib/lp_price.c
        lib/lp_lib/lp_pricePSE.c
        lib/lp_lib/lp_report.c
        lib/lp_lib/lp_scale.c
        lib/lp_lib/lp_simplex.c
        lib/lp_lib/lp_SOS.c
        lib/lp_lib/lp_utils.c
        lib/lp_lib/lp_wlp.c


        lib/quadProg_lib/Array.hh
        lib/quadProg_lib/Array.cc
        lib/quadProg_lib/QuadProg++.hh
        lib/quadProg_lib/QuadProg++.cc

        lib/rtree/collection.h
        lib/rtree/collection.cpp
        lib/rtree/filemem.h
        lib/rtree/filemem.cpp
        lib/rtree/global.h
        lib/rtree/global.cpp
        lib/rtree/header.h
        lib/rtree/hypercube.h
        lib/rtree/hypercube.cpp
        lib/rtree/mem.h
        lib/rtree/param.h
        lib/rtree/param.cpp
        lib/rtree/point.h
        lib/rtree/point.cpp
        lib/rtree/rentry.h
        lib/rtree/rentry.cpp
        lib/rtree/rnode.h
        lib/rtree/rnode.cpp
        lib/rtree/rtree.h
        lib/rtree/rtree.cpp
        lib/rtree/tgs.h
        lib/rtree/tgs.cpp
        lib/rtree/virtualRNode.h
        lib/rtree/virtualRNode.cpp
        )
find_package(Boost 1.65.1)
if(Boost_FOUND)
    include_directories(BEFORE ${Boost_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
else()
    message(STATUS "This program requires the boost library, and will not be compiled.")
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE osqp::osqpstatic)

target_link_libraries(${PROJECT_NAME} PRIVATE libqhullcpp)

