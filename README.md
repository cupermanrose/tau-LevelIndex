#   ![equation](http://www.sciweavers.org/tex2img.php?eq=%20%5Ctau%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)-levelIndex

This is a cmake project, which also means it is a C plus plus project for paper:

![equation](http://www.sciweavers.org/tex2img.php?eq=%20%5Ctau%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)-levelIndex: Towards Efficient Query Processing in Continuous Preference Space

This project require serveral libraries, e.g., 
- Qhull, a lib that computes the convex hull, Delaunay triangulation, Voronoi diagram, halfspace intersection about a point, furthest-site Delaunay triangulation, and furthest-site Voronoi diagram. Its official website is http://www.qhull.org/  
- lpsolve, a lib that can solve (mix integer) linear programming problem. Its official website is http://lpsolve.sourceforge.net/5.5/index.htm
- rtree, a data structure lib used for spatial data 
- boost, a general C++ lib. In this project we mainly use its hash lib. Its official website is https://www.boost.org/
- osqp, a Quadratic Programming lib and its offcial website is https://osqp.org/

In the next section, we introduce how to install the lib mentioned above one by one.
(to those who is familiar with cmake you can skip step 3 and 4, 
and load this project with a cmake IDE, e.g., CLion, https://www.jetbrains.com/clion/; of course, 
make sure your path is correct respecting to step 4)

Step 1: download make and cmake <br />
-------------

(usually a linux should contain make, you could try with command "make -version")
if your machine doesn't contain make, you can follow the steps shown as below:

```
sudo apt-get install build-essential
```

try "cmake -version" to see if your machine contains cmake, if doesn't, install cmake as: 
``` 
sudo snap install cmake --classic 
```


Step 2 install gdb, Qhull,  boost, and osqp<br />
------------- 
(lpsolve and rtree have been embedded in our folder /lib/lp_lib and /lib/rtree)

  To install gdb
  ```
  sudo apt-get update
  sudo apt-get install gdb
  ```
  
  To build Qhull, static libraries, shared library, and C++ interface
  ```
  git clone https://github.com/qhull/qhull.git
  cd qhull/build
  cmake --help  # List build generators, get generator for the next step, usually in linux it shall be "Unix Makefiles"
  cmake -G "<generator>" ..    # e.g.,  cmake -G "Unix Makefiles" ..   
  cmake --build .
  sudo cmake --build . --target install
  ```
  
  To install boost
  ```
  sudo apt-get update
  sudo apt-get install libboost-all-dev
  ```  

  To install osqp
  ```
  git clone --recursive https://github.com/oxfordcontrol/osqp
  cd osqp
  mkdir build
  cd build
  cmake -G "Unix Makefiles" .. #
  cmake --build .
  sudo cmake --build . --target install
  ```

  Step 3: compile this project
------------
  ```
  cd ~
  git clone https://github.com/a2xldmVs/klevel.git 
  cd klevel 
  mkdir build
  cd build
  cmake --help # List build generators, the system will tell you which generator you should use
  cmake -G "" .. # e.g., "cmake -G "Unix Makefiles" .."
  cmake --build . # it may should a lot of warnings because I set "-Wall" in CMakelists.txt
  ``` 
  
 Step 4: run this project
 ------------
 We have set some config files in config/ <br />
 You can build an index as setting in config_build1.txt and run as:
 ```
./k_level ../config/config_build1.txt
```
After building an index file, you could load that index file and 
answer some queries such as utk as setting in config_query1.txt and run as:
```
./k_level ../config/config_query1.txt
```
See more explanation details in config_build1.txt or config_query1.txt

Dataset
--------
The full dataset can be accessed by https://drive.google.com/drive/folders/1TaUDTRGgEVmkSvdTbDbYQb5LH0r6gNdm?usp=sharing
