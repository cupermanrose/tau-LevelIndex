  To build Qhull, static libraries, shared library, and C++ interface
  - git clone https://github.com/qhull/qhull.git
  - cd qhull/build
  - cmake --help  # List build generators
  - cmake -G "<generator>" ..   
  - cmake --build .
  - sudo cmake --build . --target install
  
  
  To install gdb
  - sudo apt-get update
  - sudo apt-get install gdb
  
  To install boost
  - sudo apt-get update
  - sudo apt-get install libboost-all-dev
  
  To install osqp
  - git clone --recursive https://github.com/oxfordcontrol/osqp
  - cd osqp
  - mkdir build
  - cd build
  - cmake -G "Unix Makefiles" .. #
  - cmake --build .
  - sudo cmake --build . --target install

  
  
  log file
  
  U400K4_dim4_tau20
  
  datafile+dim+tau+"".log"
