This repository is developed by Haisen Zhao for his research projects. A set of fundemental geometric process functions are well organized based on some third-party library. 
Besides that, a set of common used functions are included. Currently, we release a dynamic link library file of "carpentry_geom.dll". An example of using this library is pack-geom-test-lib. 

# Dependency

You should install [CGAL](https://github.com/CGAL/cgal) before compiling this project. For windows, recommend to use vcpkg to install CGAL. 
Our project also depends on [glm](https://github.com/g-truc/glm.git) and [Clipper](http://www.angusj.com/delphi/clipper.php), but you don't need to install them explicitly.

# Usage

- Clone this repository:

        git clone https://github.com/helm-compiler/carpentry-geometry.git

- Compiling.

        mkdir build
        cd build
        cmake ..
        make

# License
All rights about the program are reserved by the authors of this project. The programs can only be used for research purpose. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program.
