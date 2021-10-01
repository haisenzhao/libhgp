# Geometric Library for Carpentry Compiler [Project page](https://grail.cs.washington.edu/projects/carpentrycompiler/)

Our carpentry compiler converts high-level geometric designs made by users to low-level fabrication instructions that can be directly followed to manufacture parts. The compiler performs multi-objective optimization on the low-level instructions to generate Pareto-optimal candidates. This repository provides a set of fundemental geometric process functions for Carpentry Compiler, with the output of a dynamic link library file, "carpentry_geom.dll".

# Dependency

You should install [CGAL](https://github.com/CGAL/cgal) before compiling this project.
Our project also depends on [glm](https://github.com/g-truc/glm.git) and [Clipper](http://www.angusj.com/delphi/clipper.php), but you don't need to install them explicitly.

# Usage

- Clone this repository:

        git clone https://github.com/helm-compiler/carpentry-geometry.git

- Compiling.

        mkdir build
        cd build
        cmake ..
        make

# Citation
If you use this code for your research, please cite our [paper](hhttps://grail.cs.washington.edu/projects/carpentrycompiler/files/CarpentryCompiler.pdf):

```
@article {wu_siga19,
    author = {Chenming Wu and Haisen Zhao and Chandrakana Nandi and Jeffrey I. Lipton and Zachary Tatlock and Adriana Schulz},
    title = {Carpentry Compiler},
    journal = {ACM Transactions on Graphics},
    note = {presented at SIGGRAPH Asia 2019},
    volume = {38},
    number = {6},
    pages = {Article No. 195},
    year = {2019}
}
```

# License
All rights about the program are reserved by the authors of this project. The programs can only be used for research purpose. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program.
