# Develop PPGL

- Install vcpkg and use vcpkg to install CGAL.

- Clone this repository:

      git clone git@github.com:haisenzhao/personal-pack-geom-lib.git

- Use Cmake-gui to Compile：([Cmake](https://cmake.org))

  - Click "Browse Source..." and "Browse Build..." to select the project path.![image-select_path](D:\pgl\personal-pack-geom-lib\images\image-select_path.png)
  - Click the "Configue" button and select the option "Specify toolchain file for cross-compiling" and "Next".
  - Select the "vcpkg.cmake" file path in the vcpkg installation directory.(For example:"D:/vcpkg/scripts/buildsystems/vcpkg.cmake")
  - Click the "Generate" button and "Open Project".
  - Set the sub-projects as startup projects from top to bottom, and compile them one by one.(Ensure network connectivity)
  - If you encounter errors, clear the cache and try again.("File"-"delete cache")