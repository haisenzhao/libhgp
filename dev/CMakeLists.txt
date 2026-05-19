cmake_minimum_required(VERSION 3.16)

# CMP0091 必须在 project() 之前声明！
cmake_policy(SET CMP0091 NEW)

project(libhgp)
# 显式保证所有目标使用 /MD（Release）或 /MDd（Debug）
#与 x64-windows-static-md triplet 保持一致
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

set(CMAKE_CXX_STANDARD 17)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_compile_options(/utf-8)
endif()

option(USE_LOCAL_LIBS "Use local library" ON)

# 通过 vcpkg 找 CGAL
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package(CGAL REQUIRED COMPONENTS Core)

# 公共源码列表
set(LIBHGP_SOURCES
    libhgp.h
    cgal.h
    hull3D.h
    kdtree.h
    hgp2d.cpp
    hgp3d.cpp
    hgpmesh.cpp
    hgpio.cpp
    clipper/clipper.cpp
    kdtree.cpp
    hull3D.cpp
    cgal.cpp
    lgp_api.cpp
)

# 动态库（原来的 libhgp.dll）
add_library(libhgp SHARED ${LIBHGP_SOURCES})

target_include_directories(libhgp
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/local_libs
)

target_link_libraries(libhgp
    PRIVATE
        CGAL::CGAL
        CGAL::CGAL_Core
)

# 静态库版本，用于 hgp_py 静态链接
add_library(libhgp_static STATIC ${LIBHGP_SOURCES})

# ✅ 关键：为静态库启用位置独立代码（PIC）
set_target_properties(libhgp_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(libhgp_static
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/local_libs
)

target_link_libraries(libhgp_static
    PUBLIC
        CGAL::CGAL
        CGAL::CGAL_Core
)

# 通过 vcpkg 找 pybind11，并让 hgp_py 链接静态库
option(BUILD_PYTHON_BINDINGS "Build Python bindings with pybind11" ON)
if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(hgp_py hgp_py.cpp)
    
    # 链接静态库
    target_link_libraries(hgp_py PRIVATE libhgp_static)

    # ✅ 可选但推荐：显式设置符号可见性
    set_target_properties(hgp_py PROPERTIES
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
    )
endif()