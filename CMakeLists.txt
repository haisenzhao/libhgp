cmake_minimum_required (VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)

project(port)

#add _CRT_SECURE_NO_WARNINGS
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()


# Create the executable
add_executable(${PROJECT_NAME} port.cpp)


include(ExternalProject)

include(ExternalProject)
ExternalProject_Add(
    pgl
    PREFIX ${CMAKE_BINARY_DIR}/third_party/pgl
    GIT_REPOSITORY https://github.com/haisenzhao/personal-geom-lib.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(pgl source_dir)
set(PglIncludeDir ${source_dir})

include_directories(${PglIncludeDir})
add_dependencies(${PROJECT_NAME} pgl)