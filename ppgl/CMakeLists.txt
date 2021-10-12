cmake_minimum_required(VERSION 3.5)
project(ppgl)

#add _CRT_SECURE_NO_WARNINGS
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package( CGAL QUIET COMPONENTS Core)
if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()
include( ${CGAL_USE_FILE} )


include(ExternalProject)
ExternalProject_Add(
    pgl
    PREFIX ${CMAKE_BINARY_DIR}/pgl
    GIT_REPOSITORY https://github.com/haisenzhao/personal-geom-lib.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(pgl source_dir)
set(PglIncludeDir ${source_dir})

ExternalProject_Add(
    glm
    PREFIX ${CMAKE_BINARY_DIR}/glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(glm source_dir)
set(GlmIncludeDir ${source_dir})

ExternalProject_Add(
    eigen
    PREFIX ${CMAKE_BINARY_DIR}/eigen
    GIT_REPOSITORY https://github.com/libigl/eigen.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})


ExternalProject_Add(
    gte
    PREFIX ${CMAKE_BINARY_DIR}/gte
    GIT_REPOSITORY https://github.com/davideberly/GeometricTools.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(gte source_dir)
set(GteIncludeDir ${source_dir}/GTE)

add_library(${PROJECT_NAME} SHARED geom.h twoD.cpp threeD.cpp mesh.cpp clipper/clipper.cpp io.cpp cgal.h include_cgal.h)
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Export interface symbols and hide internal implementation
include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME})
add_library(sub::pglib ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir} ${PglIncludeDir} ${EigenIncludeDir} ${GteIncludeDir})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET hidden)
add_dependencies(${PROJECT_NAME} glm)
add_dependencies(${PROJECT_NAME} pgl)
add_dependencies(${PROJECT_NAME} eigen)
add_dependencies(${PROJECT_NAME} gte)

#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#        ${PROJECT_BINARY_DIR}/Release/mpir.dll
#        ${CMAKE_CURRENT_SOURCE_DIR}/dll/
#)

#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#        ${PROJECT_BINARY_DIR}/Release/ppgl.dll
#        ${CMAKE_CURRENT_SOURCE_DIR}/dll/
#)
#
#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#        ${PROJECT_BINARY_DIR}/Release/mpir.dll
#        ${PROJECT_BINARY_DIR}/../test/Release/
#)
#
#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#        ${PROJECT_BINARY_DIR}/Release/ppgl.dll
#        ${PROJECT_BINARY_DIR}/../test/Release/
#)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND start E:/Task2/personal-pack-geom-lib/build/post/Release/post.exe
)


#@echo xxxxxxxxxxx if $(PlatformTarget)==x64 (xcopy "$(ProjectDir)\Release\ppgl.dll" "$(ProjectDir)\..\..\pack-geom-lib\dll\" /y ) )