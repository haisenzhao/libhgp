cmake_minimum_required(VERSION 3.5)
project(libhgp)

#add _CRT_SECURE_NO_WARNINGS
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()


#define option to use local libs
option(USE_LOCAL_LIBS "Use local library" ON)


set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package( CGAL QUIET COMPONENTS Core)
if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()
include( ${CGAL_USE_FILE} )

if(USE_LOCAL_LIBS)
   set(PglIncludeDir ${PROJECT_SOURCE_DIR}/../local_libs)
   set(GlmIncludeDir ${PROJECT_SOURCE_DIR}/../local_libs)
   set(EigenIncludeDir ${PROJECT_SOURCE_DIR}/../local_libs)
   set(GteIncludeDir ${PROJECT_SOURCE_DIR}/../local_libs)
else()
	include(ExternalProject)
	ExternalProject_Add(
		pgl
		PREFIX ${CMAKE_BINARY_DIR}/third_party/pgl
		GIT_REPOSITORY https://github.com/haisenzhao/personal-geom-lib.git
		CONFIGURE_COMMAND ""
		UPDATE_DISCONNECTED 1
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		LOG_DOWNLOAD ON
		)
	ExternalProject_Get_Property(pgl source_dir)
	set(PglIncludeDir ${source_dir})
	
	ExternalProject_Add(
		glm
		PREFIX ${CMAKE_BINARY_DIR}/third_party/glm
		GIT_REPOSITORY https://github.com/g-truc/glm.git
		CONFIGURE_COMMAND ""
		UPDATE_DISCONNECTED 1
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		LOG_DOWNLOAD ON
		)
	ExternalProject_Get_Property(glm source_dir)
	set(GlmIncludeDir ${source_dir})
	
	ExternalProject_Add(
		eigen
		PREFIX ${CMAKE_BINARY_DIR}/third_party/eigen
		GIT_REPOSITORY https://github.com/libigl/eigen.git
		CONFIGURE_COMMAND ""
		UPDATE_DISCONNECTED 1
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		LOG_DOWNLOAD ON
		)
	ExternalProject_Get_Property(eigen source_dir)
	set(EigenIncludeDir ${source_dir})
	
	
	ExternalProject_Add(
		gte
		PREFIX ${CMAKE_BINARY_DIR}/third_party/gte
		GIT_REPOSITORY https://github.com/davideberly/GeometricTools.git
		CONFIGURE_COMMAND ""
		UPDATE_DISCONNECTED 1
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		LOG_DOWNLOAD ON
		)
	ExternalProject_Get_Property(gte source_dir)
	set(GteIncludeDir ${source_dir}/GTE)
endif()



add_library(${PROJECT_NAME} SHARED geom.h twoD.cpp threeD.cpp mesh.cpp clipper/clipper.cpp io.cpp  include_cgal.h  kdtree.h kdtree.cpp NewtonApple_hull3D.h NewtonApple_hull3D.cpp)
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Export interface symbols and hide internal implementation
include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME})
add_library(sub::pglib ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir} ${PglIncludeDir} ${EigenIncludeDir} ${GteIncludeDir})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET hidden)


if(NOT USE_LOCAL_LIBS)
	add_dependencies(${PROJECT_NAME} glm)
	add_dependencies(${PROJECT_NAME} pgl)
	add_dependencies(${PROJECT_NAME} eigen)
	add_dependencies(${PROJECT_NAME} gte)
endif()



#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND start ${CMAKE_BINARY_DIR}/post/Release/post.exe
#)


