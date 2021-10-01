project(pack-geom-test-lib)

# Create the executable
add_executable(${PROJECT_NAME} main.cpp cgal.h)

target_link_libraries(${PROJECT_NAME} sub::pglib)

include(ExternalProject)

ExternalProject_Get_Property(glm source_dir)
set(GlmIncludeDir ${source_dir})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET hidden)
add_dependencies(${PROJECT_NAME} glm)