project(post)

# Create the executable
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} sub::pglib)

include(ExternalProject)

ExternalProject_Get_Property(glm source_dir)
set(GlmIncludeDir ${source_dir})

ExternalProject_Get_Property(pgl source_dir)
set(PglIncludeDir ${source_dir})

ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir} ${PglIncludeDir} ${EigenIncludeDir})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET hidden)
add_dependencies(${PROJECT_NAME} glm)
add_dependencies(${PROJECT_NAME} pgl)
add_dependencies(${PROJECT_NAME} eigen)