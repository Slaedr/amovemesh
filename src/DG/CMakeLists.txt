add_executable(triangulate2d triangulate2d.cpp)
target_link_libraries(triangulate2d adistbinsort amesh2dh)
target_include_directories(triangulate2d PUBLIC ${CMAKE_SOURCE_DIR}/binsort/)
