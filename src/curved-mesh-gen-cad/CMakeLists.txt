
add_executable(cmg curvedmeshgen3d.cpp)
target_link_libraries(cmg amesh3d arbf alinalg amatrix2) 

add_executable(cmg2d curvedmeshgen2d.cpp)
target_link_libraries(cmg2d amesh2dh arbf alinalg amatrix2)