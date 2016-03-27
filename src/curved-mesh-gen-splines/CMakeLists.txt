
add_executable(amc curve3d.cpp)
target_link_libraries(amc amesh3d arbf alinalg ageometry3d amatrix2) 

add_executable(curveh.out curvedmeshgen2dh.cpp)
target_link_libraries(curveh.out arbf ageometryh amesh2dh adatastructures amatrix2)