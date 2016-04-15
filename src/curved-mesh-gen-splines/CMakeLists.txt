
add_executable(amc curve3d.cpp)
target_link_libraries(amc amesh3d arbf alinalg ageometry3d amatrix2) 

add_executable(curveh curvedmeshgen2dh.cpp)
target_link_libraries(curveh arbf ageometryh amesh2dh adatastructures amatrix2)

add_executable(curveelast linelast-curvedmeshgen2d.cpp)
target_link_libraries(curveelast arbf ageometry adatastructures amatrix2)