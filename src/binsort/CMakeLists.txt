
add_library(adistbinsort adistbinsort.cpp)
target_link_libraries(adistbinsort amatrix)
target_include_directories(adistbinsort PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
