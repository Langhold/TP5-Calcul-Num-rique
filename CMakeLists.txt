cmake_minimum_required(VERSION 3.15)
project(ppm LANGUAGES C)

add_executable(test src/testenv.c)

add_library(lib include/atlas_headers.h)

target_link_libraries(test PRIVATE lib)
