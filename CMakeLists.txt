cmake_minimum_required(VERSION 4.0)

project(ShroedingerVisualizer)

add_executable(main src/main.cpp)
target_compile_features(main PUBLIC cxx_std_17)