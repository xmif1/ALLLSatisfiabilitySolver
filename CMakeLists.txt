cmake_minimum_required(VERSION 3.16)
project(ALLLSatisfiabilitySolver)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(ALLLSatisfiabilitySolver main.cpp
                                        core/VariablesArray.cpp
                                        core/VariablesArray.h
                                        core/Clause.cpp
                                        core/Clause.h
                                        core/RandomBoolGenerator.h
                                        sat_instance/SATInstance.cpp
                                        sat_instance/SATInstance.h
                                        sat_instance/SubSATInstance.cpp
                                        sat_instance/SubSATInstance.h
                                        cnf_io/cnf_io.cpp
                                        cnf_io/cnf_io.h)
target_link_libraries(ALLLSatisfiabilitySolver OpenMP::OpenMP_CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG") # to specify compiler optimisations add -O3 -DNDEBUG
