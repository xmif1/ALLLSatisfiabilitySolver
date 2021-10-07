cmake_minimum_required(VERSION 3.17)
project(ALLLSatisfiabilitySolver)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(/usr/local/include)
#link_directories(/usr/local/lib)

add_executable(ALLLSatisfiabilitySolver main.cpp
                                        core/VariablesArray.cpp
                                        core/VariablesArray.h
                                        core/Clause.cpp
                                        core/Clause.h
                                        core/RandomBoolGenerator.h
                                        core/LaplacianLambdaCDP.h
                                        sat_instance/SATInstance.cpp
                                        sat_instance/SATInstance.h
                                        sat_instance/SubSATInstance.cpp
                                        sat_instance/SubSATInstance.h
                                        cnf_io/cnf_io.cpp
                                        cnf_io/cnf_io.h)
target_link_libraries(ALLLSatisfiabilitySolver OpenMP::OpenMP_CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse2 -DNDEBUG") # specify compiler optimisations
