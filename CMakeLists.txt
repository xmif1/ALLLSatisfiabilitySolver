cmake_minimum_required(VERSION 3.17)
project(ALLLSatisfiabilitySolver)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(ALLLSatisfiabilitySolver main.cpp
                                        core/VariablesArray.cpp
                                        core/VariablesArray.h
                                        core/Clause.cpp
                                        core/Clause.h
                                        core/RandomBoolGenerator.h
                                        core/LaplacianLambdaCDP.h
                                        sat_instance/SATInstance.cpp
                                        sat_instance/SATInstance.h)
target_link_libraries(ALLLSatisfiabilitySolver -lcnf_io)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse2 -DNDEBUG") # specify compiler optimisations
