cmake_minimum_required(VERSION 3.16)
project(ALLLSatisfiabilitySolver)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(ALLLSatisfiabilitySolver main.cpp
               core/VariablesArray.h
               core/Clause.h
               core/RandomBoolGenerator.h
               sat_instance/SATInstance.h
               cnf_io/cnf_io.cpp
               cnf_io/cnf_io.h)
target_link_libraries(ALLLSatisfiabilitySolver
                      OpenMP::OpenMP_CXX
                      Boost::program_options)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -DNDEBUG") # to specify compiler optimisations add -O3 -DNDEBUG
