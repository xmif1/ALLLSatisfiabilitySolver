cmake_minimum_required(VERSION 3.16)
project(ALLLSatisfiabilitySolverMain)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(ALLLSatisfiabilitySolver CONFIG REQUIRED)

add_executable(${PROJECT_NAME}
               main.cpp
               cnf_io/cnf_io.cpp
               cnf_io/cnf_io.h)
target_link_libraries(${PROJECT_NAME} Boost::program_options ALLLSatisfiabilitySolver)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -DNDEBUG") # to specify compiler optimisations add -Ofast -DNDEBUG
