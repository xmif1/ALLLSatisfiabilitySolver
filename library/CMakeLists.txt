cmake_minimum_required(VERSION 3.16)
project(ALLLSatisfiabilitySolver)
set(CMAKE_CXX_STANDARD 20)

set(DEFAULT_BUILD_TYPE "Release")
add_library(ALLLSatisfiabilitySolver INTERFACE)
target_include_directories(ALLLSatisfiabilitySolver INTERFACE
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

include(CMakeFindDependencyMacro)
find_dependency(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(ALLLSatisfiabilitySolver INTERFACE OpenMP::OpenMP_CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -DNDEBUG")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${PROJECT_BINARY_DIR}/ALLLSatisfiabilitySolverConfigVersion.cmake"
        VERSION 0.1
        COMPATIBILITY AnyNewerVersion
)

install(TARGETS ALLLSatisfiabilitySolver
        EXPORT ALLLSatisfiabilitySolverTargets
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
        RUNTIME DESTINATION bin COMPONENT Runtime
        PUBLIC_HEADER DESTINATION include COMPONENT Development
        BUNDLE DESTINATION bin COMPONENT Runtime)

include(CMakePackageConfigHelpers)
configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/ALLLSatisfiabilitySolverConfig.cmake.in"
        "${PROJECT_BINARY_DIR}/ALLLSatisfiabilitySolverConfig.cmake"
        INSTALL_DESTINATION lib/cmake/ALLLSatisfiabilitySolver
)

install(EXPORT ALLLSatisfiabilitySolverTargets DESTINATION lib/cmake/ALLLSatisfiabilitySolver)
install(FILES "${PROJECT_BINARY_DIR}/ALLLSatisfiabilitySolverConfigVersion.cmake"
        "${PROJECT_BINARY_DIR}/ALLLSatisfiabilitySolverConfig.cmake"
        DESTINATION lib/cmake/ALLLSatisfiabilitySolver)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)