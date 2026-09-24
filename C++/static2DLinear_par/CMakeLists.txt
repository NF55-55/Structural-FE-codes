cmake_minimum_required(VERSION 3.10.0)
project(Static2DLinear VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MKL_INTERFACE "lp64")
set(MKL_THREADING "intel_thread")
set(MKL_LINK "dynamic")

find_package(MKL CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

file(GLOB SOURCES
    "src/*.cpp" 
    "src/*.h" 
)

add_executable(${PROJECT_NAME} ${SOURCES})

# To switch between Release/Debug -> F1 (command palette) -> CMake select 
# variant
add_compile_options( 
    -Wall 
    -Wsign-conversion 
    -fopenmp 
)

add_subdirectory(EigenLib/eigen-5.0.0/)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE
        Eigen3::Eigen 
        MKL::MKL 
        OpenMP::OpenMP_CXX 
)

# adding necessary preprocessor macros for Eigen-MKL coupling and 
# multi-threading
target_compile_definitions(${PROJECT_NAME} 
    PRIVATE
        EIGEN_USE_MKL_ALL)
