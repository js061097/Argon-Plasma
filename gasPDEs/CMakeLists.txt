cmake_minimum_required(VERSION 3.10)
project(gasPDEs)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the paths to the SUNDIALS installation
set(SUNDIALS_DIR /Users/jstan/sundials-7.0.0)  # Adjust this path to your SUNDIALS installation directory
set(SUNDIALS_INCLUDE_DIR ${SUNDIALS_DIR}/include)
set(SUNDIALS_BUILD_INCLUDE_DIR ${SUNDIALS_DIR}/build/include)  # Include the build directory
set(SUNDIALS_LIB_DIR ${SUNDIALS_DIR}/build/src)  # Adjust this path to the actual library directory

# Include directories
include_directories(${SUNDIALS_INCLUDE_DIR} ${SUNDIALS_BUILD_INCLUDE_DIR} include)

# Source files
set(SOURCES
    src/main.cpp
    src/Inputs.cpp
    src/Processor.cpp
    src/Model.cpp
    src/Solver.cpp
    src/Reaction.cpp
)

# Add executable
add_executable(gasPDEs ${SOURCES})

# Link SUNDIALS libraries
target_link_libraries(gasPDEs ${SUNDIALS_LIB_DIR}/cvode/libsundials_cvode.a
                                 ${SUNDIALS_LIB_DIR}/nvector/serial/libsundials_nvecserial.a
                                 ${SUNDIALS_LIB_DIR}/sunlinsol/dense/libsundials_sunlinsoldense.a
                                 ${SUNDIALS_LIB_DIR}/sunmatrix/dense/libsundials_sunmatrixdense.a
                                 ${SUNDIALS_LIB_DIR}/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a
                                 ${SUNDIALS_LIB_DIR}/sunlinsol/band/libsundials_sunlinsolband.a
                                 ${SUNDIALS_LIB_DIR}/sunmatrix/band/libsundials_sunmatrixband.a
                                 ${SUNDIALS_LIB_DIR}/sundials/libsundials_core.a)

# Link additional libraries if required
target_link_libraries(gasPDEs m)  # Math library, if needed

# Additional settings
set_target_properties(gasPDEs PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    