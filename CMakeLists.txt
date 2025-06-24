CMAKE_MINIMUM_REQUIRED(VERSION 3.8 FATAL_ERROR)

# Name the plugin project.
# TODO: Set the project title to the name of your Python package.
project(hoomd-bussi-thermostat LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 17)  # Match HOOMD's C++17 standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)

# Find the installed HOOMD.
find_package(HOOMD REQUIRED)

message(STATUS "Found HOOMD ${HOOMD_VERSION}: ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}")
message(STATUS "HOOMD include path: ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}/include")
message(STATUS "HOOMD_SHORTREAL_SIZE: ${HOOMD_SHORTREAL_SIZE}")
message(STATUS "HOOMD_LONGREAL_SIZE: ${HOOMD_LONGREAL_SIZE}")

# Force installation to the HOOMD installation location.
set(CMAKE_INSTALL_PREFIX ${HOOMD_INSTALL_PREFIX} CACHE PATH "Installation prefix" FORCE)

# Add HOOMD preprocessor definitions
add_definitions(-DHOOMD_SHORTREAL_SIZE=${HOOMD_SHORTREAL_SIZE})
add_definitions(-DHOOMD_LONGREAL_SIZE=${HOOMD_LONGREAL_SIZE})

# Enable compiler warnings on gcc and clang.
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wno-sign-conversion -Wno-unknown-pragmas -Wno-deprecated-declarations -Wno-unused-result")
endif()

# Build bussi_reservoir module (thermostat only)
set(_bussi_reservoir_sources
    src/BussiReservoirThermostat.cc
    src/hoomd/bussi_reservoir/module.cc
    )

pybind11_add_module(_bussi_reservoir SHARED ${_bussi_reservoir_sources})
target_link_libraries(_bussi_reservoir PRIVATE HOOMD::_hoomd HOOMD::_md)
target_include_directories(_bussi_reservoir PRIVATE 
    ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}/include
    ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}/include/hoomd/extern/OpenRAND/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add the preprocessor definitions to the target
target_compile_definitions(_bussi_reservoir PRIVATE 
    HOOMD_SHORTREAL_SIZE=${HOOMD_SHORTREAL_SIZE}
    HOOMD_LONGREAL_SIZE=${HOOMD_LONGREAL_SIZE}
)

# Build cavitymd module (cavity force) - Re-enabled with CPU support
set(_cavitymd_sources
    src/CavityForceCompute.cc
    src/hoomd/cavitymd/module.cc
    )

# GPU support enabled for CUDA compilation
if (ENABLE_HIP)
list(APPEND _cavitymd_sources 
    src/CavityForceComputeGPU.cc
    src/CavityForceComputeGPU.cu
    )
endif()

pybind11_add_module(_cavitymd SHARED ${_cavitymd_sources})
target_link_libraries(_cavitymd PRIVATE HOOMD::_hoomd HOOMD::_md)
target_include_directories(_cavitymd PRIVATE 
    ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}/include
    ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}/include/hoomd/extern/OpenRAND/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add the preprocessor definitions to the target
target_compile_definitions(_cavitymd PRIVATE 
    HOOMD_SHORTREAL_SIZE=${HOOMD_SHORTREAL_SIZE}
    HOOMD_LONGREAL_SIZE=${HOOMD_LONGREAL_SIZE}
)

# Install the libraries
install(TARGETS _bussi_reservoir
        LIBRARY DESTINATION ${Python_SITEARCH}/hoomd/bussi_reservoir
        )

install(TARGETS _cavitymd
        LIBRARY DESTINATION ${Python_SITEARCH}/hoomd/cavitymd
        )

# Install the python packages
install(DIRECTORY src/hoomd/bussi_reservoir/
        DESTINATION ${Python_SITEARCH}/hoomd/bussi_reservoir
        FILES_MATCHING PATTERN "*.py"
        )

# Install the cavity module python packages
install(DIRECTORY src/hoomd/cavitymd/
        DESTINATION ${Python_SITEARCH}/hoomd/cavitymd
        FILES_MATCHING PATTERN "*.py"
        )

