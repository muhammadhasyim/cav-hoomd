# OpenMM CavityMD Plugin
# Cavity-molecule coupling forces for OpenMM

CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

project(OpenMMCavityMD)

# CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMM
find_package(OpenMM REQUIRED)
include_directories(${OpenMM_INCLUDE_DIR})

# Find CUDA (optional)
option(ENABLE_CUDA "Build CUDA platform" ON)
if(ENABLE_CUDA)
    find_package(CUDA QUIET)
    if(CUDA_FOUND)
        message(STATUS "CUDA found: ${CUDA_VERSION}")
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
    else()
        message(STATUS "CUDA not found, building without GPU support")
        set(ENABLE_CUDA OFF)
    endif()
endif()

# Plugin version
set(CAVITYMD_VERSION_MAJOR 1)
set(CAVITYMD_VERSION_MINOR 0)
set(CAVITYMD_VERSION_PATCH 0)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/openmmapi/include
    ${CMAKE_SOURCE_DIR}/platforms/reference/include
)

if(ENABLE_CUDA)
    include_directories(${CMAKE_SOURCE_DIR}/platforms/cuda/include)
endif()

# Source files - Core API
set(CAVITYMD_API_SOURCES
    openmmapi/src/CavityForce.cpp
    openmmapi/src/CavityForceImpl.cpp
)

set(CAVITYMD_API_HEADERS
    openmmapi/include/CavityForce.h
    openmmapi/include/CavityForceImpl.h
    openmmapi/include/CavityMDKernels.h
)

# Source files - Reference Platform
set(CAVITYMD_REF_SOURCES
    platforms/reference/src/ReferenceCavityForceKernel.cpp
    platforms/reference/src/ReferenceCavityMDKernelFactory.cpp
)

set(CAVITYMD_REF_HEADERS
    platforms/reference/include/ReferenceCavityForceKernel.h
    platforms/reference/include/ReferenceCavityMDKernelFactory.h
)

# Source files - CUDA Platform (if available)
if(ENABLE_CUDA)
    set(CAVITYMD_CUDA_SOURCES
        platforms/cuda/src/CudaCavityForceKernel.cpp
        platforms/cuda/src/CudaCavityMDKernelFactory.cpp
    )
    set(CAVITYMD_CUDA_KERNELS
        platforms/cuda/src/CudaCavityForceKernel.cu
    )
    set(CAVITYMD_CUDA_HEADERS
        platforms/cuda/include/CudaCavityForceKernel.h
        platforms/cuda/include/CudaCavityMDKernelFactory.h
    )
endif()

# Build shared library
add_library(OpenMMCavityMD SHARED
    ${CAVITYMD_API_SOURCES}
    ${CAVITYMD_API_HEADERS}
    ${CAVITYMD_REF_SOURCES}
    ${CAVITYMD_REF_HEADERS}
)

target_link_libraries(OpenMMCavityMD ${OpenMM_LIBRARIES})

# Build CUDA library if available
if(ENABLE_CUDA)
    add_library(OpenMMCavityMDCUDA SHARED
        ${CAVITYMD_CUDA_SOURCES}
        ${CAVITYMD_CUDA_HEADERS}
        ${CAVITYMD_CUDA_KERNELS}
    )
    target_link_libraries(OpenMMCavityMDCUDA OpenMMCavityMD ${OpenMM_LIBRARIES} ${CUDA_LIBRARIES})
    set_target_properties(OpenMMCavityMDCUDA PROPERTIES 
        CUDA_SEPARABLE_COMPILATION ON
    )
endif()

# Install targets
install(TARGETS OpenMMCavityMD DESTINATION lib)
install(FILES ${CAVITYMD_API_HEADERS} DESTINATION include/openmm-cavitymd)

if(ENABLE_CUDA)
    install(TARGETS OpenMMCavityMDCUDA DESTINATION lib)
endif()

# Python bindings (using SWIG)
option(BUILD_PYTHON_WRAPPERS "Build Python wrappers" ON)
if(BUILD_PYTHON_WRAPPERS)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(SWIG REQUIRED)
    include(${SWIG_USE_FILE})
    
    set(CMAKE_SWIG_FLAGS "")
    set_source_files_properties(python/cavitymd.i PROPERTIES CPLUSPLUS ON)
    
    swig_add_library(cavitymd
        TYPE SHARED
        LANGUAGE python
        SOURCES python/cavitymd.i
    )
    
    target_include_directories(cavitymd PRIVATE ${Python3_INCLUDE_DIRS})
    target_link_libraries(cavitymd OpenMMCavityMD ${Python3_LIBRARIES})
    
    install(TARGETS cavitymd DESTINATION ${Python3_SITELIB}/openmm/cavitymd)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cavitymd.py DESTINATION ${Python3_SITELIB}/openmm/cavitymd)
endif()

# Testing
option(BUILD_TESTING "Build test suite" ON)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Summary
message(STATUS "")
message(STATUS "OpenMM CavityMD Plugin Configuration:")
message(STATUS "  Version: ${CAVITYMD_VERSION_MAJOR}.${CAVITYMD_VERSION_MINOR}.${CAVITYMD_VERSION_PATCH}")
message(STATUS "  OpenMM: ${OpenMM_INCLUDE_DIR}")
message(STATUS "  CUDA: ${ENABLE_CUDA}")
message(STATUS "  Python: ${BUILD_PYTHON_WRAPPERS}")
message(STATUS "  Testing: ${BUILD_TESTING}")
message(STATUS "")
