cmake_minimum_required(VERSION 3.7)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(PROJECT_VERSION 0.0.1)

project(HFINCPP VERSION ${PROJECT_VERSION}
        DESCRIPTION "Hartree-Fock code written in full C++ standard")

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    FIND_PACKAGE(Git)
    IF(GIT_FOUND)
        EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE "HFINCPP_BUILD_VERSION"
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        MESSAGE( STATUS "Git version: ${HFINCPP_BUILD_VERSION}" )
    ELSE(GIT_FOUND)
        SET(HFINCPP_BUILD_VERSION 0)
    ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/version/version.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/include/version.h
)

include(GNUInstallDirs)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/data)

enable_testing()

set(CMAKE_CXX_STANDARD 17)

add_definitions(-DBASIS_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/basis_set_exchange/\")

add_library(libHFINCPP
        src/run.cpp

        src/basis/basis.cpp

        src/geometry/geometry.cpp
        src/geometry/resolve.cpp

        src/gradient/driver.cpp
        src/gradient/gradient.cpp
        src/gradient/numerical.cpp
        src/gradient/utils.cpp

        src/integral/integral.cpp
        src/integral/obara_saika.cpp
        src/integral/rys_quadrature.cpp

        src/optimize/optimize.cpp

        src/scf/scf.cpp
        src/scf/setup.cpp
        src/scf/occupation.cpp

        src/hf/rhf.cpp
        
        src/mixing/simple_mixing.cpp
        )


add_executable(hfincpp
        src/main.cpp
)

target_link_libraries(hfincpp libHFINCPP)
####################
#  find Libraries  #
####################

#####Armadillo######

find_package(
        Armadillo
        8.500.0
        REQUIRED
)

if (${ARMADILLO_FOUND})
    message("Found Armadillo Library,")
    message("     ARMADILLO_VERSION: " ${ARMADILLO_VERSION_STRING})
    message("     ARMADILLO_LIBRARIES   : " ${ARMADILLO_LIBRARIES})
    message("     ARMADILLO_INCLUDE_DIRS: " ${ARMADILLO_INCLUDE_DIRS})

    include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(libHFINCPP ${ARMADILLO_LIBRARIES})
endif ()

find_package(
        HDF5
)

if (${HDF5_FOUND})
    message("Found HDF5 Library,")
    message("     HDF5_VERSION: " ${HDF5_VERSION})
    message("     HDF5_LIBRARIES   : " ${HDF5_LIBRARIES})
    message("     HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS})

    include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
    target_link_libraries(libHFINCPP ${HDF5_LIBRARIES})
endif ()
######OpenMP######

FIND_PACKAGE(OpenMP)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

######GSL######

find_package(
        GSL
        2.4
        REQUIRED
)

include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(libHFINCPP ${GSL_LIBRARIES})

######Boost######

find_package(
        Boost
        1.45.0
        REQUIRED
)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
target_link_libraries(libHFINCPP ${Boost_LIBRARIES})


################################
#  Compile external libraries  #
################################


#######Catch2#######

add_subdirectory(external/Catch2)
include_directories("external/Catch2/src")
target_link_libraries(libHFINCPP Catch2WithMain)

########fmt#########
add_subdirectory(external/fmt)
include_directories("external/fmt/include")
target_link_libraries(libHFINCPP fmt)
add_definitions(-DHAVE_FMT)

########args#########
add_subdirectory(external/args)
include_directories("external/args")
target_link_libraries(libHFINCPP args)
add_definitions(-DHAVE_ARGS)

########libcint#########
add_subdirectory(external/libcint)
include_directories("external/libcint/src")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/libcint/src)
target_link_libraries(libHFINCPP cint)
add_definitions(-DHAVE_LIBCINT)


#########json#########
add_library(json SHARED external/nlohmann/json.cpp)
include_directories("external/nlohmann")
target_include_directories(libHFINCPP PUBLIC external/nlohmann/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


####################
#       Test       #
####################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

add_executable(test_basis
        src/basis/basis_unit_test.cpp
        )

target_link_libraries(test_basis libHFINCPP)

add_executable(test_integral
        src/integral/obara_saika_unit_test.cpp
        src/integral/rys_quadrature_unit_test.cpp
        )

target_link_libraries(test_integral libHFINCPP)

add_executable(test_hf
        src/hf/rhf_unit_test.cpp)

target_link_libraries(test_hf libHFINCPP)