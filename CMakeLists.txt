cmake_minimum_required(VERSION 2.8)
project(Xolik CXX)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default Build Type is Release" FORCE)

# build unittest
include_directories(src/
                    3rdparty/googletest-release-1.7.0/include/
                    3rdparty/tclap-1.2.1/include/
                    lib/)
set(gtest_disable_pthreads ON CACHE BOOL "disable pthreads")
add_subdirectory(3rdparty/googletest-release-1.7.0)

add_definitions("-std=c++1y")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if (MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
endif()

set(CMAKE_DEBUG_POSTFIX "d")

add_executable(Xolik src/cli.cpp)
if (MSVC)
    target_link_libraries(Xolik debug ${CMAKE_SOURCE_DIR}/lib/mzloaderd.lib
                                debug ${CMAKE_SOURCE_DIR}/lib/ppdatad.lib
                                optimized ${CMAKE_SOURCE_DIR}/lib/mzloader.lib
                                optimized ${CMAKE_SOURCE_DIR}/lib/ppdata.lib)
elseif(UNIX)
    target_link_libraries(Xolik debug ${CMAKE_SOURCE_DIR}/lib/libmzloader.a
                                debug ${CMAKE_SOURCE_DIR}/lib/libppdata.a
                                optimized ${CMAKE_SOURCE_DIR}/lib/libmzloader.a
                                optimized ${CMAKE_SOURCE_DIR}/lib/libppdata.a)
endif()

add_executable(unittest test/RegressionTest.cpp)
target_link_libraries(unittest gtest gtest_main)