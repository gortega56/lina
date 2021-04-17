cmake_minimum_required(VERSION 3.20)

project(lina VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

add_library(${PROJECT_NAME} INTERFACE)

set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/vector.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/matrix.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/include/quaternion.h
)

target_sources(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${SOURCE_FILES}>")
target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>)
target_include_directories(${PROJECT_NAME} SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
add_custom_target(linalib SOURCES ${SOURCE_FILES})

enable_testing()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_tests)