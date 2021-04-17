cmake_minimum_required(VERSION 3.20)

project(lina_unit_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

add_compile_definitions(GTEST_LANG_CXX11=1)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.8.1
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/common>)

target_link_libraries(${PROJECT_NAME}  gtest_main lina)

enable_testing()

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})
