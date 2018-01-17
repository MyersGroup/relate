project                     (Relate)
cmake_minimum_required      (VERSION 3.1)

#if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
  #endif()
set (CMAKE_CXX_STANDARD 11)

# Set dir
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/src)
set(PIPELINE_DIR ${CMAKE_SOURCE_DIR}/include/pipeline)
set(EVALUATE_DIR ${CMAKE_SOURCE_DIR}/include/evaluate)
set(FILEFORMATS_DIR ${CMAKE_SOURCE_DIR}/include/file_formats)
set(EXTRACT_DIR ${CMAKE_SOURCE_DIR}/include/extract)
set(MISC_DIR ${CMAKE_SOURCE_DIR}/include/misc)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/include/test)

# Set includes
include_directories(${PROJECT_SOURCE_DIR} ${TEST_DIR})

# Add subdirectories
add_subdirectory(${PROJECT_SOURCE_DIR})
add_subdirectory(${TEST_DIR})
add_subdirectory(${PIPELINE_DIR})
add_subdirectory(${EVALUATE_DIR})
add_subdirectory(${FILEFORMATS_DIR})
add_subdirectory(${EXTRACT_DIR})

# tests
add_executable(Tests ${TEST_DIR}/Tests.cpp)
target_link_libraries(Tests relateStatic UnitTest)

enable_testing()
add_test(UnitTest ${CMAKE_BINARY_DIR}/Tests)
add_custom_command(TARGET Tests COMMENT "Run tests" POST_BUILD WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND Tests)
