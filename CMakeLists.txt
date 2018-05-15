cmake_minimum_required(VERSION 2.6)
project(NERSTRAND)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if (NOT DEFINED WILDRIVER_PATH)
  set(WILDRIVER_PATH "wildriver")
endif()

if (DEFINED DEBUG OR DEFINED ASSERT)
  if (DEFINED DEBUG)
    add_definitions(-DDEBUG=${DEBUG})
    message("Debugging enabled")
  endif()
  if (DEFINED ASSERT)
    add_definitions(-DUSE_ASSERTS=${ASSERT})
    message("Assertions enabled")
  endif()
  set(CMAKE_BUILD_TYPE Debug)
  message("Building a debug release")
else()
  set(CMAKE_BUILD_TYPE Release)
  add_definitions(-DNDEBUG=1)
  message("Building an optimized release")
endif()

if (STATIC)
  set(NERSTRAND_LIBRARY_TYPE STATIC)
else()
  set(NERSTRAND_LIBRARY_TYPE SHARED)
endif()

if (DEFINED BIGEDGES AND BIGEDGES)
  add_definitions(-DNERSTRAND_64BIT_EDGES=${BIGEDGES})
  add_definitions(-DWILDRIVER_INDEX_TYPE=uint64_t)
else()
  add_definitions(-DWILDRIVER_INDEX_TYPE=uint32_t)
endif()

if (DEFINED BIGVERTICES AND BIGVERTICES)
  add_definitions(-DNERSTRAND_64BIT_VERTICES=${BIGVERTICES})
  add_definitions(-DWILDRIVER_DIMENSION_TYPE=uint64_t)
else()
  add_definitions(-DWILDRIVER_DIMENSION_TYPE=uint32_t)
endif()

if (DEFINED DOUBLEWEIGHTS AND DOUBLEWEIGHTS)
  add_definitions(-DNERSTRAND_DOUBLE_WEIGHTS=${DOUBLEWEIGHTS})
  add_definitions(-DWILDRIVER_VALUE_TYPE=double)
else()
  add_definitions(-DWILDRIVER_VALUE_TYPE=float)
endif()

set(REQ_LIBS ${REQ_LIBS} m)

if (DEFINED DEVEL AND NOT DEVEL EQUAL 0)
  message("Development mode enabled")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wmissing-prototypes") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wmissing-declarations") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-aliasing") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic") 
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99") 

# build wildriver
add_subdirectory("${WILDRIVER_PATH}")

# use gnu directories
if (NOT WIN32)
  # use gnu directories
  include(GNUInstallDirs)
else()
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()



set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${DOMLIB_PATH}")
include_directories("${WILDRIVER_PATH}/include")
add_subdirectory("include")
add_subdirectory("src")

