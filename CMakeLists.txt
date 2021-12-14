cmake_minimum_required(VERSION 3.10)

project(c3c4free_proj VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(c3c4free ${PROJECT_SOURCES}
  src/main.cpp
  src/c3c4free.cpp
  src/io.cpp
  src/vertex_impl.cpp
  )

target_include_directories(c3c4free PRIVATE src)

target_include_directories(c3c4free PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/lib/boost/boost/
  src/
  )

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -O3 -pthread")

