# Require a minimum version of cmake
cmake_minimum_required(VERSION 3.1)
# set the name of the project
set(PROJECT_NAME_STR rminimaxlib)
project(${PROJECT_NAME_STR})

set(CMAKE_BUILD_TYPE Release)

file(GLOB PROJECT_SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../rminimax/)

include_directories(${PROJECT_INCLUDE_DIR} ${DEP_INCLUDE_DIRS})
add_library(rminimax SHARED ${PROJECT_SRC_FILES})

target_link_libraries(rminimax PRIVATE ${LIB_DEPS})

install(TARGETS rminimax DESTINATION lib)
install(DIRECTORY ${PROJECT_INCLUDE_DIR} DESTINATION include)