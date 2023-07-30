cmake_minimum_required(VERSION 3.26.4)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

project(tests)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.55 REQUIRED COMPONENTS unit_test_framework)

file(GLOB test_SRC "*.cpp")

add_executable(${PROJECT_NAME} ${test_SRC})
# add_executable(${PROJECT_NAME} 
# ./test_fm_index.cpp
# )
target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
# target_link_libraries(${PROJECT_NAME} example ${Boost_LIBRARIES})

# add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${PROJECT_NAME})