cmake_minimum_required(VERSION 3.0)
project(cell_pattern)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

# Boost
find_package(Boost 1.50.0 COMPONENTS numpy REQUIRED)
list(APPEND cell_pattern_LIBRARIES ${Boost_LIBRARIES})
list(APPEND cell_pattern_INCLUDES ${Boost_INCLUDE_DIRS})

# Eigen
add_definitions(-DHAVE_EIGEN)
list(APPEND cell_pattern_INCLUDES eigen)

# Python
add_subdirectory(pybind11)

# Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND cell_pattern_INCLUDES Aboria/src)
list(APPEND cell_pattern_INCLUDES Aboria/third-party)

add_library(cell_pattern MODULE cell_pattern.cpp GayBernePotential.cpp GayBernePotential.h)

target_link_libraries(cell_pattern PRIVATE pybind11::module)
target_include_directories(cell_pattern PRIVATE ${cell_pattern_INCLUDES})
set_target_properties(cell_pattern PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                           SUFFIX "${PYTHON_MODULE_EXTENSION}")
