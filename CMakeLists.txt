cmake_minimum_required(VERSION 3.16)
project(nano_geo_matrix VERSION 0.1.0 LANGUAGES CXX)

add_library(nano_geo_matrix INTERFACE)
# target_compile_definitions(
#   nano_geo_matrix
#   INTERFACE
#   NGM_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\"
# )
target_compile_features(nano_geo_matrix INTERFACE cxx_std_17)

target_include_directories(nano_geo_matrix INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prefer config mode (provides Eigen3::Eigen). Fall back to module mode for portability.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(Eigen3_FOUND)
  target_link_libraries(nano_geo_matrix INTERFACE Eigen3::Eigen)
else()
  find_package(Eigen3 3.3 REQUIRED)
  target_include_directories(nano_geo_matrix INTERFACE ${EIGEN3_INCLUDE_DIRS})
endif()

find_package(Armadillo REQUIRED)
target_include_directories(nano_geo_matrix INTERFACE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(nano_geo_matrix INTERFACE ${ARMADILLO_LIBRARIES})

find_package(GSL REQUIRED)

target_link_libraries(nano_geo_matrix
  INTERFACE
    GSL::gsl
    GSL::gslcblas
)

# ---- Bessel backend selection (after target exists!) ----
set(NGM_BESSEL_BACKEND "BL" CACHE STRING "Bessel backend: BL, COMPLEX_BESSEL, SPECIAL")
set_property(CACHE NGM_BESSEL_BACKEND PROPERTY STRINGS BL COMPLEX_BESSEL SPECIAL)

if(NGM_BESSEL_BACKEND STREQUAL "BL")
  target_compile_definitions(nano_geo_matrix INTERFACE USE_MYBESSEL_BL)

elseif(NGM_BESSEL_BACKEND STREQUAL "COMPLEX_BESSEL")
  target_compile_definitions(nano_geo_matrix INTERFACE USE_MYBESSEL)

  find_path(COMPLEX_BESSEL_INCLUDE_DIR complex_bessel.h PATHS /usr/local/include /usr/include)
  find_library(COMPLEX_BESSEL_LIBRARY complex_bessel PATHS /usr/local/lib /usr/lib /usr/lib/x86_64-linux-gnu)

  if(NOT COMPLEX_BESSEL_INCLUDE_DIR OR NOT COMPLEX_BESSEL_LIBRARY)
    message(FATAL_ERROR "complex_bessel not found. Set COMPLEX_BESSEL_INCLUDE_DIR and COMPLEX_BESSEL_LIBRARY.")
  endif()

  target_include_directories(nano_geo_matrix INTERFACE ${COMPLEX_BESSEL_INCLUDE_DIR})
  target_link_libraries(nano_geo_matrix INTERFACE ${COMPLEX_BESSEL_LIBRARY})

elseif(NGM_BESSEL_BACKEND STREQUAL "SPECIAL")
  target_compile_definitions(nano_geo_matrix INTERFACE USE_MYBESSEL_SPECIAL)

else()
  message(FATAL_ERROR "Unknown NGM_BESSEL_BACKEND=${NGM_BESSEL_BACKEND}")
endif()

option(NGM_BUILD_TESTS "Build nano_geo_matrix tests" ON)
if(NGM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
