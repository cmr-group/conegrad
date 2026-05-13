cmake_minimum_required(VERSION 3.26...3.30)
project(conegrad LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# MSVC: make C source more forgiving toward the EPIC-derived code.
if(MSVC)
  add_compile_options(/wd4244 /wd4267 /wd4305 /wd4101 /wd4018 /wd4715 /wd4477)
  # _CRT_SECURE_NO_WARNINGS silences sprintf/fopen deprecation warnings.
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
  add_compile_options(-Wno-unused-result -Wno-unused-but-set-variable
                      -Wno-unused-variable -Wno-implicit-function-declaration
                      -Wno-format)
endif()

find_package(Python 3.11 COMPONENTS Interpreter Development.Module Development.SABIModule REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
  _conegrad
  STABLE_ABI
  NB_STATIC
  src/conegrad/_bindings.cpp
  src/conegrad/_c/scones_design.c
)

target_compile_definitions(_conegrad PRIVATE STANDALONE=1)
target_include_directories(_conegrad PRIVATE src/conegrad/_c)

install(TARGETS _conegrad LIBRARY DESTINATION conegrad)
