cmake_minimum_required(VERSION 3.28)

project(chem-engine VERSION 0.2.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG        v2.5.3
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(doctest)

function(chem_apply_strict_warnings target)
  target_compile_options(${target} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:
      -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Werror>
  )
endfunction()

function(chem_apply_sanitizers target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(${target} PRIVATE -fsanitize=address,undefined)
  endif()
endfunction()

add_library(chem)

target_sources(chem
  PRIVATE
    src/chem/core/element.cpp
    src/chem/core/molecular_graph.cpp
    src/chem/core/composition.cpp
    src/chem/parsing/formula_parser.cpp
    src/chem/parsing/smiles_parser.cpp
    src/chem/canonical/canonical_smiles.cpp
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS src
    FILES
      src/chem/core/errors.hpp
      src/chem/core/element.hpp
      src/chem/core/molecular_graph.hpp
      src/chem/core/composition.hpp
      src/chem/parsing/formula_parser.hpp
      src/chem/parsing/smiles_parser.hpp
      src/chem/canonical/canonical_smiles.hpp
)

target_compile_features(chem PUBLIC cxx_std_20)
set_target_properties(chem PROPERTIES
  CXX_EXTENSIONS OFF
  CXX_STANDARD_REQUIRED ON
)

configure_file(data/elements.csv data/elements.csv COPYONLY)
target_compile_definitions(chem PRIVATE
  CHEM_ELEMENTS_CSV_PATH="${CMAKE_CURRENT_BINARY_DIR}/data/elements.csv"
)

configure_file(data/golden/corpus.csv data/golden/corpus.csv COPYONLY)

chem_apply_strict_warnings(chem)
chem_apply_sanitizers(chem)

add_executable(test-runner
  tests/main.cpp
  tests/core/element_test.cpp
  tests/core/molecular_graph_test.cpp
  tests/parsing/formula_parser_test.cpp
  tests/parsing/smiles_parser_test.cpp
  tests/canonical/canonical_smiles_test.cpp
  tests/canonical/golden_corpus_test.cpp
)

target_compile_features(test-runner PRIVATE cxx_std_20)
set_target_properties(test-runner PROPERTIES
  CXX_EXTENSIONS OFF
  CXX_STANDARD_REQUIRED ON
)

target_link_libraries(test-runner PRIVATE chem doctest::doctest)

target_compile_definitions(test-runner PRIVATE
  CHEM_GOLDEN_CORPUS_CSV_PATH="${CMAKE_CURRENT_BINARY_DIR}/data/golden/corpus.csv"
)

chem_apply_strict_warnings(test-runner)
chem_apply_sanitizers(test-runner)

enable_testing()
add_test(NAME chem_tests COMMAND test-runner)
