cmake_minimum_required(VERSION 3.28)

project(chem-engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_STANDARD AND CMAKE_CXX_STANDARD LESS 17)
  message(FATAL_ERROR "Minimum supported C++ standard is C++17")
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

file(GLOB DATA_FILES "data/*")
foreach(data_file ${DATA_FILES})
  get_filename_component(filename ${data_file} NAME)
  configure_file(
    "data/${filename}"
    "data/${filename}"
    COPYONLY
  )
endforeach(data_file ${DATA_FILES})

file(GLOB CHEM_SRC "src/*.cpp")
file(GLOB CHEM_HEADERS "src/*.hpp")

add_library(chem)

target_sources(chem
  PRIVATE
    ${CHEM_SRC}
  
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS
      src
    FILES
      ${CHEM_HEADERS}
)

add_executable(test-runner)

file(GLOB CHEM_TESTS "test/*.cpp")

target_sources(test-runner
  PRIVATE
    ${CHEM_TESTS}
)

target_link_libraries(test-runner
  PRIVATE
    chem
)
