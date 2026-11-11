cmake_minimum_required(VERSION 3.20)
project(SRothmanEECs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(EEC_BUILD_TESTS "Build unit tests" ON)
option(EEC_BUILD_BENCHMARKS "Build benchmarks" ON)

# Inside a cmsenv'd shell, default external locations come from the CMSSW toolbox.
find_program(SCRAM_EXECUTABLE scram)
function(eec_scram_tag tool tag out)
  if(SCRAM_EXECUTABLE AND DEFINED ENV{CMSSW_BASE})
    execute_process(COMMAND ${SCRAM_EXECUTABLE} tool tag ${tool} ${tag}
                    WORKING_DIRECTORY $ENV{CMSSW_BASE}
                    OUTPUT_VARIABLE v OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    set(${out} "${v}" PARENT_SCOPE)
  endif()
endfunction()

eec_scram_tag(fastjet FASTJET_BASE _fj)
eec_scram_tag(boost_header INCLUDE _boost)
eec_scram_tag(catch2 INCLUDE _catch2)
eec_scram_tag(google-benchmark GOOGLE_BENCHMARK_BASE _gbench)

set(FASTJET_ROOT "${_fj}" CACHE PATH "fastjet install prefix")
set(BOOST_INCLUDE_DIR "${_boost}" CACHE PATH "Boost headers")
set(CATCH2_INCLUDE_DIR "${_catch2}" CACHE PATH "Directory containing catch.hpp (Catch2 v2)")
set(GBENCH_ROOT "${_gbench}" CACHE PATH "google-benchmark install prefix")

find_path(FASTJET_INCLUDE_DIR fastjet/PseudoJet.hh HINTS ${FASTJET_ROOT}/include REQUIRED)
find_library(FASTJET_LIBRARY fastjet HINTS ${FASTJET_ROOT}/lib REQUIRED)
find_path(BOOST_HISTOGRAM_INCLUDE boost/histogram.hpp HINTS ${BOOST_INCLUDE_DIR} REQUIRED)

# CMSSW include convention: "SRothman/EECs/interface/X.h".
set(EEC_STAGE ${CMAKE_BINARY_DIR}/include/SRothman)
file(MAKE_DIRECTORY ${EEC_STAGE})
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR} ${EEC_STAGE}/EECs SYMBOLIC)

add_library(EECs SHARED src/Binning.cc)
target_include_directories(EECs PUBLIC ${CMAKE_BINARY_DIR}/include)
target_include_directories(EECs SYSTEM PUBLIC ${FASTJET_INCLUDE_DIR} ${BOOST_HISTOGRAM_INCLUDE})
target_link_libraries(EECs PUBLIC ${FASTJET_LIBRARY})
target_compile_options(EECs PRIVATE -Wall -Wextra -Werror=missing-braces)

if(EEC_BUILD_TESTS)
  find_path(CATCH2_HEADER catch.hpp HINTS ${CATCH2_INCLUDE_DIR} REQUIRED)
  enable_testing()
  file(GLOB EEC_TEST_SOURCES CONFIGURE_DEPENDS test/*.cc)
  add_executable(testEECs ${EEC_TEST_SOURCES})
  target_include_directories(testEECs SYSTEM PRIVATE ${CATCH2_HEADER})
  target_link_libraries(testEECs PRIVATE EECs)
  target_compile_options(testEECs PRIVATE -Wall -Wextra -Werror=missing-braces)
  add_test(NAME testEECs COMMAND testEECs)

  # Each case must fail to compile with output matching its regex (GCC or clang wording).
  set(EEC_CF_FLAGS -std=c++20 -fsyntax-only -I${CMAKE_BINARY_DIR}/include
      -isystem ${FASTJET_INCLUDE_DIR} -isystem ${BOOST_HISTOGRAM_INCLUDE})
  set(EEC_CF_CASES  # pairs: case name, pass regex
      pairbins_unbinned "precompute::PairBins requires a binned output"
      update_extra_missing "Extra has no update\\(\\) accepting these arguments"
      update_extra_badargs "Extra has no update\\(\\) accepting these arguments"
      concept_bad_weight "satisfaction of 'WeightPolicy<|does not satisfy 'WeightPolicy'"
      concept_bad_frame "satisfaction of 'FramePolicy<|does not satisfy 'FramePolicy'"
      concept_bad_extra_runtime "satisfaction of 'ExtraPolicy<Extra, 3>'|'ExtraPolicy<PairOnlyExtra, 3>' evaluated to false")
  list(LENGTH EEC_CF_CASES n)
  math(EXPR last "${n} - 1")
  foreach(i RANGE 0 ${last} 2)
    math(EXPR j "${i} + 1")
    list(GET EEC_CF_CASES ${i} name)
    list(GET EEC_CF_CASES ${j} msg)
    add_test(NAME compile_fail_${name}
             COMMAND ${CMAKE_CXX_COMPILER} ${EEC_CF_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/test/compile_fail/${name}.cc)
    set_tests_properties(compile_fail_${name} PROPERTIES PASS_REGULAR_EXPRESSION "${msg}")
  endforeach()
endif()

if(EEC_BUILD_BENCHMARKS)
  find_package(benchmark REQUIRED HINTS ${GBENCH_ROOT}/lib64/cmake/benchmark ${GBENCH_ROOT}/lib/cmake/benchmark)
  add_executable(eec_bench bench/bench_penc.cc)
  target_link_libraries(eec_bench PRIVATE EECs benchmark::benchmark)
  add_executable(eec_bench_cov bench/bench_cov.cc)
  target_link_libraries(eec_bench_cov PRIVATE EECs benchmark::benchmark)
endif()
