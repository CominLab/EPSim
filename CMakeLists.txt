cmake_minimum_required(VERSION 2.8)
project(Thesis)
#find_package(seqan)

set (SEQAN_INCLUDE_DIR 
    /home/ciompin/PairwiseEP/seqan-trunk/core/include 
    /home/ciompin/PairwiseEP/seqan-trunk/extras/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(include
                   ${SEQAN_INCLUDE_DIR})

add_library(${PROJECT_NAME}
  src/stat_models.cpp include/stat_models.h
  src/kernel_gen.cpp include/kernel_gen.h
  src/entropies.cpp include/entropies.h
  src/score_filler.cpp include/score_filler.h
)

# TESI
add_executable (ep_sim src/ep_sim.cpp)
target_link_libraries(ep_sim
                      ${PROJECT_NAME})

# per impiantare motivi nelle seq reali
add_executable (real_seq src/real_seq.cpp)

# per estrarre sequenze
add_executable (seq_splitter src/seq_splitter.cpp)

# per i grafici
add_executable (ppv_plotter src/ppv_plotter.cpp)
add_executable (ppv_compare src/ppv_compare.cpp)