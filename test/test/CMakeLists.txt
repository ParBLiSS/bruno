set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")



if (BUILD_EXAMPLE_APPLICATIONS)

#----------------
# COMMENTED OUT because the template parameters are not yet updated.
# also BenchmarkKmerIndex is serving the same purpose.
#add_executable(testKmerIndex_ST2 TestKmerIndexAndQuerySingleThread_DistMap.cpp)
#target_link_libraries(testKmerIndex_ST2 ${EXTRA_LIBS})

#add_executable(testKmerIndex_read TestKmerIndex_ParseOnly.cpp)
#target_link_libraries(testKmerIndex_read ${EXTRA_LIBS})

#add_executable(testKmerIndex_read_fasta TestKmerIndex_FASTAParseOnly.cpp)
#target_link_libraries(testKmerIndex_read_fasta ${EXTRA_LIBS})
#-----------------

#add_executable(test_de_bruijn_graph_DNA4 test_de_bruijn_graph_construction.cpp)
#SET_TARGET_PROPERTIES(test_de_bruijn_graph_DNA4
#     PROPERTIES COMPILE_FLAGS "-DpDNA=4")
#target_link_libraries(test_de_bruijn_graph_DNA4 ${EXTRA_LIBS})

add_executable(test_de_bruijn_graph_DNA5 test_de_bruijn_graph_construction.cpp)
SET_TARGET_PROPERTIES(test_de_bruijn_graph_DNA5
     PROPERTIES COMPILE_FLAGS "-DpDNA=5")
target_link_libraries(test_de_bruijn_graph_DNA5 ${EXTRA_LIBS})

add_executable(test_de_bruijn_graph_DNA16 test_de_bruijn_graph_construction.cpp)
SET_TARGET_PROPERTIES(test_de_bruijn_graph_DNA16
     PROPERTIES COMPILE_FLAGS "-DpDNA=16")
target_link_libraries(test_de_bruijn_graph_DNA16 ${EXTRA_LIBS})


add_executable(compact_debruijn_graph_fastq compact_debruijn_graph.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ")
target_link_libraries(compact_debruijn_graph_fastq ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta compact_debruijn_graph.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTA")
target_link_libraries(compact_debruijn_graph_fasta ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_low_mem_fastq compact_debruijn_graph_low_mem.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_low_mem_fastq
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ")
target_link_libraries(compact_debruijn_graph_low_mem_fastq ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_low_mem_fasta compact_debruijn_graph_low_mem.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_low_mem_fasta
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTA")
target_link_libraries(compact_debruijn_graph_low_mem_fasta ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_refactored_fastq compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_refactored_fastq
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ")
target_link_libraries(compact_debruijn_graph_refactored_fastq ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_refactored_fasta compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_refactored_fasta
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTA")
target_link_libraries(compact_debruijn_graph_refactored_fasta ${EXTRA_LIBS})

#benchmarkers

add_executable(compact_debruijn_graph_fastq_A4_K15 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A4_K15
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=15")
target_link_libraries(compact_debruijn_graph_fastq_A4_K15 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A16_K15 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A16_K15
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ -DpK=15")
target_link_libraries(compact_debruijn_graph_fastq_A16_K15 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A4_K21 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A4_K21
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=21")
target_link_libraries(compact_debruijn_graph_fastq_A4_K21 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A16_K21 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A16_K21
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ -DpK=21")
target_link_libraries(compact_debruijn_graph_fastq_A16_K21 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A4_K31 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A4_K31
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=31")
target_link_libraries(compact_debruijn_graph_fastq_A4_K31 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A16_K31 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A16_K31
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ -DpK=31")
target_link_libraries(compact_debruijn_graph_fastq_A16_K31 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A4_K63 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A4_K63
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=63")
target_link_libraries(compact_debruijn_graph_fastq_A4_K63 ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_A16_K63 compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A16_K63
     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ -DpK=63")
target_link_libraries(compact_debruijn_graph_fastq_A16_K63 ${EXTRA_LIBS})



add_executable(test_get_file_size test_get_file_size.cpp)
target_link_libraries(test_get_file_size ${EXTRA_LIBS})

endif()

