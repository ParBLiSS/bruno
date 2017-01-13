set(TEST_NAME bruno_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")



if (BUILD_TEST_APPLICATIONS)

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

endif(BUILD_TEST_APPLICATIONS)

#benchmarkers
if (BUILD_EXAMPLE_APPLICATIONS)


# 51 added because hipmer uses that
foreach(k 21 31 51)

#no DNA5 - fails when searching for matching kmers.
foreach(dna 4)


add_executable(compact_debruijn_graph_fastq_A${dna}_K${k} compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}
     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k}")
target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k} ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_A${dna}_K${k} compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}
     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k}")
target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k} ${EXTRA_LIBS})

endforeach(dna)

add_executable(compact_debruijn_graph_fastq_A4_K${k}_minmem compact_debruijn_graph_block_construct.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A4_K${k}_minmem
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=${k}")
target_link_libraries(compact_debruijn_graph_fastq_A4_K${k}_minmem ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_A4_K${k}_minmem compact_debruijn_graph_block_construct.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A4_K${k}_minmem
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA -DpK=${k}")
target_link_libraries(compact_debruijn_graph_fasta_A4_K${k}_minmem ${EXTRA_LIBS})


endforeach(k)


endif(BUILD_EXAMPLE_APPLICATIONS)

