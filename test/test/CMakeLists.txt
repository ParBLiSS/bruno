set(TEST_NAME bruno_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")



if (BUILD_TEST_APPLICATIONS)

# need support for gap or unknown characters.
# add_executable(test_de_bruijn_graph_DNA test_de_bruijn_graph_construction.cpp)
# SET_TARGET_PROPERTIES(test_de_bruijn_graph_DNA
#      PROPERTIES COMPILE_FLAGS "-DpDNA=4")
# target_link_libraries(test_de_bruijn_graph_DNA ${EXTRA_LIBS})
#
#add_executable(test_de_bruijn_graph_DNA5 test_de_bruijn_graph_construction.cpp)
#SET_TARGET_PROPERTIES(test_de_bruijn_graph_DNA5
#     PROPERTIES COMPILE_FLAGS "-DpDNA=5")
#target_link_libraries(test_de_bruijn_graph_DNA5 ${EXTRA_LIBS})
#
#add_executable(test_de_bruijn_graph_DNA16 test_de_bruijn_graph_construction.cpp)
#SET_TARGET_PROPERTIES(test_de_bruijn_graph_DNA16
#     PROPERTIES COMPILE_FLAGS "-DpDNA=16")
#target_link_libraries(test_de_bruijn_graph_DNA16 ${EXTRA_LIBS})


#add_executable(compact_debruijn_graph_fastq compact_debruijn_graph.cpp)
#SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq
#     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ")
#target_link_libraries(compact_debruijn_graph_fastq ${EXTRA_LIBS})
#
#add_executable(compact_debruijn_graph_fasta compact_debruijn_graph.cpp)
#SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta
#     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTA")
#target_link_libraries(compact_debruijn_graph_fasta ${EXTRA_LIBS})

#add_executable(compact_debruijn_graph_low_mem_fastq compact_debruijn_graph_low_mem.cpp)
#SET_TARGET_PROPERTIES(compact_debruijn_graph_low_mem_fastq
#     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTQ")
#target_link_libraries(compact_debruijn_graph_low_mem_fastq ${EXTRA_LIBS})
#
#add_executable(compact_debruijn_graph_low_mem_fasta compact_debruijn_graph_low_mem.cpp)
#SET_TARGET_PROPERTIES(compact_debruijn_graph_low_mem_fasta
#     PROPERTIES COMPILE_FLAGS "-DpDNA=16 -DpPARSER=FASTA")
#target_link_libraries(compact_debruijn_graph_low_mem_fasta ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fastq_refactored compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_refactored
    PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ")
target_link_libraries(compact_debruijn_graph_fastq_refactored ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_refactored compact_debruijn_graph_refactor.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_refactored
    PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA")
target_link_libraries(compact_debruijn_graph_fasta_refactored ${EXTRA_LIBS})


add_executable(compact_debruijn_graph_fastq_blocked compact_debruijn_graph_block_construct.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_blocked
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ")
target_link_libraries(compact_debruijn_graph_fastq_blocked ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_blocked compact_debruijn_graph_block_construct.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_blocked
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA")
target_link_libraries(compact_debruijn_graph_fasta_blocked ${EXTRA_LIBS})


add_executable(compact_debruijn_graph_fastq_blocked_incr compact_debruijn_graph_block_construct.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_blocked_incr
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DMIN_MEM")
target_link_libraries(compact_debruijn_graph_fastq_blocked_incr ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_blocked_incr compact_debruijn_graph_block_construct.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_blocked_incr
     PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA -DMIN_MEM")
target_link_libraries(compact_debruijn_graph_fasta_blocked_incr ${EXTRA_LIBS})


add_executable(compact_debruijn_graph_fastq_freq compact_dbg_freq_construct_clean.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_freq
    PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=31 -DNCLEAN")
target_link_libraries(compact_debruijn_graph_fastq_freq ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_freq compact_dbg_freq_construct_clean.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_freq
    PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA -DpK=31 -DNCLEAN")
target_link_libraries(compact_debruijn_graph_fasta_freq ${EXTRA_LIBS})


add_executable(compact_debruijn_graph_fastq_freq_incr compact_dbg_freq_construct_clean.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_freq_incr
    PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=31 -DMIN_MEM -DNCLEAN")
target_link_libraries(compact_debruijn_graph_fastq_freq_incr ${EXTRA_LIBS})

add_executable(compact_debruijn_graph_fasta_freq_incr compact_dbg_freq_construct_clean.cpp)
SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_freq_incr
    PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA -DpK=31 -DMIN_MEM -DNCLEAN")
target_link_libraries(compact_debruijn_graph_fasta_freq_incr ${EXTRA_LIBS})


endif(BUILD_TEST_APPLICATIONS)

#benchmarkers
if (BUILD_EXAMPLE_APPLICATIONS)


# 51 added because hipmer uses that
foreach(k 31)

    #no DNA5 - fails when searching for matching kmers.
    foreach(dna 4 16)
        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k} compact_debruijn_graph_block_construct.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k}")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k} ${EXTRA_LIBS})

        # add_executable(compact_debruijn_graph_fasta_A${dna}_K${k} compact_debruijn_graph_block_construct.cpp)
        # SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}
        #     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k}")
        # target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k} ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DNCLEAN")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq ${EXTRA_LIBS})

        # add_executable(compact_debruijn_graph_fasta_A${dna}_K${k}_freq compact_debruijn_graph_block_construct.cpp)
        # SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}_freq
        #     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k}")
        # target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k}_freq ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_minimizer compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_minimizer
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMINIMIZER -DNCLEAN")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_minimizer ${EXTRA_LIBS})



        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_incr compact_debruijn_graph_block_construct.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_incr ${EXTRA_LIBS})

        # add_executable(compact_debruijn_graph_fasta_A${dna}_K${k}_incr compact_debruijn_graph_block_construct.cpp)
        # SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}_incr
        #     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k} -DMIN_MEM")
        # target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k}_incr ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_incr compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM -DNCLEAN")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_incr ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_minimizer_incr compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_minimizer_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM -DMINIMIZER -DNCLEAN")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_minimizer_incr ${EXTRA_LIBS})


    endforeach(dna)

endforeach(k)

foreach(k 21 51 55 63)

    #no DNA5 - fails when searching for matching kmers.
    foreach(dna 4)
        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k} compact_debruijn_graph_block_construct.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k}")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k} ${EXTRA_LIBS})

        # add_executable(compact_debruijn_graph_fasta_A${dna}_K${k} compact_debruijn_graph_block_construct.cpp)
        # SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}
        #     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k}")
        # target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k} ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DNCLEAN")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq ${EXTRA_LIBS})

        # add_executable(compact_debruijn_graph_fasta_A${dna}_K${k}_freq compact_debruijn_graph_block_construct.cpp)
        # SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}_freq
        #     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k}")
        # target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k}_freq ${EXTRA_LIBS})



        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_incr compact_debruijn_graph_block_construct.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_incr ${EXTRA_LIBS})

        # add_executable(compact_debruijn_graph_fasta_A${dna}_K${k}_incr compact_debruijn_graph_block_construct.cpp)
        # SET_TARGET_PROPERTIES(compact_debruijn_graph_fasta_A${dna}_K${k}_incr
        #     PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTA -DpK=${k} -DMIN_MEM")
        # target_link_libraries(compact_debruijn_graph_fasta_A${dna}_K${k}_incr ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_incr compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM -DNCLEAN")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_incr ${EXTRA_LIBS})

    endforeach(dna)

endforeach(k)

foreach(k 31)

    #no DNA5 - fails when searching for matching kmers.
    foreach(dna 4)
        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k}")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_recompact compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_recompact
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DRECOMPACT")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_recompact ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_incr compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_incr ${EXTRA_LIBS})

        add_executable(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_recompact_incr compact_dbg_freq_construct_clean.cpp)
        SET_TARGET_PROPERTIES(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_recompact_incr
            PROPERTIES COMPILE_FLAGS "-DpDNA=${dna} -DpPARSER=FASTQ -DpK=${k} -DRECOMPACT -DMIN_MEM")
        target_link_libraries(compact_debruijn_graph_fastq_A${dna}_K${k}_freq_clean_recompact_incr ${EXTRA_LIBS})


    endforeach(dna)

endforeach(k)


endif(BUILD_EXAMPLE_APPLICATIONS)

