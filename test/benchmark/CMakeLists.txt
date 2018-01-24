set(TEST_NAME bruno_benchmark_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


# available params:


if (ENABLE_BENCHMARKING)



# 51 added because hipmer uses that
foreach(k 21 31 51 55 63)

    add_executable(benchmark_compact_dbg_fastq_A4_K${k}_minmem benchmark_compact_dbg_block_construct.cpp)
    SET_TARGET_PROPERTIES(benchmark_compact_dbg_fastq_A4_K${k}_minmem
        PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTQ -DpK=${k} -DMIN_MEM")
    target_link_libraries(benchmark_compact_dbg_fastq_A4_K${k}_minmem ${EXTRA_LIBS})

    add_executable(benchmark_compact_dbg_fasta_A4_K${k}_minmem benchmark_compact_dbg_block_construct.cpp)
    SET_TARGET_PROPERTIES(benchmark_compact_dbg_fasta_A4_K${k}_minmem
        PROPERTIES COMPILE_FLAGS "-DpDNA=4 -DpPARSER=FASTA -DpK=${k} -DMIN_MEM")
    target_link_libraries(benchmark_compact_dbg_fasta_A4_K${k}_minmem ${EXTRA_LIBS})


endforeach(k)



## EXECUTABLES
add_executable(profile_mem_usage profile_mem_usage.cpp) 
target_link_libraries(profile_mem_usage ${EXTRA_LIBS})


if (Boost_FOUND)
	add_executable(cust_alloc test_custom_allocator.cpp)
	target_link_libraries(cust_alloc ${EXTRA_LIBS})
endif(Boost_FOUND)

add_executable(test_prealloc test_prealloc.cpp) 
target_link_libraries(test_prealloc ${EXTRA_LIBS})


endif(ENABLE_BENCHMARKING)

