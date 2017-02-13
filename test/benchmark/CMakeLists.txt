set(TEST_NAME bruno_benchmark_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


# available params:


if (BUILD_EXAMPLE_APPLICATIONS)

## EXECUTABLES
add_executable(profile_mem_usage profile_mem_usage.cpp) 
target_link_libraries(profile_mem_usage ${EXTRA_LIBS})


if (Boost_FOUND)
	add_executable(cust_alloc test_custom_allocator.cpp)
	target_link_libraries(cust_alloc ${EXTRA_LIBS})
endif(Boost_FOUND)

add_executable(test_prealloc test_prealloc.cpp) 
target_link_libraries(test_prealloc ${EXTRA_LIBS})


endif(BUILD_EXAMPLE_APPLICATIONS)

