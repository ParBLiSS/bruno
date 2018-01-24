###### GCCFilt support, using my fixed version of gccfilter.
OPTION(ENABLE_COMPILER_PROFILING "Measure compilation time and memory usage.  Requires /usr/bin/time" OFF)


if (ENABLE_COMPILER_PROFILING)

    add_definitions(-ftime-report)

    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "/usr/bin/time --format 'COMPILE [TIME]: %e s. [MAX_RSS]: %M kB'")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "/usr/bin/time --format 'LINK [TIME]: %e s. [MAX_RSS]: %M kB'")

endif (ENABLE_COMPILER_PROFILING)
 
