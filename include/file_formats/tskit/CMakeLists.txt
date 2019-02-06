project (tskit)

add_library(kastoreStatic STATIC kastore)
add_library(kastoreShared SHARED kastore)

add_library(tskitStatic STATIC core tables trees genotypes convert stats tskit)

