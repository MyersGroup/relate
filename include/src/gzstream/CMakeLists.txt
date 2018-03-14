project (gzstream)

add_library(gzstreamStatic STATIC gzstream)
target_link_libraries(gzstreamStatic ${ZLIB_LIBRARIES})
add_library(gzstreamShared SHARED gzstream)
target_link_libraries(gzstreamShared ${ZLIB_LIBRARIES})
