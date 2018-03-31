
###### compile code for pipeline

add_executable(Relate ${PIPELINE_DIR}/Relate.cpp)
target_link_libraries(Relate relateStatic gzstreamStatic)

