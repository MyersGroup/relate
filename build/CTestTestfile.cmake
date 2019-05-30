# CMake generated Testfile for 
# Source directory: /data/desertfinch/speidel/relate
# Build directory: /data/desertfinch/speidel/relate/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/data/desertfinch/speidel/relate/bin/Tests")
set_tests_properties(UnitTest PROPERTIES  _BACKTRACE_TRIPLES "/data/desertfinch/speidel/relate/CMakeLists.txt;52;add_test;/data/desertfinch/speidel/relate/CMakeLists.txt;0;")
subdirs("include/src")
subdirs("include/src/gzstream")
subdirs("include/file_formats/tskit")
subdirs("include/test")
subdirs("include/pipeline")
subdirs("include/evaluate")
subdirs("include/treeview")
subdirs("include/file_formats")
subdirs("include/extract")
