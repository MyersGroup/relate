# CMake generated Testfile for 
# Source directory: /users/myers/speidel/Documents/relate
# Build directory: /users/myers/speidel/Documents/relate/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/users/myers/speidel/Documents/relate/bin/Tests")
subdirs("include/src")
subdirs("include/src/gzstream")
subdirs("include/test")
subdirs("include/pipeline")
subdirs("include/evaluate")
subdirs("include/treeview")
subdirs("include/file_formats")
subdirs("include/extract")
