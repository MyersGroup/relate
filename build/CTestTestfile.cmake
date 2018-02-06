# CMake generated Testfile for 
# Source directory: /homes/speidel/speidel/genomics/relate
# Build directory: /homes/speidel/speidel/genomics/relate/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/homes/speidel/speidel/genomics/relate/bin/Tests")
subdirs("include/src")
subdirs("include/test")
subdirs("include/pipeline")
subdirs("include/evaluate")
subdirs("include/file_formats")
subdirs("include/extract")
