# CMake generated Testfile for 
# Source directory: /Users/leo/Documents/genomics/relate
# Build directory: /Users/leo/Documents/genomics/relate/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/Users/leo/Documents/genomics/relate/bin/Tests")
subdirs(include/src)
subdirs(include/src/gzstream)
subdirs(include/test)
subdirs(include/pipeline)
subdirs(include/evaluate)
subdirs(include/treeview)
subdirs(include/file_formats)
subdirs(include/extract)
