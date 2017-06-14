# CMake generated Testfile for 
# Source directory: /home/sotnem/CLionProjects/biodynamo
# Build directory: /home/sotnem/CLionProjects/biodynamo/cmake-build-debug
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(runBiodynamoTests "/home/sotnem/CLionProjects/biodynamo/cmake-build-debug/runBiodynamoTests")
add_test(valgrind "valgrind" "--tool=memcheck" "--leak-check=full" "--show-leak-kinds=all" "--show-reachable=no" "--suppressions=/home/sotnem/CLionProjects/biodynamo/cmake-build-debug/../valgrind-biod.supp" "--error-exitcode=1" "./runBiodynamoTests")
