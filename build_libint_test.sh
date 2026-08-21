#!/bin/bash
set -e
cd /home/jieyang/apps/scfcxx
INC="-Isrc -Ithird_party/eigen $(pkg-config --cflags libint2)"
LIBS="$(pkg-config --libs libint2)"
B=/tmp/scfbuild
mkdir -p $B
g++ -std=c++17 -O1 -Wall -Wextra $INC -DSCFCXX_BASIS_PATH=\"$PWD/third_party/libint2_basis\" \
  -c src/Libint2IntegralProvider.cpp -o $B/Libint2IntegralProvider.o
g++ -std=c++17 -O1 -Wall -Wextra $INC -c src/basis.cpp -o $B/basis.o
g++ -std=c++17 -O1 $INC -DSCFCXX_BASIS_PATH=\"$PWD/third_party/libint2_basis\" \
  tests/test_libint.cpp src/SzaboHeHIntegral.cpp $B/Libint2IntegralProvider.o $B/basis.o $LIBS -o $B/test_libint
echo '=== running test_libint ==='
$B/test_libint
