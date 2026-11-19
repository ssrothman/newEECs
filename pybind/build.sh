#!/bin/bash
# Builds the eecs_py Python module (pybind/eecs_py.cc and the library sources) into OUTDIR, with libfastjet beside it.
#   build.sh OUTDIR PYTHON_INCLUDE_DIR PYBIND11_INCLUDE_DIR EXT_SUFFIX
# gcc 11 and its fastjet, so that the module needs at most GLIBCXX_3.4.29, the libstdc++ of the coffea AlmaLinux 9
# container (the CMSSW_15 gcc 12 toolchain needs 3.4.30).
set -e
OUT=$1 PYINC=$2 PBINC=$3 SUFFIX=$4
[ -n "$SUFFIX" ] || { echo "usage: build.sh OUTDIR PYTHON_INCLUDE_DIR PYBIND11_INCLUDE_DIR EXT_SUFFIX" >&2; exit 2; }
HERE=$(cd "$(dirname "$0")/.." && pwd)
GCC=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gcc/11.4.1-30ebdc301ebd200f2ae0e3d880258e65
FJ=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fastjet/3.4.1-b5a7b930eb5755ed2b2b87b323687b41
BOOST=$(ls -d /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-*/include | head -1)
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/SRothman" "$OUT"
ln -s "$HERE" "$STAGE/SRothman/EECs"   # CMSSW include convention: SRothman/EECs/interface/X.h
export LD_LIBRARY_PATH=$GCC/lib64:$GCC/lib:$LD_LIBRARY_PATH   # cc1plus needs the isl, mpc, mpfr, gmp in $GCC/lib
"$GCC/bin/g++" -std=c++20 -O2 -fPIC -shared -fvisibility=hidden -Wall -Wextra -Werror=missing-braces \
    -I"$STAGE" -isystem "$FJ/include" -isystem "$BOOST" -isystem "$PYINC" -isystem "$PBINC" \
    "$HERE/pybind/eecs_py.cc" "$HERE/src/Binning.cc" "$HERE/src/Config.cc" \
    -L"$FJ/lib" -lfastjet -Wl,-rpath,'$ORIGIN' -o "$OUT/eecs_py$SUFFIX"
cp -L "$FJ/lib/libfastjet.so.0" "$OUT/"
echo "built $OUT/eecs_py$SUFFIX"
