#!/bin/sh
CONDA_PREFIX=$(conda info --base)/envs/planarfold

g++ \
  -I$CONDA_PREFIX/include/python3.9 \
  -I$CONDA_PREFIX/lib/python3.9/site-packages/numpy/core/include \
  -I$CONDA_PREFIX/include/ \
  -L$CONDA_PREFIX/lib \
  -fpic -shared -Wl,-soname,dsflyext.so \
  -o dsflyext.so dsflyext.cpp libnetcdf.a planarfold_core.a \
  -lpython3.9  -lnetcdf -lboost_python39
