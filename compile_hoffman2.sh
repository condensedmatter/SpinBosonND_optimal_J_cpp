#!/bin/bash
FFTW3_HOME="/u/local/apps/fftw3/current"

echo $FFTW_HOME
. /u/local/Modules/default/init/modules.sh

module load intel/13.cs
module load gcc/4.9.3

g++ -static -std=c++11 main_proj_1.cpp -o main_proj_1 -I$FFTW3_HOME/include -L$FFTW3_HOME/lib -lfftw3 -lm
