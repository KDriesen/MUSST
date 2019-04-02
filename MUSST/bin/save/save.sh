#!/bin/bash

yea=$(date +%Y)
mon=$(date +%m)
day=$(date +%d)
hou=$(date +%H)
min=$(date +%M)
sec=$(date +%S)
fich=$yea$mon$day'_'$hou$min$sec

cd ..
MUSST/bin/7z/7z a -t7z -r -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on MUSST/$fich.7z                   \
            MUSST/src/*.f90 MUSST/makefile MUSST/src/inc_doc/*.md MUSST/css/* MUSST/MST.md \
            DIGIS/src/*.f90 DIGIS/makefile                                                 \
            MSOLV/src/*.f90 MSOLV/makefile MSOLV/inc/*.f90                                 \
            SPLIN/src/*.f90 SPLIN/makefile                                                 \
            UTILS/src/*.f90 UTILS/makefile                                                 \
            > MUSST/out/7z.log
cd MUSST

        


