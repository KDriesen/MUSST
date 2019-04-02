#!/bin/bash

#~ PROFILAGE : après avoir lancé le prog, gmon.out est créé, lancer : gprof -b prg gmon.out > gprof.txt
python gprof2dot.py gprof.txt | dot -Tpng -o output.png
