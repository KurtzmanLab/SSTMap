#!/bin/bash
gcc -O3 -lm -bundle -I /Users/kamranhaider/anaconda2/include/python2.7/ -I /Users/kamranhaider/anaconda2/lib/python2.7/site-packages/numpy/core/include/ _sstmap_ext.c -o _sstmap_ext.so -undefined dynamic_lookup
#gcc -O3 -lm -bundle -I /Users/kamranhaider/anaconda2/include/python2.7/ -I /Users/kamranhaider/anaconda2/lib/python2.7/site-packages/numpy/core/include/ _sstmap_ext_old.c -o _sstmap_ext_old.so -undefined dynamic_lookup
