#!/bin/bash

rm bin/plambda; make bin/plambda CFLAGS='-g -O2 -DNDEBUG -DDONT_USE_TEST_MAIN'
plambda zero:512x512 "x 255 +" -o /tmp/k.png  # this runs smoothly

rm bin/plambda; make bin/plambda CFLAGS='-g -O2 -DNDEBUG -DDONT_USE_TEST_MAIN -fopenmp'
plambda zero:512x512 "x 255 +" -o /tmp/k.png  # this segfaults 
