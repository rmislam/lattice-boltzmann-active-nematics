#!/bin/bash
echo "Compiling..."
gcc -o prog main.c -lm -lpthread -march=native -std=gnu99 -fopenmp -Wall
echo "Running simulation..."
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=test.prof ./prog
echo "Generating visualizations..."
python visualization.py
echo "Done"