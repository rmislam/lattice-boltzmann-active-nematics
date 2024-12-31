#!/bin/bash
echo "Compiling..."
nvcc -o prog main.cu -lm -lpthread -Xcompiler -march=native
echo "Running simulation..."
./prog
echo "Profiling..."
nvprof ./prog
echo "Generating visualizations..."
python visualization.py
echo "Done"