#!/bin/bash
echo "Compiling..."
nvcc -o prog main.cu -lm -lpthread -Xcompiler -march=native
echo "Running simulation..."
./prog
#echo "Profiling..."
#ncu -o profile -f --launch-skip 100 --launch-count 32 --section SchedulerStats --section WarpStateStats ./prog
echo "Generating visualizations..."
python visualization.py
echo "Done"