#!/bin/bash
echo "Compiling..."
gcc -o prog main.c -lm -lpthread -march=native -std=gnu99 -fopenmp -Wall
echo "Running simulation..."
./prog
echo "Generating visualizations..."
python visualization.py
echo "Done"