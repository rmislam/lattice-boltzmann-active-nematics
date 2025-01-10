# Lattice Boltzmann Active Nematics

C implementation of 2D Beris-Edwards nematodynamics with activity patterns, as conducted in Zhang et. al.'s paper "[Logic operations with active topological defects](https://www.science.org/doi/10.1126/sciadv.abg9060)".

This code is an extension of the wonderful [nematodynamics tutorial](https://zenodo.org/records/4737814) provided by Ziga Kos and Miha Ravnik.

## Environment
Using a Python virtual environment is recommending for managing dependencies for plotting. After creating a `venv` environment (simply called `venv` below), run these commands to install the Python dependencies:

```
source venv/bin/activate
pip install -r requirements.txt
```

## How to run
This command will compile the simulation, run it, and produce two GIFs (one for the velocity field, one for the director field) visualizing the simulation output.

```
./run_all.sh
```

## CUDA
The `cuda` branch provides a CUDA implementation compiled with `nvcc` that is much faster than the CPU-only OpenMP `master` branch. The same `./run_all.sh` script is used to run it.
