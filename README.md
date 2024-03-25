# Single Point Simulator
### Computational Mechanics and Materials Laboratory @ Mississippi State University
#### For Questions: ch3136 AT msstate DOT edu, clarkhensley AT duck DOT com

## Description
A single-element finite element compatability layer in Python 3. This Single Point Simulator (SPS) uses the NumPy F2PY module to compile Fortran based (V)UMAT material models from Abaqus/Standard or Abaqus Explicit. Then, with a simple material property input file (in .toml or .json format), this tool quickly computes these single-element models and provides tools to visualize the results.

## Installation
```sh
python -m pip install singlepointsim
```

## Basic Usage
```sh
python -m singlepointsim <inputfile>
```

See example\_input.toml for an overview of how to format the input file.
(.json input format is also accepted by the singlepointsim)

## Recommended Usage
Consider a directory containing the desired (V)UMAT fortran file, `foo.f90`
In the same directory, create `input.toml` using the same format as `example_input.toml` as shown in this github repository.
Then, simply running:
```sh
python -m singlepointsim input.toml
```
should compile the (V)UMAT to a .so or .dll if it is not already compiled and should call on this object to run analysis.
Results will be stored in a .csv file in the desired "results" directory, as given by `input.toml` (or the present working directory if not specified).
