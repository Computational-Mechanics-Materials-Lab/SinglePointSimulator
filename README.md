# Single Point Simulator
### Computational Mechanics and Materials Laboratory @ Mississippi State University
#### For Questions: ch3136 AT msstate DOT edu, clarkhensley AT duck DOT com

## Description
An Abaqus FEA single-element finite element compatability layer in Python 3. This Single Point Simulator (SPS) uses the NumPy f2py module to compile Fortran based UMAT material models from Abaqus/Standard or Abaqus Explicit. Then, with a simple material property input file (in .toml or .json format), this tool quickly computes these single-element models and provides tools to visualize the results.

## Installation
```sh
python -m pip install singlepointsim
```

## Basic Usage
```sh
# the umat can be specified in the input file or passed as a command-line argument
python -m singlepointsim <inputfile> (<umat>)
```

See example\_input.toml for an overview of how to format the input file.
(.json input format is also accepted by the singlepointsim)

## Recommended Usage
Consider a directory containing the desired UMAT fortran file, `foo.f`
In the same directory, create `input.toml` using the same format as `example_input.toml` as shown in this github repository.
Then, simply running:
```sh
python -m singlepointsim input.toml foo.f
```
should compile the UMAT to a .so or .dll if it is not already compiled and should call on this object to run analysis.
Results will be stored in a .csv file in the desired "results" directory, as given by `input.toml` (or the present working directory if not specified).

## Issues?
Fortran, especially Fixed-Form/F77 Fortran, can be finnicky. Please don't hesitate to reach out for support on GitHub, but, ensure that the meson patch built into singlepointsim has successfully run (may require reinstalling meson) and that your file is using either the '.f' extension (for fixed-format Fortran) or the '.f90' extension (for free-format Fortran).
