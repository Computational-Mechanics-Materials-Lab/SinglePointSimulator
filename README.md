# Single Point Simulator
### Computational Mechanics and Materials Laboratory @ Mississippi State University
#### For Questions: ch3136 AT msstate DOT edu, clarkhensley AT duck DOT com

## Description
A single-element finite element compatability layer in Python 3. This Single Point Simulator (SPS) uses the NumPy F2PY module to compile Fortran based (V)UMAT material models from Abaqus/Standard or Abaqus Explicit. Then, with a simple material property input file (in .toml or .json format), this tool quickly computes these single-element models and provides tools to visualize the results.
