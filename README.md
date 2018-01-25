# py-NetHypGeom
Python wrapper for embedding networks in hyperbolic space and evaluating their navigability, with several performance improvements over the original package. (The original package, NetHypGeom, is written for R and available here https://github.com/galanisl/NetHypGeom.)

This script wraps certain methods available in NetHypGeom, provides convenience functions in Python, and offers an exhaustive but fast implementation for the evaluation of navigability, which can be used on any arbitrary distance matrix (the original function is slow and only does a random pairwise sampling of nodes).


Since it uses the R package NetHypGeom, R and the following R packages need to be installed:
- NetHypGeom (https://github.com/galanisl/NetHypGeom)
- igraph-R

Python 2.7 and the following packages are also required (not tested for other versions):
- rpy2
- numpy
- pandas
- Plotly (for plotting networks in hyperbolic space)
