# svd-protein-dynamics

Python script applying SVD (Singular Value Decomposition) to an MD (Molecular Dynamics) trajectory to determine the modes of motion of a protein. See "main.pdf" for my report and a more detailed description.

## Requirements
**MDAnalysis** (https://www.mdanalysis.org/) to read from and write to molecular dynamics file formats such as PSF, PDB, and DCD.

**NumPy** (https://numpy.org/) to manipulate matrices and calculate the SVD.

**Matplotlib** (https://matplotlib.org/) to make simple 2-D plots.

**VMD** (https://www.ks.uiuc.edu/Research/vmd/) to produce 3-D visualizations of the simulation
trajectory and its low rank approximation.

