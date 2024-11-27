# Warranties and Liability 

THE LATTICE MATCHING SCRIPT IS PROVIDED AS IS WITHOUT ANY WARRANTIES OR SUPPORT. YOU ARE SOLELY RESPONSIBLE FOR ITS USE.
THE CONTRIBUTORS HAVE NO OBLIGATIONS OR LIABILITY FOR ANY CONSEQUENCES FROM THE USE OF THIS SCRIPT OR ANY DERIVED WORKS.

# Citations

If you use this script or its concepts, please consider citing our article: 
Yizhen Song, Paulo C. D. Mendes, Sergey Kozlov. J. Mater. Chem. A, 2023, 11, 13665-13676 (https://doi.org/10.1039/d3ta01940c).

The function get_r_list() is inspired by ideas originally published by other authors. If you use the function, please cite their article:
A. Zur and T. C. McGill, J. Appl. Phys., 1984, 55, 378–386 (https://doi.org/10.1063/1.333084).

# Usage

Lattice matching is done by the function scan_composite_combinations(), which transforms and compares two matrices.
The set of candidates is controled by NSearch and tol. Selection is based on the user-defined deformation limit. 

The script has built in the example of ZnO/Cu with parameters that identify the 4x4/3x3 match. 

The output file results_explored_matrices.txt contains the list of combinations found from the scan. 
If you find no matches for your system, you can consider adjusting the parameters in the script.

# Contributors

Programming, conceptualization, and testing: Paulo C. D. Mendes
Revisions and testing: Yizhen Song
Conceptualization and supervision: Sergey Kozlov

Published upon request by Prof. Sergey Kozlov.
