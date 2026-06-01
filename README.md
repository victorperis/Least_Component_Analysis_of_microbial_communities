# Least Component Analysis (LCA) of microbial communities
This repository includes the code needed to produce all figure panels in the manuscript titled "Least Component Analysis reveals ecological constraints in microbial communities" (biorxiv DOI: https://doi.org/10.64898/2026.05.22.727221). 

LCA is just like Principal Component Analysis (PCA), but you consider the components of lowest variance, instead of the highest, while ensuring that they are statistically significant.

## How the repository is organized

Directories named "FIGX_description/" contain, as the name suggests, info regarding the panels in the corresponding figure in the manuscript. This might include the scripts that produce the panels or further information on how to obtain them. Each directory includes a README file with specific instructions.

The directory named "SI/" is organized in the same manner, but for figures in the supplementary material of the manuscript.

The "datasets" directory includes the different experimental and natural datasets analyzed in the manuscript. It includes the true data (obtained following instructions in the manuscripts we got them from), as well as null data obtained from one of two numerical null models (as described in the manuscript) and the names, based on taxonomy, for the taxa in each dataset.

The python notebook "full_data_analysis.ipynb" contains the basic code necessary to analyze a dataset, given the matrix of taxon abundances across samples. It takes the matrix of abundances, diagonalizes its correlation matrix, and gets statistically significant LCs and PCs. If no null data is provided, it performs basic shuffling. It also computes the similarity between LCs and calculates the coupling matrix between taxa. Results are output to a directory with the dataset name in "out/". 

The notebook "plot_eigenvector.ipynb" takes a component from an analyzed dataset and plots it. For the datasets used in the manuscript, it plots components in the same format as in the figure panels (unlabeled, dark-edged bars for experimental datasets; labeled, no-edged bars for natural datasets).

--- 

If you have any questions, feel free to contact me at: victor.peris@phys.ens.fr
