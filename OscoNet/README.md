# OscoNet
Bootstrap-based OscoNet method: 
Method to infer sinusoidal oscillations in single cell data.

## Installation
1. Create new environment `conda create --name osconet python=3.8` and to activate it with 
`conda activate osconet` before proceeding.
2. Install required packages using `pip install -r requirements.txt`
3. Install package `pip install -e .`
4. Install numba `conda install numba`
5. Verify your installation by running `pytest` from the project root directory `osconet`.
Note this can take around 1-2 minutes.


## Synthetic data (Optional)
This is a simple run of the core Oscope script on synthetic data as a test. Not essential to do this. This will not output the Results and Plots of a full OscoNet run.

1. Run ` python OscopeBootstrap/oscope_tf.py --test_mode
` for a simple demonstration of the method on synthetic data. This will run under a quick configuration
to demonstrate the capabilities of the method. This should take 10-20 seconds.
2. Remove the `--test_mode` flag for a 1000-sample bootstrap test on the exact synthetic run configuraiton used in the paper
(1000 genes, 100 cells with 3 clusters of co-oscillating genes).

## Notebooks (Optional)
We have compiled a series of notebooks which can be worked through to provide more undestanding of the method. Some notebooks can be used to test reproduction of our data and ensure functionality of your installation.

1. `notebooks/OscoNet introduction.ipynb`: provides an introduction to the hypothesis test on a simple synthetic example.
2. `notebooks/Reproduce_figures_5_7.ipynb` : pseudotime on Whitfield microarray data. To see how the spectral embedding
pseudotime method can be applied.

## Data 
Ensure you have your scRNA-seq gene expression matrix downloaded into /OscoNet/YourDataset/Data. Now is a good time to rename YourDataset to something referring to your specific experiment. The data can be provided as .csv, .tsv or .txt, however you will need to uncomment the correct line (42-44) in /OscoNet/YourDataset/Code/1-Recluster.R.

## Runinng Osconet
The OscoNet pipeline consists of a series of scripts which should be ran in a specified sequence. 

### 1. Reclustering

The 1-Recluster.R script can be performed regardless of whether you have pre-clustered your data. OscoNet looks for repeating relationships in expression between pairs of genes over time to infer oscillations. If your data includes heterogenous cells then it may misconstrue inter-cell type differences in gene expression as periodic fluctuations. To this end, input data should be as homogenous as possible. This script will use Seurat to identify clusters in your dataset so OscoNet can be run separately on each cluster. 

If you have pre-clustered your data and are happy with each cluster's homogeneity, then you can skip this step. You will have to provide a metafile named clusters.csv with the list of cell IDs in the 1st column and the corresponding cluster ID for that cell in a 2nd column, this file should include a header. These cell IDs should correspond to those in your expression matrix file, which should be named counts.csv. Examples of each of these can be found in the 'examples' directory.

There is an option to run this reclustering on a subset of your data. Say you have pre-determined distinct cell types with different IDs then you can choose to run this on only cells matching a certain criteria such as those named "CD44+" cells. To do so, just uncomment the relevant sections of the code.

Remember to modify the code to reflect the pathname of your local OscoNet directory.

This script will output an expression data file named counts and a metafile named clusters in the Data sub-directory.

### 2. Data filtering

The 2-NormFilter.R script should be run next. This script performs two filtering steps to prepare for the OscoNet run. 

The first is the alpha filter, which removes genes which do not meet a defined threshold for percentage of cells it is expressed in. The threshold can be modified by changing the alpha parameter. For example, an alpha value of 0.2 (recommended) will filter out genes which are not expressed in >= 20% of cells (or expressed in at least 80% of cells).

The second is the variance filter which will remove genes which do not exhibit variant expression across the dataset, as these are unlikley to be oscillators.

Remember to modify the code to reflect the pathname of your local OscoNet directory. If you have run the code on a subpopulation of the data, then make sure filenames are appropriately named and specified. Also edit the 'case' value to reflect the dataset or subpopulation.

This will return filtered files to the OsconetInput folder ready to run OscoNet in the next step.

### 3. Oscope Bootstrap

The 3-Run_All.sh script is run next. This is a shell script which will run collaboration.py on each individual cluster. This will do the bootstrapping of and calculate the distance between each gene pair. 

First edit the experiment variable on line 6 to match the 'case' value set in the previous script. Next open the terminal and activate your conda environment with `conda activate osconet`. Then execute the shell script.

This will create a directory named OsconetOutput with an out.csv and psi.csv for each cluster. This is now ready for the community extraction.

### 4. Community Extraction

This final step should be performed with 4-CommunityExtract.R. This script runs diagnosis and clustering on the graphs produced as a result of OscoNet. For each cluster the number of co-oscillating genes as well as the total number of genes which passed filtering will be outputted. Communities are extracted using the Walktrap algorithm from the igraph package. A plot of the Network with overlaid communities will be created in the Plots sub-directory. 

OscoNet cannot distinguished between pairs of genes which co-oscillate and those which are simply co-expressed (linear). To handle this, a linearity test is run on each individual community. Communities of potentially linearly co-expressed genes are flagged with a 1 (default value non-linear communities is 0). 

Remember to modify the code to reflect the pathname of your local OscoNet directory. If you have run the code on a subpopulation of the data, then make sure filenames are appropriately named and specified. Also edit the 'experiment' value to reflect the dataset or subpopulation. Modify the number of 'case' list variable to reflect the number of clusters identified in your data. Please also ensure the alpha variable mathches that performed in the filtering step.



