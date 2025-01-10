# scTEL reproducibility

This repository contains the code to reproduce the analyses performed in the manuscript "scTEL: A Joint Analysis Framework of CITE-seq and scRNA-seq Using Transformer".
This repository contains the code to reproduce the analyses performed in the manuscript "A joint analysis of single cell transcriptomics and proteomics using transformer".

<img src="/figures/scTEL.png" alt="scTEL" style="zoom: 10%;" />

# General Flow

It is recommended the user proceeds as follows.

1. Clone this repository to their local machine.
2. Download the data.
3. Install necessary packages.
4. Run scTEL notebooks.

## Directory

- `Data` stores the datasets
- `scTEL` contains the main script of scTEL
- `test` contains the testing script of scTEL on the datasets in the manuscript and running script of baseline methods. 

## Data source and reference

(1) Seurat 4 human peripheral blood mononuclear cells (GEO: GSE164378). (Hao Y et al.(2021))

(2) Mucosa-Associated Lymphoid Tissue (MALT) dataset. ( Justin Lakki et al.(2021))

(3) H1N1 influenza PBMC dataset. (Kotliarov et al. (2020))

(4)  Human blood monocyte and dendritic cell CITE-seq dataset (Monocytes
dataset). ( Justin Lakki et al.(2021))

The University of Pennsylvania has put these data sets together for the convenience of downloading. Download [here](https://upenn.app.box.com/s/1p1f1gblge3rqgk97ztr4daagt4fsue5) and place them into the currently empty data folder.

## Install necessary packages

Recomended installation procedure is as follows.

1. Install [Anaconda](https://www.anaconda.com/products/individual) if you do not already have it, so that you can access conda commands in terminal.

2. Create a conda environment, and then activate it as follows in terminal. The code here requires the versions of the packages specified in the `scTEL_env.yml` file. 

   ```
   $ conda env create -f scTEL_env.yml
   $ python -m ipykernel install --user --name=scTEL
   ```

### Run scTEL notebooks

It is recommended that the user first run the scTEL notebooks. Simply, open each of the following notebooks in jupyter. Make sure to set the active conda kernel in jupyter to "scTEL" and then run all cells. Repeat this for every notebook listed below.

1. [Monocyte 1Monocyte 2.ipynb](https://github.com/)
2. [PBMC_H1N1.ipynb](https://github.com)
3. [pbmc_Malt.ipynb](https://github.com/)
4. [PBMC_PBMC.ipynb](https://github.com)
