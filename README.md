# Guide
## Data
Raw and processed data are available at GSE226502 (reviewer token required)
## Genomes and annotations
Data from the 129 and NOD mouse strains were aligned to custom reference genomes. These genomes are located in the `fasta` directory.

The *Ly49* and *MAP1/8* annotations used in this study are located in the `annotations` directory.
## Code
### non-R code
The [snakeATAC](https://github.com/ChangxuFan/snakeATAC) pipeline, written by me, was used for ATAC-seq processing where indicated. 
These are mostly scenarios where the functions of the [AIAP](https://github.com/Zhang-lab/ATAC-seq_QC_analysis) pipeline did not suffice.

The [WGBS](https://github.com/ChangxuFan/wgbs/tree/fanc) pipeline was used for bisulfite data processing. 
It was originally written by [Hyung Joo Lee](https://github.com/hyungjoo-lee/wgbs). 
Xiaoyu Zhuo and I made it into a docker+snakemake pipeline.
### R code
Most analyses code were written as functions and organized into R packages. 
The `analyses` directory of this repo contains R scripts that call these functions. 
Therefore, to recapitulate the analyses, the following R packages (written by me, hence the suffix "Fanc") need to be installed. 
Please note that these packages do not load any dependencies. If you have trouble installing them as packages, 
you could simply source the scripts under the `R/` directory.
* [abaFanc](https://github.com/ChangxuFan/abaFanc)
* [abaFanc2](https://github.com/ChangxuFan/abaFanc2)
* [bamFanc](https://github.com/ChangxuFan/bamFanc)
* [cageFanc](https://github.com/ChangxuFan/cageFanc)
* [common](https://github.com/ChangxuFan/common)
* [liteRnaSeqFanc](https://github.com/ChangxuFan/liteRnaSeqFanc)
* [scFanc](https://github.com/ChangxuFan/scFanc)
* [v4c](https://github.com/ChangxuFan/v4c)
* [utilsFanc](https://github.com/ChangxuFan/utilsFanc)

Some bash scritps and R scripts were also used. They are collected at:
* [scripts](https://github.com/ChangxuFan/scripts)
* [R_for_bash](https://github.com/ChangxuFan/R_for_bash)

Please note that the goal of this study is not to make R packages. 
These code are released to show how I performed each analyses exactly, as an effort to facilitate open science. 
Therefore, you might need to adjust the code for it to run properly on your computer. 
Most notably, the R code calls bash tools according to where they are installed on our server. 
It is likely that you would have to modify these paths according to where these tools are installed on your computer. 
These code are not likely to run on a Windows machine due to function calls such as `system("ln -s `realpath source` destination")`
