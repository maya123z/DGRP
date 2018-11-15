# Clean and Process GWAS Results for the Drosophila Genetic Reference Panel

This package takes the raw results output after performing a genome-wide association study (GWAS) using the website for the Drosophila Genetic Reference Panel (DGRP, http://dgrp2.gnets.ncsu.edu/). It organizes the results into a more easily readable format. It can also accept a binary phenotype file and use this to perform additionary analyses such as odds ratios.

References:
[Mackay et al. 2012](https://www.nature.com/articles/nature10811)

## General Usage Guide 

The DGRP-Cleanup package was created so that Drosophila scientists without a programming background can still make use of the powerful GWAS tools provided by the DGRP. While the website performs the primary analysis for you, it can still be difficult to parse the results that it outputs. This general usage guide serves as a foolproof guide for scientists to use the DGRP package regardless of their experience level with R.

NOTE: This is my very first R package, so if you have any suggestions for improvement, please don't hesitate to contact me using the Issues tab!

### Step 1: Install R

This program runs on R. I'd recommend also installing RStudio, which provides a nicer user interface. Click here for instructions on how to install both these programs on your computer.

If you're completely new to R, I'd recommend using free courses such as DataCamp's Introduction to R to get a handle on the basics. You definitely don't need to be an R expert to use this package, but it helps to have a basic understanding of syntax and usage.

### Step 2: Install the DGRP Package

First, open the RStudio program. In the Console, type the following:

`install.packages("devtools")`

This will install another package that you need to access GitHub. Once the package is installed you need to load it:

`library(devtools)`

Finally, use the following command to install the DGRP package:

`install_github("maya123z/DGRP-Cleanup")`

To begin using the package, just load it with this command:

`library(DGRP)`

### Step 3: Load Your Data

The DGRP tool allows you to run GWAS by simply uploading your phenotypic data to their website. A few minutes later, they'll send you an email with a link to your results.

[in progress--more coming soon!]
