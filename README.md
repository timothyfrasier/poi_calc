# poi_calc
This repository contains an R script that can implement several functions for calculating the probability of identity. Notes describing these functions are included as comments within the code.

It has several dependencies on the [gstudio](https://github.com/dyerlab/gstudio) R package, and therefore users must have that installed and loaded into R as well (prior to running these scripts). A walk-through of the whole process is given below, using the example data set.

Make sure that your genotype file, and the **poiCalc.R** script are in R's working directory.

Also not that the first columns of your genotype file should have the header "ID", and should contain the names of your samples/individuals.

First, load the **gstudio** package, and source the **poiCalc.R** script.    
`library(gstudio)`    
`source("poiCalc.R")`    

Read the data into **gstudio**. Note that this genotype file contains genotypes of 50 individuals genotyped at 18 microsatellite loci. The sample number is in column 1, and the genotypes are in columns 2-37. Note that missing data should first be indicated as NA.    
`data <- read_population("genotypes.csv", type = "column", locus.columns = 2:37, sep = ",", header = TRUE)`    

Get the allele frequency data.    
`freqs <-  frequencies(data)`    

Calculate the poi for the entire population using these allele frequency data.    
`poi.all(freqs)`    

Calculate the poi_sib for the entire population using these allele frequency data.    
`poi.all.sib(freqs)`    

Calculate poi for each individual genotype.    
`out <- poi.ind(microdata = data, freqs = freqs)`    
`out`    

Calculate the poi_sib for each individual genotype.    
`out_sib <- poi.ind.sib(microdata = data, freqs = freqs)`    
`out_sib`    
