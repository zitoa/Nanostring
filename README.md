## Background ##
microRNAs (miRNAs) are non-coding RNA molecules that play important roles in gene expression regulation. 
Several assays have been developed to study the miRNAs profiles. However, due to technical limitations 
and difficulties to accurately measure their expression levels, miRNAs functions are largely unexplored. 
Nanostring nCounter assay allows to assess the miRNAs expression levels. 

## R code ##
Here a simple R code to preprocess miRNA expression profiles according to the nanostring's guidelines. 
Preprocessing includes background adjustment and subsequent, ofter required, normalization of the data.

Background adjustment is essential for accurate quantification of expression data. It is performed using 
the platform's negative control beads.
Background-noise adjustments could be done using the average of negative control signals, the average and the 
standard deviation of negative control signals, or the highest signal between the negative controls. 
Background subtraction is followed by data normalization. Finally, the expression values can be 
log2 transformed for downstream analysis. 

