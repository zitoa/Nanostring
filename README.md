## Background ##
microRNAs (miRNAs) are non-coding RNA molecules that play important roles in gene expression regulation. 
Several assays have been developed to study the miRNAs profiles. However, due to technical limitations 
and difficulties to accurately measure their expression levels, miRNAs functions are largely unexplored. 
Nanostring nCounter assay allows to assess the miRNAs expression levels.

## R code ##
Here a simple R code to preprocess miRNA expression profiles according to the nanostring's guidelines. 
Preprocessing includes background adjustment and subsequent, ofter required, normalization of the data.

Background adjustment is essential for accurate quantification of expression levels. Lane-specific background 
noise is quantified using the negative control beads. Negative controls are probes for which no target is expected 
to be present. The more negative probes are used the more accurate will be the estimation of background noise.
Background-noise adjustments can be done using the mean, the mean plus a multiple of standard deviation or the 
highest signal of the negative control probes. These methods differ each other with respect to the stringency 
of selection and FDR cutoffs. After quantification, background noise is substracted from the probe signals in 
a lane specific manner. For probes exhibiting counts near background levels and for which background adjustment 
results in negative or zero values, expression values can be set to arbitrary values (like 1). These probes could
be discarded or marked when fold-changes will be calculated.Background adjustment is followed by data normalization.

Normalization involves an initial step of lane-specific techinical standardization using positive controls. Then,
two types of normalization can be performed, which are Housekeeping (reference) genes or Global normalization. 
Housekeeping (reference) genes normalization relies on the assumption that housekeeping genes expression levels 
are similar between samples or replicates. Their target sequences are then assumed to have consistent expression 
levels. Global normalization relies on the assumption that most probes exhibit similar expression levels between 
samples. Only a minority of them will significantly deviate from the overall expression levels between samples. 
However, if bigger then expected fraction of probes show differential expression between samples this normalization 
type may not be appropriate. This method works best when comparing samples with similar overall miRNAs expression 
profiles. 
    
Background adjusted and normalized data can then be log2 transformed for downstream analysis. 

## Usage ##
The code is straightforward to use. You will need to download the mirnaNanostring.R script and load into
your working directory via source("path/to_the_script"). Afterwards, on the raw data as input and by 
specifying the metrics you want, you can implement the functions nanostring.bgAdjust() and 
nanostring.normalize() to respectively preprocess and then normalize the raw data. 
