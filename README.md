microRNAs (miRNAs) are non-coding RNA molecules that play important roles in gene expression regulation. 
Several assays have to date been developed to study the miRNAs profiles but due to technical and difficulty to accurately detect them, are largely under explored.

miRNAs high-throughput data still lacks a well-suited pipeline describing the best methodology to be applied depending on the experimental
design

MiRNA expression profile data were treated as recommended by the Nanostring’s
guidelines. Data pre-processing consists mainly of background adjustment steps which
can be carried out using the negative control beads provided by technology. Specif-
ically, it is possible to perfom the background correction using the mean of negative
controls, the mean plus standard deviation of negative controls, or the highest value
of negative controls. Subsequently to background substraction, it is recommended to
normalize the data. A detailed description of the implemented methods is provided by
R code in the appendix. After having perfomed these transformation, the data were
log2 transformed. 

Here a simple R script to pre-process and normalize nanostring miRNA expression profiles according to the nanostring's 
guidelines (2014).

