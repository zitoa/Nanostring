###############################################
##### R code for background adjustment and ####
### normalization of miRNA expression data ####
##### based on the Nanostring's guidelines ####
###############################################

nanostring.bgAdjust <- function(rawData, bg.method=c("mean.neg.ctrl", "mean.sd.neg.ctrl", "highest.neg.ctrl")) {

    if(!class(rawData) %in% c("matrix","data.frame"))
        stop("Input data must be a matrix or dataframe")
    
    if(!bg.method %in% c("mean.neg.ctrl","mean.sd.neg.ctrl","highest.neg.ctrl"))
        stop("Background adjustment must be one of the following: mean.neg.ctrl, mean.sd.neg.ctrl, highest.neg.ctrl")

    ## rawData is a matrix/dataframe containing raw unprocessed miRNA expression data. ##
    ## rows specifies the probes and columns the samples. ##
    ## Negative controls (probes for which no target is expected to be present) are ##
    ## here used to quantify the lane's background noise. The more negative probes ###
    ## are used the more accurate will be the estimation of systemic background noise ##
    ## Several methods exist to quantify the background noise. They differ with  #######
    ## respect to the stringency of selection and FDR cutoffs. Background subtraction ##
    ## will amplify estimations of fold-change in differential expression analysis #####
    ## For probes exhibiting counts near background levels, subtraction may result in ##
    ## negative or zero values being obtained. If fold change determination for such a ##
    ## target is desired, the value for this target should be set equal to the background ##
    ## level in that lane (or some other arbitrary value, such as 1). Fold change ##########
    ## determinations of genes at or near background levels should be considered an estimate.
    ## neg in a vector of integers with the position of negative controls in the matrix. ##
    Data <- apply(rawData, 2, as.numeric)  
    rownames(Data) <- rownames(rawData)

    if(bg.method == "mean.neg.ctrl")  {
        bg.cutoffs <- apply(Data[neg,], 2, mean)
    }
    
    if(bg.method == "mean.sd.neg.ctrl")  {
        i=1;
        bg.cutoffs <- c();
        for(i in 1:ncol(Data))  {
            bg.cutoffs <- c(bg.cutoffs, mean(Data[neg,i]) + 2*sd(Data[neg,i]))
            names(bg.cutoffs)[i] <- colnames(Data)[i]
        }
    }

    if(bg.method == "highest.neg.ctrl")  {
        i=1;
        bg.cutoffs <- c();
        for(i in 1:ncol(Data))  {
            jj <- which.max(Data[neg,i]);
            bg.cutoffs <- c(bg.cutoffs, Data[neg,i][jj]);
        }
    }
            
    ## matrix with data adjusted for the lane-specific background thresholds ##
    Data.bg.Adjusted <- matrix(nrow=nrow(Data),ncol=ncol(Data))
    rownames(Data.bg.Adjusted) <- rownames(Data)
    colnames(Data.bg.Adjusted) <- colnames(Data)
    i=1;
    for(i in 1:length(bg.cutoffs))  {
         Data.bg.Adjusted[,i] <- Data[,i] - bg.cutoffs[i]
     }
    
    return(Data.bg.Adjusted)
    
} 


nanostring.normalize <- function(Data.bg.Adjusted, norm.method=c("Housekeeping.ctrl", "global.ctrl"),
                                 method.global=c("total.count","average.count","geom.mean"))  {     

    if(!norm.method %in% c("Housekeeping.ctrl","global.ctrl"))
        stop("Normalization method must be one of the following: Housekeeping.ctrl, global.ctrl")

    if(norm.method == "global.ctrl") {
        if(!method.global %in% c("total.count","average.count","geom.mean"))
            stop("When normalization is global, the metrics must be one of total.count, average.count, geom.mean")
    }

    ## Lane specific technical standardization using positive controls: ##
    ## Here, we adjust for platform-associated sources of variation  ##
    ## Variation due to technical/biological replicates is not taken into account. ##
    ## Lane-specific positive control scaling factors (sc) are calculated. ##
    ## Values of sc outside the range of 0.3-3, may indicate under-performance of a lane. ##
    ## Below, pos is a vector of integers corresponding to the ##
    ## position of positive controls in the matrix Data. ##
    Data <- data.bg.subtracted
    Data1 <- apply(Data[pos,], 2, mean)
    global.mean <- mean(Data1) 
    
    i=1;
    scaling_factors <- c(); 
    for(i in 1:length(Data1)) {
        scaling_factors <- c(scaling_factors, global.mean/Data1[i])
    }

    jj <- which(scaling_factors<0.3 | scaling_factors>3)
    if(length(jj) !=0) 
        message("Under performed lanes (samples) detected in the platform")
    

    ## Technical standardized matrix ##
    Data2 <- matrix(nrow=nrow(Data),ncol=ncol(Data))
    rownames(Data2) <- rownames(Data)
    colnames(Data2) <- colnames(Data)
    i=1;
    for(i in 1:ncol(Data2))    {
        Data2[,i] <-  Data[,i] * scaling_factors[i]
    }
    
    geo.mean <- function(x) {
        prodx <- prod(x)
        n <- length(x)
        gm <- prodx ^ (1/n) ## or prodx**(1/n) 
    }

    ## Normalization: #################################################################
    ## After the technical standardization, further normalizations may be required. ###
    ## The best normalization strategy will ultimately depend on experimental setup. ## 
    ############### #################### ################# ################# ##########
    ## Housekeeping (reference) genes normalization relies on the assumption that #####
    ## housekeeping genes expression levels are similar between samples or replicates #
    ## Their target sequences are then assumed to have consistent expression levels. ##
    ## The more reference genes are used the more accurate the normalization will be. #
    ## house and endo are vectors of integers with the positions of housekeeping genes#
    ## and endogenous genes in the matrix, respectively. ##############################
    ############### #################### ################# ################# ##########
    ## Global normalization relies on the assumption that most probes exhibit similar #
    ## expression levels between samples. Only a minority of them will significantly ##
    ## deviate from the overall expression levels between samples. However, if bigger #
    ## then expected fraction of probes show differential expression between samples ##
    ## this normalization type may not be appropriate. This method works best when ####
    ## comparing samples with similar overall miRNAs expression profiles. #############
    ############### #################### ################# ################# ##########
    if(norm.method == "Housekeeping.ctrl") {  

        data.norm <- apply(Data2[house,], 2, geo.mean)     
        global.geomean <- mean(data.norm) 
        
        normalization_factors <- c()     
        i=1;
        for(i in 1:length(data.norm))   {
            normalization_factors <- c(normalization_factors, global.geomean/data.norm[i])
        }

        ## Normalized matrix ##
        res <- matrix(nrow=nrow(Data2[endo,]),ncol=ncol(Data2))
        rownames(res) <- rownames(Data2[endo,])
        colnames(res) <- colnames(Data2)
        i=1;
        for(i in 1:ncol(Data2))    {
            res[,i] <- Data2[endo,i] * normalization_factors[i]
        }

    }
        
    if(norm.method == "global.ctrl")  {

        if(method.global == "total.count") 
            data.norm <- apply(Data2[endo,], 2, sum)

        if(method.global == "average.count") 
            data.norm <- apply(Data2[endo,], 2, mean)

        if(method.global == "geo.mean") {
            jj <- apply(Data2, 1, mean)
            jj1 <- jj[order(jj,decreasing=TRUE)]
            data.norm0 <- Data2[names(jj1),]
            data.norm <- apply(Data.norm0[1:100,], 2, geo.mean)     
        }
            
        global.metric <- mean(data.norm) 
        normalization_factors <- c()
        i=1;
        for(i in 1:length(data.norm) )   {
            normalization_factors <- c(normalization_factors, global.metric/data.norm[i])
        }
        
        ## Normalized matrix ##
        res <- matrix(nrow=nrow(Data2),ncol=ncol(Data2))
        rownames(res) <- rownames(Data2)
        colnames(res) <- colnames(Data2)
        i=1;
        for(i in 1:ncol(Data2))    {
            res[,i] <- Data2[,i] * normalization_factors[i]
        }
        
    }

    return(res)
              
}   
 



