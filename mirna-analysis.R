#############################################
###### Code to preprocess  miRNA data #######
########## Nanostring guidelines ############
#############################################


nanostring.normalize <- function( data.bg.subtracted, norm.method=c("Housekeeping.ctrl",
                                                                    "global.ctrl") )  

{     

### Technical standardization using positive control
### In this step, we adjust for all platform associated source of variation
### It will not account for variation due to technical-biological replicates  
### If the calculated scaling factor is outside a range of 0.3-3, it may be
### under-performance of lane or lanes.
    #data <- apply( raw0, 2, as.numeric )  # ? already did in backgr. function??
    #rownames(data) <- rownames(raw0)      # ? already did in backgr. function??
    data <- data.bg.subtracted
    data1 <- apply( data[pos,], 2, mean )
    global.mean <- mean(data1) 

    scaling.factors <- c()    ## line specific
    i=1;
    for( i in 1:length(data1) )   {
        scaling.factors <- c( scaling.factors,global.mean/data1[i] )
    }

    if( all.equal(names(scaling.factors),colnames(data)) != TRUE ) { stop("error1") } 
    if( all.equal(length(scaling.factors),dim(data)[2]) != TRUE ) { stop("error1.1") } 

    f <- which( scaling.factors < 0.3 )  
    f1 <-  which( scaling.factors > 3 )   
    if( length(f) | length(f1) !=0 )    { message("Exist under performed samples") } 
                                        #ff <- names(scaling.factors[c(f,f1)]
                                        #cat( "under performed lines:", ff, sep=" ") 

    data2 <- matrix( nrow=nrow(data), ncol=ncol(data) )
    rownames(data2) <- rownames(data)
    colnames(data2) <- colnames(data)

     if( all.equal(names(scaling.factors),colnames(data2))==TRUE &&
       all.equal(length(scaling.factors),dim(data2)[2])==TRUE )   { 

        i=1;
        for(i in 1:ncol(data2))    {
            s <- data[,i]
            ss <- s*scaling.factors[i]
            data2[,i] <- ss
        }
    } else { stop("error2") }
            

             if(norm.method=="Housekeeping.ctrl")    {  
### Housekeeping genes are not expected to vary between samples or replicates.
### This normalization assume that the Housekeeping target sequences are consistent
### in their expression levels.
                 geo.mean <- function(x){
                                        # geometric mean calculation
                     prodx <- prod(x)
                     n <- length(x)
                     gm <- prodx ^ (1/n)    ## or prodx**(1/n)
                     gm
                 }
                 
                 data.norm <- apply( data2[house,], 2, geo.mean )     
                 global.geomean <- mean(data.norm) 

                 normalization.factors <- c()     ## line specific
                 i=1;
                 for( i in 1:length(data.norm) )   {
                     normalization.factors <- c( normalization.factors,global.geomean/data.norm[i] )
                 }

                 TT <- matrix( nrow=nrow(data2[endo,]), ncol=ncol(data2) )
                 rownames(TT) <- rownames(data2[endo,])
                 colnames(TT) <- colnames(data2)

                 if( all.equal(names(normalization.factors),colnames(data2))==TRUE &&
                     all.equal(length(normalization.factors),dim(data2)[2])==TRUE )   { 

                     i=1;
                     for(i in 1:ncol(data2))    {
                         s <- data2[endo,i]
                         ss <- s*normalization.factors[i]
                         TT[,i] <- ss
                     } 
                 }   else { stop("error3") }

                 return(TT)

             }   ## end if cycle

               if(norm.method=="global.ctrl")  {
### Also if some targets may increase or decrease in any given sample, these are likely
### to represent a small portion of the total, and so that the overall level of expression
### within a sample will be the same. In this cases a global normalization it is useful.
### But, If a significant fraction of probes exhibit differential expression from sample
### to sample, this method may not be appropriate

                 data.norm <- apply( data2[endo,], 2, sum )
                 global.sum <- mean(data.norm) 

                 normalization.factors <- c()     ## line specific
                 i=1;
                 for( i in 1:length(data.norm) )   {
                     normalization.factors <- c( normalization.factors,global.sum/data.norm[i] )
                 }

                 TT <- matrix( nrow=nrow(data2), ncol=ncol(data2) )
                 rownames(TT) <- rownames(data2)
                 colnames(TT) <- colnames(data2)

                 if( all.equal(names(normalization.factors),colnames(data2))==TRUE &&
                     all.equal(length(normalization.factors),dim(data2)[2])==TRUE )   { 

                     i=1;
                     for(i in 1:ncol(data2))    {
                         s <- data2[,i]
                         ss <- s*normalization.factors[i]
                         TT[,i] <- ss
                     }

                  }  else { stop("error3") }    

                 return(TT)
              
              }   ## end if cycle

   

}   ## end function      



nanostring.bg.adj<- function( raw0, bg.method=c("mean.neg.ctrl", "mean.sd.neg.ctrl",
                                                "highest.neg.ctrl") )
                                                
{

### In the dataset, there are several probes for which don't exist the relative targets
### The negative control can be used to estimate systematic background counts.This estimation
### can be performed in many ways, each of which change from the other in terms of
### stringency and FDR, and so FNR.
    data <- apply( raw0, 2, as.numeric )  
    rownames(data) <- rownames(raw0)

    if(bg.method=="mean.neg.ctrl")  {
        bg.thresholds <- apply( data[neg,], 2, mean )     ## line specific
    }

    
    if(bg.method=="mean.sd.neg.ctrl")  {
        bg.mean <- apply( data[neg,], 2, mean )     ## line specific
        bg.sd <- apply( data[neg,], 2, sd )         ## line specific
        bg.thresholds <- c()
        i=1;
        for(i in 1:dim(data)[2])  {
             bg.thresholds <- c( bg.thresholds, bg.mean[i]+bg.sd[i] )  ## line specific
         }
    }
            
        
    if(bg.method=="highest.neg.ctrl")  {
        j <- apply( data[neg,], 2, which.max )    
        bg.thresholds <- c()
        i=1;
        for(i in 1:dim(data)[2])  {
            bg.thresholds <- c( bg.thresholds, data[neg,i][j][i] )   ## line specific
         }
    }
            
            
    data.bg.subtracted <- matrix( nrow=nrow(data), ncol=ncol(data) )
    rownames(data.bg.subtracted) <- rownames(data)
    colnames(data.bg.subtracted) <- colnames(data)
    i=1;
    for(i in 1:length(bg.thresholds))  {
         data.bg.subtracted[,i] <- data[,i] - bg.thresholds[i] 
     }

    return(data.bg.subtracted)


}  ## end function  

######################################################        
