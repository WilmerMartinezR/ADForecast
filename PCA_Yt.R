
################################################
#   Principal components analysis (PCA)
################################################

############################################################################
###   Functions
############################################################################
# Test unit root
DF.Pval <- function(Xs){
  ADCheck <- NULL
  for(indi in 2:(ncol(Xs)-2)){
    #indi <- 2
    #plot(as.ts(Xs[,indi]))
    adPval <-  tseries::adf.test(na.omit(as.ts(Xs[,indi]) ))$p.value
    ADCheck <- rbind(ADCheck, cbind(names(Xs)[indi], adPval))
  }
  ADCheck <- data.frame(ADCheck)
  colnames(ADCheck) <- c("Ind", "DFpval")
  ADCheck$DFpval <- as.numeric(ADCheck$DFpval)
  return(ADCheck)
}
############################################################################
# Differenciation
Diff.Xs <- function(Xs, ADCheck){
  Xs2 <- data.frame(Xs)
  Xs.d <- Xs2[-1,] 
  for(indi in 2:(ncol(Xs)-2)){ 
    #indi=2
    if(ADCheck$DFpval[indi-1]>0.05){
      Xs.d[,indi] <- diff(Xs2[,indi])
    }
  }
  return(Xs.d)
}
############################################################################


################################################
#   Reading data
################################################

indiAgre <- which(GB1.VA$Codigo<10)
GB2 <- GB1.VA[-indiAgre,] # excluye infla total
fechaFin0 <- which(names(GB2) == FechaFinSample) 
if(fechaFin0 < ncol(GB2)) GB2 <- GB2[,-c((fechaFin0+1):ncol(GB2))] 

GB2 <- data.frame(t( GB2[,-c(1:5)] ) )
GB2 <- data.frame( GB2[,1], GB2)
GB2$GB2...1. <- rownames(GB2)
GB2$AO <- as.numeric(substr(rownames(GB2), 5,8))
GB2$Mes <- substr(rownames(GB2), 1,3)

############################################################################
# Differentiation transformation
############################################################################

# First running
Xs <- GB2
ADCheck <- DF.Pval(Xs)
any(ADCheck$DFpval > 0.05)
length(which(ADCheck$DFpval > 0.05))
Xs.d <- Diff.Xs(Xs, ADCheck)

# Second running
#class(Xs.d) <- class(Xs)
ADCheck2 <- DF.Pval(Xs.d)
any(ADCheck2$DFpval > 0.05)
length(which(ADCheck2$DFpval > 0.05))
Xs.d2 <- Diff.Xs(Xs.d, ADCheck2)

# Third running
#class(Xs.d2) <- class(Xs)
ADCheck3 <- DF.Pval(Xs.d2)
any(ADCheck3$DFpval > 0.05)
length(which(ADCheck3$DFpval > 0.05))
#Xs.d3 <- Diff.Xs(Xs.d2, ADCheck3)


#############################################################################
# Replace last NA with the first observed value
Xs.d3 <- data.frame(Xs.d2)
for(indi in 2:(ncol(Xs.d2)-2)){ 
  #indi=188
  idNA <- which(is.na(Xs.d2[,indi]))
  if(any(idNA < 50)){ 
    idNA <- idNA[idNA < 50]
    Xs.d3[idNA,indi] <- Xs.d2[max(idNA)+1,indi]  
  }
}

Xs.d4 <- na.omit(Xs.d3)
pca1 <- prcomp(Xs.d4[,-c(1,ncol(Xs.d4)-1, ncol(Xs.d4))], scale = TRUE)
AA <- summary(pca1)

# Componentes

# Chequeo valor propio > 1
#eigVals <- eigen(cor(Xs.d4[,-c(1,ncol(Xs.d4)-1, ncol(Xs.d4))]))
#eigVals$values > 1
#ProVarAcum <- cumsum((pca1$sdev*pca1$sdev)/sum(pca1$sdev*pca1$sdev))
#ProVarAcum[eigVals$values > 1]

PcX <-pca1$x


# biplot(pca1, scale = 0)
# #calculate total variance explained by each principal component
# results <- pca1
# var_explained = results$sdev^2 / sum(results$sdev^2)
# 
# #create scree plot
# qplot(c(1:53), var_explained) +
#   geom_line() +
#   xlab("Principal Component") +
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot") +
#   ylim(0, 1)

###################################
# Define names and dates

# The first components that acummulate the CritACP% (defines in M0.R) of the variability
ProVarAcum <- cumsum((pca1$sdev*pca1$sdev)/sum(pca1$sdev*pca1$sdev))

YT <- data.frame(t(PcX[,c(which(ProVarAcum < CritACP)) ]))  
YT <- data.frame(Codigo=1:nrow(YT), YT)  

Maos <- paste(Xs.d4$Mes, Xs.d4$AO,sep = "_")

colnames(YT) <- c("Codigo",Maos)

YT1 <- YT

