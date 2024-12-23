
#########################################################################################
### Function to calculate year-on.year inflation

vAnual <- function(varia, orden){
vari <- NULL
for(i in (orden+1):length(varia)){
	if(varia[(i-orden)] > 0 & varia[(i)] > 0){
		vari <- rbind(vari, (varia[i]/varia[(i-orden)] - 1)*100) 
	}else vari <- rbind(vari, NA)
}
return(vari)
}
#########################################################################################
###   This function re-scale, calculate year-on.year inflation and standardize ##########
# escala = 1 if re-scale, 0 o.w.
# vAnual = 1 if year-on.year inflation, 0 o.w.
# orden = to define the differential order 4 quarterly, 12 monthly,
# estand = 1 standardize, 0 o.w.
vA <- function(vari, escala, vAnual, orden, estand){
  if(escala == 1) fact <- vari+abs(min(vari))+100 else fact <- vari
  if(vAnual == 1) fact <- vAnual(fact, orden) else fact <- fact
  if(estand == 1) fact <- (fact-mean(fact))/sd(fact) else fact <- fact
  fact <- data.frame(as.vector(fact))
  return(fact)
}

########################################################################################
# Non-autocorrelation test
########################################################################################
No.correla <- function(x, lags){
  Ljung_Box <- NULL
  for(i in 1:lags){
    Ljung_Box1 <- Box.test(x, lag = i, type="Ljung")$p.value
    Ljung_Box <- rbind(Ljung_Box, cbind(i,Ljung_Box1))
  }
  colnames(Ljung_Box) <- c("Rezago","p-valor");
  return(Ljung_Box)
}


##############################################################################
##### Function to calculate Hodrick and Prescott                         #####
##############################################################################

## The ""standard-values"" for Lambda are 
## 100 for yearly data, 1600 for quarterly data  and 14400 for monthly data.

hpfilter <- function(x,lambda){
  eye <- diag(length(x))
  result <- solve(eye+lambda*crossprod(diff(eye,lag=1,d=2)),x)	
  return(result)
}

############################################################################
###   Inverse differenciation transformation 

Integra <- function(yyd, yy0){
  yy1 <- yyd[1] + yy0
  if(length(yyd) >= 2) for(i in 2:length(yyd)) yy1[i] <- yyd[i] + yy1[i-1]
  return(yy1)
}

##########################################################
# Test unit root
DF.Pval <- function(Xs){
  ADCheck <- NULL
  for(indi in 1:(nrow(Xs))){
    #indi <- 1
    #plot(as.ts(Xs[,indi]))
    adPval <- adf.test(na.omit(as.ts(t(Xs[indi,-c(1:5)]) ) ))$p.value
    ADCheck <- rbind(ADCheck, cbind(Xs[indi,1], adPval))
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
  for(indi in 1:(nrow(Xs))) { 
    #indi=1
    if(ADCheck$DFpval[indi]>0.05){
      #Xs.d[indi,(6+diFF):ncol(Xs.d)] <- diff(Xs2[indi,-c(1:(5+diFF-1) )])
      #Xs.d[indi,6:(6+diFF-1)] <- Xs.d[indi,6+diFF]
      Xs2[indi,7:ncol(Xs)] <- diff(t(Xs2[indi,-c(1:5)]))    
    }
  }
  Xs.d <- Xs2[,-6]
  return(Xs.d)
}


#################################################
# Transform time series to I(0)
Estaciona <- function(GB1){
  
  GB2 <- GB1[-which(GB1$Codigo<10),] # exclude aggregates
  
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
  #Xs.d2 <- cbind(GB2[,c(1:6)], Xs.d2)
  return(Xs.d2)
}


######################################################################
##   Eval functions
######################################################################

Resport1 <- function(outputForest, fechaPand0 = TRUE, fecPand = "Dic_2019", 
                     antes=TRUE, i=1, filterMod=TRUE,
                     GuardaPlot=FALSE, file.out, tab=TRUE){
  
  
  #outputForest = AB2
  if(filterMod){
    outputForestE0 <- outputForest[outputForest$Mod1 == 2,]
  }else outputForestE0 <- outputForest 
  
  outputForestE0 <- outputForestE0[order(outputForestE0$Codigo),]
  outputForestE0 <- outputForestE0 %>% dplyr::filter(Codigo == 0 | Codigo > 5)
  
  
  if(fechaPand0){
    fechaPand <- min(outputForestE0[outputForestE0$fechaf == fecPand,"ordenf"])
    if(antes){
      outputForestE0 <- outputForestE0 %>% dplyr::filter(ordenf <= fechaPand)
    }else outputForestE0 <- outputForestE0 %>% dplyr::filter(ordenf > fechaPand)
    perform.Ys.Full <- CheckResiFores2(outputForestE = outputForestE0)
  }else perform.Ys.Full <- CheckResiFores2(outputForestE = outputForestE0)
  
  # Agregando errores
  # perform.Ys.Full$PerformanceAD
  # Errores desde los agregados
  #  perform.Ys.Full$PerformanceAyAfD
  ForestW <- perform.Ys.Full$ForesAyAfD
  
  tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(ForestW)))
  tvar1 <- "^h"; name1 <- c(grep(tvar1, names(ForestW)))
  #i=1
  if(GuardaPlot){
    
    png(file.out)
    
    op <- par(mfrow=c(1,2),las=3,mar=c(4,3.5,2,0)+.5, mgp=c(2.2,0.3,0), cex.axis=1.1, cex.lab=1.2)
    
    ## Filtra el pronostico usando el agregado
    ForestW.F <- ForestW[ForestW$Codigo == 0,]
    plot(ForestW.F$ordenf, ForestW.F[,name0[i]], 
         ylim=range(ForestW[,c(name0[i],name1[i])],na.rm = T),pch=20,
         ylab="Inflación", xlab="", axes=F)
    axis(2)
    axis(1, ForestW.F$ordenf, substr(unique(ForestW.F$fechaf),5,8) )
    # pronósticos del total
    lines(ForestW.F$ordenf, ForestW.F[,name1[i]], col="red")
    
    # pronósticos del agregado
    ForestW.AD <- ForestW[ForestW$Codigo == 1,]
    lines(ForestW.AD$ordenf, ForestW.AD[,name1[i]], col="blue",lty=2,lwd=2)
    
    ### Residuals
    tvar2 <- "^Resi"; name2 <- c(grep(tvar2, names(ForestW)))
    i=1
    ForestW.F <- ForestW[ForestW$Codigo == 0,]
    plot(ForestW.F$ordenf, ForestW.F[,name2[i]], 
         ylim=range(ForestW[,c(name2[i])],na.rm = T), col="red",
         ylab="Errores de pronósticos", xlab="", axes=F, type = "o")
    axis(2)
    axis(1, ForestW.F$ordenf, substr(unique(ForestW.F$fechaf),5,8) )
    
    ForestW.AD <- ForestW[ForestW$Codigo == 1,]
    lines(ForestW.AD$ordenf, ForestW.AD[,name2[i]], col="blue",lty=2,lwd=2, type = "o")
    #name2 <- c(grep(tvar2, names(ForestWBest)))
    #lines(ForestW.ADBest$ordenf, ForestW.ADBest[,name2[i]], col="green")
    
    par(op)
    dev.off()
  }
  
  Perform <- merge(perform.Ys.Full$PerformanceAD,  perform.Ys.Full$PerformanceAyAfD,
                   by.x = "h", by.y = "h")
  Perform <- Perform[,-ncol(Perform)]
  colnames(Perform) <- c("h","RMSEA.0","RMSED.0","Acc.0","PvalDM.0","SsFM.0",     
                         "RMSEA.1","RMSED.1","Acc.1","PvalDM.1")
  if(tab) Perform
  
}

# This function calculates the forecast from the disagregates for the aggregate
#  the aggregate should define as 0 for the code value
Resport1A <- function(outputForest){
  
  outputForestE0 <- outputForest 
  
  outputForestE0 <- outputForestE0[order(outputForestE0$Codigo),]
  outputForestE0 <- outputForestE0 %>% dplyr::filter(Codigo == 0 | Codigo > 5)
  
  perform.Ys.Full <- CheckResiFores2(outputForestE = outputForestE0)
  
  ForestW <- perform.Ys.Full$ForesAyAfD
  
  tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(ForestW)))
  tvar1 <- "^h"; name1 <- c(grep(tvar1, names(ForestW)))
  tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(ForestW)))
  #i=1
  
  Out <- ForestW[,c(which(names(ForestW) == "fechaf"),
                    which(names(ForestW) == "ordenf"),
                    name0, name1, name2)]
  return(Out)
}
