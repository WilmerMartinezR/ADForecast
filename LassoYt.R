
###############################################
#    Lasso Regression
###############################################

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


################################################
#   Data
################################################

# The object GB1.VA comes from GB0.R
datos1 <- GB1.VA

# Take dates according to the PCA run
Ini0 <- which( names(datos1)  == names(YT1)[2])
Fin0 <- which( names(datos1)  == names(YT1)[ncol(YT1)])

datos1 <- datos1[,c(1:5, Ini0:Fin0)]

# Crea columnas vacias para la evaluacion de pronosticos
# Create empty columns to store forecast
datosF0 <- datos1[,c(6:(6+11))]
datosF0[,1:12] <- NA

months <- c("Ene","Feb","Mar","Abr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dic")
last.year <- as.numeric(substr(names(datos1)[ncol(datos1)],5,8) )
last.month <- names(datos1)[ncol(datos1)]

aos <- rep(last.year:(last.year+1),each=12)
Maos <- paste(months,aos,sep = "_")

lastMObs <- which(Maos == last.month) 
names(datosF0) <- c(Maos[(lastMObs+1):(lastMObs+12)])

datos <- cbind(datos1, datosF0)

datosX <- Estaciona(GB1.VA)

##################################################
# Corre modelos y pronostica
# Adjust models and forecast

outputForestE <- outputForestFul <- NULL
IniT <- Sys.time()
for(mW in 0:(ncol(datos)-6-W.Wide-h)){
  
  fecha.ini <- names(datos)[6]
  fecha.fin <- names(datos)[6+W.Wide+mW]
  
  output <- T_ModResiForest_RL3(datos, fecha.ini, fecha.fin, 
                                h, KindRegre = 1, datosX, datosX2=NULL)
  
  IniPronos <- which(names(datos)==fecha.fin)
  datosFM <- data.frame(datos[,(IniPronos+1):(IniPronos+h)])
  datosFM$Codigo <- datos[,1]
  names(datosFM) <- c(paste("Obs",1:12,sep = ""),"Codigo")
  
  output2 <- merge(output$BestModel, datosFM, by.x="Codigo", by.y="Codigo")
  
  outputForestE <- rbind(outputForestE, output2)
  
}  
FinT <- Sys.time()
FinT - IniT

fechaf = unique(outputForestE$fechaf)
fechaOrd <- data.frame(fechaf = fechaf, ordenf = 1:length(fechaf))

outputForestE <- merge(datos[,c("Codigo","Pesos")], outputForestE,
                       by.x = "Codigo", by.y = "Codigo")

outputForestE <- merge(fechaOrd, outputForestE,
                       by.x = "fechaf", by.y = "fechaf")
outputForestE <- outputForestE[order(outputForestE$ordenf),]

# TO store the results
file.out <- paste("Fores",Country,"_LassoYt_", outputForestE$fechaf[1], "_",                   
                  outputForestE$fechaf[nrow(outputForestE)],".txt",sep = "")
write.table(outputForestE, paste(path.out, file.out, sep = "") 
            ,sep = "\t")