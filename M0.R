
### Run PathUS

# #################################################
# CPI data

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
#months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dec")
last.year <- as.numeric(substr(names(datos1)[ncol(datos1)],5,8) )
last.month <- names(datos1)[ncol(datos1)]

aos <- rep(last.year:(last.year+1),each=12)
Maos <- paste(months,aos,sep = "_")

lastMObs <- which(Maos == last.month) 
names(datosF0) <- c(Maos[(lastMObs+1):(lastMObs+12)])

datos <- cbind(datos1, datosF0)

##################################################
# Corre modelos y pronostica
# Adjust models and forecast

outputForestE <- outputForestFul <- NULL
IniT <- Sys.time()
for(mW in 0:(ncol(datos)-6-W.Wide-h)){

  fecha.ini <- names(datos)[6]
  fecha.fin <- names(datos)[6+W.Wide+mW]
  
  output <- T_ModResiForestARX3_selectTAR(datos, fecha.ini, fecha.fin, h, 
                                          datosX=NULL, datosX2=YT1,
                                          datosX3=NULL,
                                          TARMod=FALSE, IPPindv = FALSE,
                                          CompSeaso=FALSE)
  
  IniPronos <- which(names(datos)==fecha.fin)
  datosFM <- data.frame(datos[,(IniPronos+1):(IniPronos+h)])
  datosFM$Codigo <- datos[,1]
  names(datosFM) <- c(paste("Obs",1:12,sep = ""),"Codigo")
  
  #output2 <- merge(output$BestModel, datosFM, by.x="Codigo", by.y="Codigo")
  outputF <- merge(output$FullModels, datosFM, by.x="Codigo", by.y="Codigo")
  
  #outputForestE <- rbind(outputForestE, output2)
  outputForestFul <- rbind(outputForestFul, outputF)
  
}  
FinT <- Sys.time()
FinT - IniT

##################################################
# Organize output
##################################################

outputForestE <- outputForestFul

fechaf = unique(outputForestE$fechaf)
fechaOrd <- data.frame(fechaf = fechaf, ordenf = 1:length(fechaf))

outputForestE <- merge(datos[,c("Codigo","Pesos")], outputForestE,
                       by.x = "Codigo", by.y = "Codigo")

outputForestE <- merge(fechaOrd, outputForestE,
                       by.x = "fechaf", by.y = "fechaf")
outputForestE <- outputForestE[order(outputForestE$ordenf),]

file.out <- paste("Fores",Country,"_M0_PCA_TARYt_", outputForestE$fechaf[1], "_",                  
                  outputForestE$fechaf[nrow(outputForestE)],".txt",sep = "")

write.table(outputForestE, paste(path.out, file.out, sep = "") 
            ,sep = "\t")
