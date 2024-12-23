

######################################################################
##   Functions
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

# This function calculates the forecast from the disagregates for the total
#  the total should define as 0 for the code value
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


###################################################################
#     The desaggregate performance in plots
###################################################################


file.read0 <- c("Fores_M0_PCA_TARYt_Mar_2016_Jun_2023.txt",
                "Fores_ARp_direct_Mar_2016_Jun_2023.txt",
                "Fores_M0_PCA_TARYt_Mar_2016_Jun_2023.txt",
                "Fores_M0_PCA_TARYt_Mar_2016_Jun_2023.txt",
                "Fores_RidgeYt_Mar_2016_Jun_2023.txt",
                "Fores_LassoYt_Mar_2016_Jun_2023.txt",
                "Fores_RFYt_Mar_2016_Jun_2023.txt")

file.read <- c(paste(substr(file.read0,1,5),Country,substr(file.read0,6,1000000L),sep = ""))

NomVar <- c("M0",
            "ARp", 
            "TAR",
            "PY0", 
            "RidY",
            "LasY",
            "RFY")

Jtas <- c(1:length(NomVar) )


#############################################################################  
#############################################################################

Agreg <- c(0,6) 
Letra <- list("", c("b","s")) 
NameOut <- c("Headline", "Core")
Rango <- RMSERango <- matrix(NA, nrow = length(Agreg), ncol = 2)
ValMax <- 50
ValMaxRMSE <- 10

for(ji in 1:length(Agreg)){
  #ji <- 3
  NumSG <- Agreg[ji]
  LetraSG <- Letra[[ji]]

  if(NumSG>0){
    filtraSG <- GB1.VA[,c("Codigo","Tipo")] %>% dplyr::filter(Tipo %in% LetraSG)
  }
  
  for(j in Jtas){
    # the first run different
    #j=1 
    
    if(j==1){
      #M0
      outputForest <- read.table(paste(path.out, file.read[j], sep = ""))
      AB2 <- outputForest[outputForest$Mod1 == 1,]
      colnames(AB2)[which(colnames(AB2) == "Pesos")] <- c("Pesos2018")
    }else if(j==2){
      #ARp
      outputForest <- read.table(paste(path.out, file.read[j], sep = ""))
      AB2 <- outputForest
      colnames(AB2)[which(colnames(AB2) == "Pesos")] <- c("Pesos2018")
    }else if(j==3){
      #TAR
      outputForest <- read.table(paste(path.out, file.read[j], sep = ""))
      AB2 <- outputForest[outputForest$Mod1 == 3,]
      colnames(AB2)[which(colnames(AB2) == "Pesos")] <- c("Pesos2018")
    }else if(j==4){
      #PY0
      outputForest <- read.table(paste(path.out, file.read[j], sep = ""))
      AB2 <- outputForest[outputForest$Mod1 == 2,]
      colnames(AB2)[which(colnames(AB2) == "Pesos")] <- c("Pesos2018")
    }else{
      outputForest <- read.table(paste(path.out, file.read[j], sep = ""))
      AB2 <- outputForest
      colnames(AB2)[which(colnames(AB2) == "Pesos")] <- c("Pesos2018")
    }
    
    if(NumSG>0){
      AB2 <- AB2 %>% dplyr::filter(Codigo %in% filtraSG$Codigo | Codigo == NumSG)
      AB2$Codigo <- ifelse(AB2$Codigo == NumSG, 0, AB2$Codigo)
    }
    
    ##################################################################
    ###### Joint obs, h, and Fores from disagr
    if(j==1){
      Obsh1 <- Resport1A(AB2)
      tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(Obsh1)))
      tvar1 <- "^h"; name1 <- c(grep(tvar1, names(Obsh1)))
      tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(Obsh1)))
      
      
      #colnames(Obsh1)[c(name1,name2)] <- paste(colnames(Obsh1)[c(name1,name2)], 
      #                                         "_",NomVar[j],sep="") 
      Obsh1$Mod <- j
    }else{
      Obsh1a <- Resport1A(AB2)
      tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(Obsh1a)))
      tvar1 <- "^h"; name1 <- c(grep(tvar1, names(Obsh1a)))
      tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(Obsh1a)))
      
      #colnames(Obsh1a)[c(name1,name2)] <- paste(colnames(Obsh1a)[c(name1,name2)], 
      #                                          "_",NomVar[j],sep="") 
      Obsh1a$Mod <- j
      Obsh1 <- rbind(Obsh1, Obsh1a)
      
    }
    
    ### Performance
    AC1 <- Resport1(AB2, fechaPand0 = FALSE, fecPand = "Dic_2019",  
                    antes = TRUE, i=hh, filterMod = FALSE, 
                    GuardaPlot = F, fileN)
    
    if(j==1){
      AC1A <- AC1[,c("h", "RMSEA.0","RMSED.1")]
      AC1A$Mod <- j
    }else{
      AC1Aa <- AC1[,c("h","RMSEA.0","RMSED.1")]
      AC1Aa$Mod <- j
      AC1A <- rbind(AC1A, AC1Aa)
    }
    
  }
  
  ############################################################
  ####  Here we create a graphical report for each item
  ####  comparing the forecast (h=1,3,6,12), and the observed
  ####  Also, we create a plot with the RMSE per item 
  ############################################################
  
  AB3F <- Obsh1
  AB3F$ModT <- as.factor(AB3F$Mod)
  
  HH <- c(1,3,6,12)
  
  for(hh in HH){
    #hh=1
    
    # define title
    tit <- paste(#NumSG," - ", 
      round(GB1.VA[GB1.VA$Codigo==NumSG,3],2)," - ",
      NameOut[ji],sep="" )
    
    # define name of the plot to store
    nom1Val <- paste(pathplotStore, NameOut[ji],"_h_",hh,"Fdesa_",Country,".png",sep = "")
    
    tvar20 <- "^Forest_AfD"; name20 <- c(grep(tvar20, names(AB3F)))
    tvar21 <- "^h"; name21 <- c(grep(tvar21, names(AB3F)))
    
    if(hh==3){
      colnames(AB3F)[name20[1]] <- "Forest_AfD13"
      colnames(AB3F)[name20[3]] <- "Forest_AfD1"
      
      colnames(AB3F)[name21[1]] <- "h13"
      colnames(AB3F)[name21[3]] <- "h1"
    }else if(hh==6){
      colnames(AB3F)[name20[3]] <- "Forest_AfD14"
      colnames(AB3F)[name20[6]] <- "Forest_AfD1"
      
      colnames(AB3F)[name21[3]] <- "h14"
      colnames(AB3F)[name21[6]] <- "h1"
    }else if(hh==12){
      colnames(AB3F)[name20[6]] <- "Forest_AfD15"
      colnames(AB3F)[name20[12]] <- "Forest_AfD1"
      AB3F[AB3F$Mod==2,"Forest_AfD1"] <- NA
      
      colnames(AB3F)[name21[6]] <- "h15"
      colnames(AB3F)[name21[12]] <- "h1"
      AB3F[AB3F$Mod==2,"h1"] <- NA
    }
    
    AB3F[,"Forest_AfD1"] <- ifelse(abs(AB3F$Forest_AfD1) > ValMax,NA, AB3F$Forest_AfD1)
    Rango[ji,] <- round(range(AB3F$Forest_AfD1, AB3F$h1, 
                              AB3F$Obs1,na.rm = T),0)  
    
    AB3Fa <- AB3F %>% dplyr::filter(Mod == 1)
    AB3Fa$Obs1 <- lag(AB3Fa$Obs1, hh)
    pp1<- ggplot(AB3F[AB3F$ordenf>hh,], aes(x=ordenf, y=Forest_AfD1, group=ModT)) + 
      geom_line(aes(color=ModT)) +
      geom_line(aes(y=rep(AB3Fa$Obs1[-c(1:hh)],max(AB3F$Mod)) )) +
      scale_color_manual(name='Model', labels=c(c(NomVar,"Best","Obs")),
                         values=c( 'darkgreen', 'red', 'steelblue',
                                   "gray", "green","blue","orange",
                                   "pink","brown","yellow",
                                   "#66CDAA","purple","#C6E2FF","#00FFFF",
                                   '#FFD39B',
                                   "black"))+
      labs(y=paste("h = ",hh, sep = ""), x="Time (Months)")
    
    ### Dates in x axes
    Ms <- substr(AB3Fa$fechaf[-c(1:hh)],1,3)
    Year <- substr(AB3Fa$fechaf[-c(1:hh)],5,8)
    Mes <- c()
    for(i in 1:length(Ms)) Mes[i] <- which(Ms[i]==months)
    fechas <- paste(Mes,"/",Year,sep = "")
    posi <- which(Mes==12)
    
    pp1+theme_bw()+ggtitle(tit)+
      scale_x_discrete(labels=c(fechas[posi]),
                       limits=c(posi+hh))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"),
            legend.title = element_text(size=25,face="bold"),
            legend.text = element_text(size=18),
            plot.title = element_text(size=28,face="bold"))+
      ylim(-10,20)
    ggsave(nom1Val, width = 10, height = 7)
  
    ##########################################################
    ## Aggregate plot
    ##########################################################
    
    nom1Agre <- paste(pathplotStore, NameOut[ji],"_h_",hh,"_",Country,".png",sep = "")
    
    AB3F[,"h1"] <- ifelse(abs(AB3F$h1) > ValMax,NA, AB3F$h1)
    
    pp1<- ggplot(AB3F[AB3F$ordenf>hh,], aes(x=ordenf, y=h1, group=ModT)) + 
      geom_line(aes(color=ModT)) +
      geom_line(aes(y=rep(AB3Fa$Obs1[-c(1:hh)],max(AB3F$Mod)) ) ) +
      scale_color_manual(name='Model', labels=c(c(NomVar,"Best","Obs")),
                         values=c('darkgreen', 'red', 'steelblue',
                                  "gray", "green","blue","orange",
                                  "pink","brown","yellow",
                                  "#66CDAA","purple","#C6E2FF","#00FFFF",
                                  '#FFD39B',
                                  "black"))+
      labs(y=paste("h = ",hh, sep = ""), x="Time (Months)")
    
    ### Defines dates on x axis
    Ms <- substr(AB3Fa$fechaf[-c(1:hh)],1,3)
    Year <- substr(AB3Fa$fechaf[-c(1:hh)],5,8)
    Mes <- c()
    for(i in 1:length(Ms)) Mes[i] <- which(Ms[i]==months)
    fechas <- paste(Mes,"/",Year,sep = "")
    posi <- which(Mes==12)
    
    pp1+theme_bw()+ggtitle(tit)+
      scale_x_discrete(labels=c(fechas[posi]),
                       limits=c(posi+hh))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"),
            legend.title = element_text(size=25),
            legend.text = element_text(size=18),
            plot.title = element_text(size=28,face="bold"))+
      ylim(-10,20)
    ggsave(nom1Agre, width = 10, height = 7)
    
  }
  
}
###################################################

