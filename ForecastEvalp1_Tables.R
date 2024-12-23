
##############################################################################
###  Files
######################################################################


### Calculate the relative RMSFE 
RMSE_rela <- TRUE
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

hh = 1
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
    
  ##################################################################
  ###### Joint obs, h, and Fores from disagr
  
  if(j==1){
    Obsh1 <- Resport1A(AB2)
    tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(Obsh1)))
    tvar1 <- "^h"; name1 <- c(grep(tvar1, names(Obsh1)))
    tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(Obsh1)))
    
    colnames(Obsh1)[c(name1,name2)] <- paste(colnames(Obsh1)[c(name1,name2)], 
                                             "_",NomVar[j],sep="") 
  }else{
    Obsh1a <- Resport1A(AB2)
    tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(Obsh1a)))
    tvar1 <- "^h"; name1 <- c(grep(tvar1, names(Obsh1a)))
    tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(Obsh1a)))
    
    colnames(Obsh1a)[c(name1,name2)] <- paste(colnames(Obsh1a)[c(name1,name2)], 
                                              "_",NomVar[j],sep="") 
    
    Obsh1 <- cbind(Obsh1, Obsh1a[,-c(1:2,name0)])
    
  }
  #}
  ##################################################################
  ###### Performance in terms of RMSE
  
  fileN <- NULL
  AC1 <- Resport1(AB2, fechaPand0 = FALSE, fecPand = "Dic_2019",  
                  antes = TRUE, i=hh, filterMod = FALSE, 
                  GuardaPlot = F, fileN)
  
  if(j==1){
    AC1A <- AC1[,c("h","RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1A) <- c("h", 
                        paste(c("RMSEA.0","RMSED.0","RMSED.1"), NomVar[j],sep = ""))
  }else{
    AC1Aa <- AC1[,c("RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1Aa) <- paste(c("RMSEA.0","RMSED.0","RMSED.1"), NomVar[j],sep = "")
    AC1A <- cbind(AC1A, AC1Aa)
  }
  
  ##########################################################
  ######  Results before COVID-19
  
  #nk=2
  #fileN = paste(pathplot,Namesplots[j+1,nk],".png",sep = "")
  fileN = NULL
  AC1 <- Resport1(AB2, fechaPand0 = T, fecPand = "Dic_2019",  
                  antes = TRUE, i=hh, filterMod = FALSE, 
                  GuardaPlot = F, fileN)
  
  # Performance
  if(j==1){
    AC1B <- AC1[,c("h","RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1B) <- c("h", 
                        paste(c("RMSEA.0","RMSED.0","RMSED.1"), NomVar[j],sep = ""))
  }else{
    AC1Ba <- AC1[,c("RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1Ba) <- paste(c("RMSEA.0","RMSED.0","RMSED.1"), NomVar[j],sep = "")
    AC1B <- cbind(AC1B, AC1Ba)
  }
  
  ################################################################
  ######  Results after COVID-19

  AC1 <- Resport1(AB2, fechaPand0 = T, fecPand = "Dic_2019",  
                  antes = F, i=hh, filterMod = FALSE, 
                  GuardaPlot = F, fileN)
  
  # Performance
  if(j==1){
    AC1C <- AC1[,c("h","RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1C) <- c("h", 
                        paste(c("RMSEA.0","RMSED.0","RMSED.1"), NomVar[j],sep = ""))
  }else{
    AC1Ca <- AC1[,c("RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1Ca) <- paste(c("RMSEA.0","RMSED.0","RMSED.1"), NomVar[j],sep = "")
    AC1C <- cbind(AC1C, AC1Ca)
  }
  
}

if(RMSE_rela){
  
  Namestab1 <- c("All_A", "Antes_A", "Desp_A")
  
  tvar1 <- "^RMSEA"; name1 <- c(grep(tvar1, names(AC1A)))
  AC1A.0 <- AC1A[,name1]
  colnames(AC1A.0) <- NomVar

  tvar1 <- "^RMSEA"; name1 <- c(grep(tvar1, names(AC1B)))
  AC1B.0 <- AC1B[,name1]
  colnames(AC1B.0) <- NomVar
  
  tvar1 <- "^RMSEA"; name1 <- c(grep(tvar1, names(AC1C)))
  AC1C.0 <- AC1C[,name1]
  colnames(AC1C.0) <- NomVar
  
  ## Names of the results
  Namestab2 <- c("All_D", "Antes_D", "Desp_D")
  
  tvar1 <- "^RMSED.1"; name1 <- c(grep(tvar1, names(AC1A)))
  AC1A.0D <- AC1A[,name1]
  colnames(AC1A.0D) <- NomVar
  
  tvar1 <- "^RMSED.1"; name1 <- c(grep(tvar1, names(AC1B)))
  AC1B.0D <- AC1B[,name1]
  colnames(AC1B.0D) <- NomVar
  
  tvar1 <- "^RMSED.1"; name1 <- c(grep(tvar1, names(AC1C)))
  AC1C.0D <- AC1C[,name1]
  colnames(AC1C.0D) <- NomVar
}

########################################################
#    Functions to evaluation
### calcula el DM test

Eval_Agre <- TRUE
if(Eval_Agre){
  
  Eval_RMSE.DM <- function(AC1A.0, AC1A.0D, h, Obsh1, NomVar,
                           typeFores = "_A",
                           saveTab = TRUE, path = pathplot, 
                           NomSale = Namestab2[3],
                           ListLags = c(1,3,6,9,12),
                           Samplef=FALSE, antes=FALSE,fecPand = "Dic_2019"){
    
    # AC1A.0 is the RMSE table for the models 
    # ref_i is the reference position in AC1A.0 (RMSE table 
    #        of the models by column), this order depend on file.read object   
    # Ex:ref_i <- 1
    # Obsh1 is the data frame with the observations and the forecast
    # NomVar is the list of models
    # typeFores = "_A" if we create agregate table typeFores = "_D" otherwise
    
    if(typeFores == "_A"){
      AC1A.0A <- AC1A.0
      mods <- 1:ncol(AC1A.0A)

    }else{
      AC1A.0A <- AC1A.0D
      mods <- 1:ncol(AC1A.0A)
    }
    ## Calculate DM test, potentially we could use Martinez et.al. test
    
    if(Samplef){
      fechaPand <- min(Obsh1[Obsh1$fechaf == fecPand,"ordenf"])
      if(antes){
        Obsh1 <- Obsh1 %>% dplyr::filter(ordenf <= fechaPand)
      }else Obsh1 <- Obsh1 %>% dplyr::filter(ordenf > fechaPand)
    }
    
    tvar0 <- paste("*","Obs",sep=""); name0 <- c(grep(tvar0, names(Obsh1)))
    A00 <- Obsh1[,name0]
    
    # We calculate which model is the one that has minimum 
    #  RMSE to use as reference to compare performance
    
    ref_i <- c()
    for(j in 1:nrow(AC1A.0)) ref_i <- c(ref_i, 
                                        which(AC1A.0[j,] == min(AC1A.0[j,]))[1])
    
    ### DM test comparing with the desagregates of the reference
    PvalDM <- c()
    for(fw in 1:h){ 
      #fw <- 3
      tvar1 <- paste("*",NomVar[ref_i[fw]],"$",sep=""); 
      name1 <- c(grep(tvar1, names(Obsh1)))
      A0 <- Obsh1[,name1]
      PvalDM <- c(PvalDM, 
                  dm.test(na.omit(A00[,c(fw)]-A0[,c(fw)]), 
                          na.omit(A00[,c(fw)]-A0[,c(fw+h)]), 
                          #h=fw, alternative ="greater",
                          h=fw, alternative ="two.sided",
                          varestimator = c("bartlett"))$p.value)
    }
    DM_comp <- data.frame(PvalDM)
    #Noms0 <- c(paste(NomVar[ref_i],c("D"),sep="_"))
    
    Noms0 <- c()
    #ij in mods
    for(ij in mods){
      #ij = 7
      tvar2 <- paste("*",NomVar[ij],"$",sep=""); 
      name2 <- c(grep(tvar2, names(Obsh1)))
      A1 <- Obsh1[,name2]
      #if(length(name2) > 24) break
      ### DM test comparing agregates
      
      PvalDM <- c()
      for(fw in 1:h){
        # This is to choose the best model per horizon
        tvar1 <- paste("*",NomVar[ref_i[fw]],"$",sep=""); 
        name1 <- c(grep(tvar1, names(Obsh1)))
        A0 <- Obsh1[,name1]
        
        if(mean(A0[,c(fw)] == A1[,c(fw)]) != 1){
          PvalDM <- c(PvalDM, 
                      dm.test(na.omit(A00[,c(fw)]-A0[,c(fw)]), 
                              na.omit(A00[,c(fw)]-A1[,c(fw)]), 
                              #h=fw, alternative ="greater",
                              h=fw, alternative ="two.sided",
                              varestimator = c("bartlett"))$p.value)
        }else PvalDM <- c(PvalDM, 1) 
      }
      
      DM_comp <- cbind(DM_comp, PvalDM)
      
      ### DM test comparing with the desagregates
      PvalDM <- c()
      for(fw in 1:h) PvalDM <- c(PvalDM, 
                                 dm.test(na.omit(A00[,c(fw)]-A0[,c(fw)]), 
                                         na.omit(A00[,c(fw)]-A1[,c(fw+h)]), 
                                         #h=fw, alternative ="greater",
                                         h=fw, alternative ="two.sided",
                                         varestimator = c("bartlett"))$p.value)
      DM_comp <- cbind(DM_comp, PvalDM)
      Noms0 <- c(Noms0, paste(NomVar[ij],c("A","D"),sep="_"))
    }
    
    # This is remove because it is duplicate
    DM_comp <- DM_comp[,-1]
    colnames(DM_comp) <- Noms0
    
    DM_comp2 <- DM_comp
    for(i in 1:nrow(DM_comp))
      for(j in 1:ncol(DM_comp))
        DM_comp2[i,j] <- ifelse(DM_comp2[i,j]>0.1,"",
                                ifelse(DM_comp2[i,j]<=0.1 & DM_comp2[i,j]>0.05,"*",
                                       ifelse(DM_comp2[i,j]<=0.05 & DM_comp2[i,j]>0.01,"**",
                                              #ifelse(DM_comp2[i,j]<=0.01 & DM_comp2[i,j]>0.001,"**",
                                              "***")))#)
    #typeFores <- "_A"
    tvar3 <- paste("*",typeFores,sep=""); 
    name3 <- c(grep(tvar3, names(DM_comp2)))
    DM_comp3 <- DM_comp2[,name3]
    for(i in 1:nrow(DM_comp3))
      for(j in 1:ncol(DM_comp3))
        DM_comp3[i,j] <- paste(round(AC1A.0A[i,j],2),DM_comp2[i,name3[j]],sep="")
    
    colnames(DM_comp3) <- names(AC1A.0A)
    #NomSale <- Namestab2[3]
    #Relativo <- 
    #ListLags <- c(1,3,6,9,12)
    if(saveTab){
      print(xtable(t(DM_comp3[ListLags,])), 
            file = paste(path, NomSale,".tex",sep = ""), 
            compress = FALSE)
    }else return(list(DM_comp = DM_comp3, BestMod_h = ref_i))
  }
  
  JointTabNew <- function(A, B, path, NomSale, 
                          ListLags = c(1,3,6,9,12)){
    A <- t(A[ListLags,]); B <- t(B[ListLags,])
    
    AB <- cbind(A, B)
    #rownames(AB)[ref_i] <- rownames(B)[ref_i]
    print(xtable(AB), 
          file = paste(path, NomSale,".tex",sep = ""), 
          compress = FALSE)
  }
  
  ########################################################
  ## Full evaluation set vs Agregates
  
  saveTab <- F
  A11 <- Eval_RMSE.DM(AC1A.0, AC1A.0D,h, Obsh1, NomVar,
                      typeFores = "_A",
                      saveTab, path = pathplot, 
                      NomSale = Namestab1[1])
  A1 <- A11$DM_comp
  BestA1 <- A11$BestMod_h[c(1,3,6,9,12)]
  ## Before pandemic evaluation set vs Agregates
  A22 <- Eval_RMSE.DM(AC1B.0, AC1B.0D, h, Obsh1, NomVar,
                      typeFores = "_A",
                      saveTab, path = pathplot, 
                      NomSale = Namestab1[2],
                      Samplef=T, antes=T,fecPand = "Dic_2019")
  A2 <- A22$DM_comp
  BestA2 <- A22$BestMod_h[c(1,3,6,9,12)]
  ## After pandemic evaluation set vs Agregates
  A33 <- Eval_RMSE.DM(AC1C.0, AC1C.0D,h, Obsh1, NomVar,
                      typeFores = "_A",
                      saveTab, path = pathplot, 
                      NomSale = Namestab1[3],
                      Samplef=T, antes=F,fecPand = "Dic_2019")
  A3 <- A33$DM_comp
  BestA3 <- A33$BestMod_h[c(1,3,6,9,12)]
  #########################################
  ## Full evaluation set vs Desagregates
  B1 <- Eval_RMSE.DM(AC1A.0, AC1A.0D,h, Obsh1, NomVar,
                     typeFores = "_D",
                     saveTab, path = pathplot, 
                     NomSale = Namestab2[1])$DM_comp
  ## Before pandemi evaluation set vs Desagregates
  B2 <- Eval_RMSE.DM(AC1B.0, AC1B.0D,h, Obsh1, NomVar,
                     typeFores = "_D",
                     saveTab, path = pathplot, 
                     NomSale = Namestab2[2],
                     Samplef=T, antes=T,fecPand = "Dic_2019")$DM_comp
  ## After pandemi evaluation set vs Desagregates
  B3 <- Eval_RMSE.DM(AC1C.0, AC1C.0D,h, Obsh1, NomVar,
                     typeFores = "_D",
                     saveTab, path = pathplot, 
                     NomSale = Namestab2[3],
                     Samplef=T, antes=F,fecPand = "Dic_2019")$DM_comp
  
  JointTabNew(A1, B1, path= pathplot, NomSale = Namestab1[1])
  JointTabNew(A2, B2, path= pathplot, NomSale = Namestab1[2])
  JointTabNew(A3, B3, path= pathplot, NomSale = Namestab1[3])
  
  
}

## Best models at horizons 1, 3, 6, 9, and 12

# all evaluation period
NomVar[BestA1]
# Before COVID-19
NomVar[BestA2]
# after COVID-19
NomVar[BestA3]

#############################################################################  
###################################
#   Checking subsets: SinAli_niReg (core). 
#############################################################################

NumSG <- 6; LetraSG <- c("b","s") 

filtraSG <- GB1.VA[,c("Codigo","Tipo")] %>% dplyr::filter(Tipo %in% LetraSG)

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
  
  AB2 <- AB2 %>% dplyr::filter(Codigo %in% filtraSG$Codigo | Codigo == NumSG)
  AB2$Codigo <- ifelse(AB2$Codigo == NumSG, 0, AB2$Codigo)
  
  ##################################################################
  ###### Joint obs, h, and Fores from disagr
  if(j==1){
    Obsh1 <- Resport1A(AB2)
    tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(Obsh1)))
    tvar1 <- "^h"; name1 <- c(grep(tvar1, names(Obsh1)))
    tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(Obsh1)))
    
    colnames(Obsh1)[c(name1,name2)] <- paste(colnames(Obsh1)[c(name1,name2)], 
                                             "_",NomVar[j],sep="") 
  }else{
    Obsh1a <- Resport1A(AB2)
    tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(Obsh1a)))
    tvar1 <- "^h"; name1 <- c(grep(tvar1, names(Obsh1a)))
    tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(Obsh1a)))
    
    colnames(Obsh1a)[c(name1,name2)] <- paste(colnames(Obsh1a)[c(name1,name2)], 
                                              "_",NomVar[j],sep="") 
    
    Obsh1 <- cbind(Obsh1, Obsh1a[,-c(1:2,name0)])
    
  }
  
  
  ##################################################################
  ###### Performance in terms of RMSE
  AC1 <- Resport1(AB2, fechaPand0 = FALSE, fecPand = "Dic_2019",  
                  antes = TRUE, i=hh, filterMod = FALSE, 
                  GuardaPlot = F, fileN)
  
  if(j==1){
    AC1A <- AC1[,c("h","RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1A) <- c("h", 
                        paste(c("RMSEA.0","RMSED.0","RMSED.1"), 
                              NomVar[j],sep = ""))
  }else{
    AC1Aa <- AC1[,c("RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1Aa) <- paste(c("RMSEA.0","RMSED.0","RMSED.1"), 
                             NomVar[j],sep = "")
    AC1A <- cbind(AC1A, AC1Aa)
  }
  
  ##################################################################
  # Performance before Ene2020 (Jan2020)
  
  AC1 <- Resport1(AB2, fechaPand0 = T, fecPand = "Dic_2019",  
                  antes = TRUE, i=hh, filterMod = FALSE, 
                  GuardaPlot = F, fileN)
  
  # Accumulate Performance
  if(j==1){
    AC1B <- AC1[,c("h","RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1B) <- c("h", 
                        paste(c("RMSEA.0","RMSED.0","RMSED.1"), 
                              NomVar[j],sep = ""))
  }else{
    AC1Ba <- AC1[,c("RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1Ba) <- paste(c("RMSEA.0","RMSED.0","RMSED.1"), 
                             NomVar[j],sep = "")
    AC1B <- cbind(AC1B, AC1Ba)
  }
  
  ###################################################################
  # Performance from Ene2020 (Jan2020)
  
  AC1 <- Resport1(AB2, fechaPand0 = T, fecPand = "Dic_2019",  
                  antes = F, i=hh, filterMod = FALSE, 
                  GuardaPlot = F, fileN)
  
  # Accumulate Performance
  if(j==1){
    AC1C <- AC1[,c("h","RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1C) <- c("h", 
                        paste(c("RMSEA.0","RMSED.0","RMSED.1"), 
                              NomVar[j],sep = ""))
  }else{
    AC1Ca <- AC1[,c("RMSEA.0","RMSED.0","RMSED.1")]
    colnames(AC1Ca) <- paste(c("RMSEA.0","RMSED.0","RMSED.1"), 
                             NomVar[j],sep = "")
    AC1C <- cbind(AC1C, AC1Ca)
  }
  
}


### calcula el RMSE relativo
RMSE_rela <- TRUE
if(RMSE_rela){
  ## Exporta resumen de resultados del desempeño desde agregados
  Namestab2 <- c("All_A", "Antes_A", "Desp_A")
  
  if(length(LetraSG)>1) LetraSG <- c("bs")
  tvar1 <- "^RMSEA"; name1 <- c(grep(tvar1, names(AC1A)))
  AC1A.0 <- AC1A[,name1]
  colnames(AC1A.0) <- NomVar
  # print(xtable(AC1A.0), 
  #       file = paste(pathplot, Namestab2[1],LetraSG,".tex",sep = ""), 
  #       compress = FALSE)
  
  tvar1 <- "^RMSEA"; name1 <- c(grep(tvar1, names(AC1B)))
  AC1B.0 <- AC1B[,name1]
  colnames(AC1B.0) <- NomVar
  # print(xtable(AC1B.0),  
  #       file = paste(pathplot, Namestab2[2],LetraSG,".tex",sep = ""), 
  #       compress = FALSE)
  
  tvar1 <- "^RMSEA"; name1 <- c(grep(tvar1, names(AC1C)))
  AC1C.0 <- AC1C[,name1]
  colnames(AC1C.0) <- NomVar
  # print(xtable(AC1C.0), 
  #       file = paste(pathplot, Namestab2[3],LetraSG,".tex",sep = ""), 
  #       compress = FALSE)
  
  
  ## Exporta resumen de resultados del desempeño agregando
  Namestab2 <- c("All_D", "Antes_D", "Desp_D")
  
  tvar1 <- "^RMSED.1"; name1 <- c(grep(tvar1, names(AC1A)))
  AC1A.0D <- AC1A[,name1]
  colnames(AC1A.0D) <- NomVar
  # print(xtable(AC1A.0D), 
  #       file = paste(pathplot, Namestab2[1],LetraSG,".tex",sep = ""), 
  #       compress = FALSE)
  
  tvar1 <- "^RMSED.1"; name1 <- c(grep(tvar1, names(AC1B)))
  AC1B.0D <- AC1B[,name1]
  colnames(AC1B.0D) <- NomVar
  # print(xtable(AC1B.0D),  
  #       file = paste(pathplot, Namestab2[2],LetraSG,".tex",sep = ""), 
  #       compress = FALSE)
  
  tvar1 <- "^RMSED.1"; name1 <- c(grep(tvar1, names(AC1C)))
  AC1C.0D <- AC1C[,name1]
  colnames(AC1C.0D) <- NomVar
  # print(xtable(AC1C.0D),  
  #       file = paste(pathplot, Namestab2[3],LetraSG,".tex",sep = ""), 
  #       compress = FALSE)
}

########################################################
#    Functions to evaluation
### calcula el DM test

Eval_Agre <- TRUE
if(Eval_Agre){
  
  Eval_RMSE.DM <- function(AC1A.0, AC1A.0D, h, Obsh1, NomVar,
                           typeFores = "_A",
                           saveTab = TRUE, path = pathplot, 
                           NomSale = Namestab2[3],
                           ListLags = c(1,3,6,9,12),
                           Samplef=FALSE, antes=FALSE,fecPand = "Dic_2019"){
    
    # AC1A.0 is the RMSE table for the models 
    # ref_i is the reference position in AC1A.0 (RMSE table 
    #        of the models by column), this order depend on file.read object   
    # Ex:ref_i <- 1
    # Obsh1 is the data frame with the observations and the forecast
    # NomVar is the list of models
    # typeFores = "_A" if we create agregate table typeFores = "_D" otherwise
    
    if(typeFores == "_A"){
      AC1A.0A <- AC1A.0
      mods <- 1:ncol(AC1A.0A)
      # mods <- mods[-ref_i]
      # for(i in mods) AC1A.0A[,i] <- AC1A.0[,i]/AC1A.0[,ref_i] 
      # AC1A.0A <- AC1A.0A[,-ref_i]
    }else{
      AC1A.0A <- AC1A.0D
      mods <- 1:ncol(AC1A.0A)
      #mods <- mods[-ref_i]
      #for(i in 1:ncol(AC1A.0A)) AC1A.0A[,i] <- AC1A.0D[,i]/AC1A.0[,ref_i]
    }
    ## Calculate DM test, potentially we could use Martinez et.al. test
    
    #names(Obsh1)
    #head(Obsh1)
    #Obsh <- Obsh1
    if(Samplef){
      fechaPand <- min(Obsh1[Obsh1$fechaf == fecPand,"ordenf"])
      if(antes){
        Obsh1 <- Obsh1 %>% dplyr::filter(ordenf <= fechaPand)
      }else Obsh1 <- Obsh1 %>% dplyr::filter(ordenf > fechaPand)
    }
    
    tvar0 <- paste("*","Obs",sep=""); name0 <- c(grep(tvar0, names(Obsh1)))
    A00 <- Obsh1[,name0]
    
    # We calculate which model is the one that has minimum 
    #  RMSE to use as reference to compare performance
    
    ref_i <- c()
    for(j in 1:nrow(AC1A.0)) ref_i <- c(ref_i, 
                                        which(AC1A.0[j,] == min(AC1A.0[j,]))[1])
    
    ### DM test comparing with the desagregates of the reference
    PvalDM <- c()
    for(fw in 1:h){ 
      #fw <- 3
      tvar1 <- paste("*",NomVar[ref_i[fw]],"$",sep=""); 
      name1 <- c(grep(tvar1, names(Obsh1)))
      A0 <- Obsh1[,name1]
      PvalDM <- c(PvalDM, 
                  dm.test(na.omit(A00[,c(fw)]-A0[,c(fw)]), 
                          na.omit(A00[,c(fw)]-A0[,c(fw+h)]), 
                          #h=fw, alternative ="greater",
                          h=fw, alternative ="two.sided",
                          #h=fw, alternative ="less",
                          varestimator = c("bartlett"))$p.value)
    }
    DM_comp <- data.frame(PvalDM)
    #Noms0 <- c(paste(NomVar[ref_i],c("D"),sep="_"))
    
    Noms0 <- c()
    #ij in mods
    for(ij in mods){
      #ij = 10
      tvar2 <- paste("*",NomVar[ij],"$",sep=""); 
      name2 <- c(grep(tvar2, names(Obsh1)))
      A1 <- Obsh1[,name2]
      #if(length(name2) > 24) break
      ### DM test comparing agregates
      
      PvalDM <- c()
      for(fw in 1:h){
        # This is to choose the best model per horizon
        tvar1 <- paste("*",NomVar[ref_i[fw]],"$",sep=""); 
        name1 <- c(grep(tvar1, names(Obsh1)))
        A0 <- Obsh1[,name1]
        
        if(mean(A0[,c(fw)] == A1[,c(fw)]) != 1){
          PvalDM <- c(PvalDM, 
                      dm.test(na.omit(A00[,c(fw)]-A0[,c(fw)]), 
                              na.omit(A00[,c(fw)]-A1[,c(fw)]), 
                              #h=fw, alternative ="greater",
                              h=fw, alternative ="two.sided",
                              #h=fw, alternative ="less",
                              varestimator = c("bartlett"))$p.value)
        }else PvalDM <- c(PvalDM, 1) 
      }
      
      DM_comp <- cbind(DM_comp, PvalDM)
      
      ### DM test comparing with the desagregates
      PvalDM <- c()
      for(fw in 1:h) PvalDM <- c(PvalDM, 
                                 dm.test(na.omit(A00[,c(fw)]-A0[,c(fw)]), 
                                         na.omit(A00[,c(fw)]-A1[,c(fw+h)]), 
                                         #h=fw, alternative ="greater",
                                         h=fw, alternative ="two.sided",
                                         #h=fw, alternative ="less",
                                         varestimator = c("bartlett"))$p.value)
      DM_comp <- cbind(DM_comp, PvalDM)
      Noms0 <- c(Noms0, paste(NomVar[ij],c("A","D"),sep="_"))
    }
    
    # This is remove because it is duplicate
    DM_comp <- DM_comp[,-1]
    colnames(DM_comp) <- Noms0
    
    DM_comp2 <- DM_comp
    for(i in 1:nrow(DM_comp))
      for(j in 1:ncol(DM_comp))
        DM_comp2[i,j] <- ifelse(DM_comp2[i,j]>0.1,"",
                                ifelse(DM_comp2[i,j]<=0.1 & DM_comp2[i,j]>0.05,"*",
                                       ifelse(DM_comp2[i,j]<=0.05 & DM_comp2[i,j]>0.01,"**",
                                              #ifelse(DM_comp2[i,j]<=0.01 & DM_comp2[i,j]>0.001,"**",
                                              "***")))#)
    #typeFores <- "_A"
    tvar3 <- paste("*",typeFores,sep=""); 
    name3 <- c(grep(tvar3, names(DM_comp2)))
    DM_comp3 <- DM_comp2[,name3]
    for(i in 1:nrow(DM_comp3))
      for(j in 1:ncol(DM_comp3))
        DM_comp3[i,j] <- paste(round(AC1A.0A[i,j],2),DM_comp2[i,name3[j]],sep="")
    
    colnames(DM_comp3) <- names(AC1A.0A)
    #NomSale <- Namestab2[3]
    #Relativo <- 
    #ListLags <- c(1,3,6,9,12)
    if(saveTab){
      print(xtable(t(DM_comp3[ListLags,])), 
            file = paste(path, NomSale,".tex",sep = ""), 
            compress = FALSE)
    }else return(list(DM_comp = DM_comp3, BestMod_h = ref_i))
  }
  
  JointTabNew <- function(A, B, path, NomSale, 
                          ListLags = c(1,3,6,9,12)){
    A <- t(A[ListLags,]); B <- t(B[ListLags,])
    
    AB <- cbind(A, B)
    #rownames(AB)[ref_i] <- rownames(B)[ref_i]
    print(xtable(AB), 
          file = paste(path, NomSale,".tex",sep = ""), 
          compress = FALSE)
  }
  
  ########################################################
  ## Full evaluation set vs Agregates
  #pathplot <- "C:/Users/wmartiri/DesagregaIPC/RMSE_results2/"
  saveTab <- F
  A11 <- Eval_RMSE.DM(AC1A.0, AC1A.0D,h, Obsh1, NomVar,
                      typeFores = "_A",
                      saveTab, path = pathplot, 
                      NomSale = Namestab1[1])
  A1 <- A11$DM_comp
  BestA1 <- A11$BestMod_h[c(1,3,6,9,12)]
  ## Before pandemi evaluation set vs Agregates
  A22 <- Eval_RMSE.DM(AC1B.0, AC1B.0D, h, Obsh1, NomVar,
                      typeFores = "_A",
                      saveTab, path = pathplot, 
                      NomSale = Namestab1[2],
                      Samplef=T, antes=T,fecPand = "Dic_2019")
  A2 <- A22$DM_comp
  BestA2 <- A22$BestMod_h[c(1,3,6,9,12)]
  ## After pandemi evaluation set vs Agregates
  A33 <- Eval_RMSE.DM(AC1C.0, AC1C.0D,h, Obsh1, NomVar,
                      typeFores = "_A",
                      saveTab, path = pathplot, 
                      NomSale = Namestab1[3],
                      Samplef=T, antes=F,fecPand = "Dic_2019")
  A3 <- A33$DM_comp
  BestA3 <- A33$BestMod_h[c(1,3,6,9,12)]
  #########################################
  ## Full evaluation set vs Desagregates
  B1 <- Eval_RMSE.DM(AC1A.0, AC1A.0D,h, Obsh1, NomVar,
                     typeFores = "_D",
                     saveTab, path = pathplot, 
                     NomSale = Namestab1[1])$DM_comp
  ## Before pandemi evaluation set vs Desagregates
  B2 <- Eval_RMSE.DM(AC1B.0, AC1B.0D,h, Obsh1, NomVar,
                     typeFores = "_D",
                     saveTab, path = pathplot, 
                     NomSale = Namestab1[2],
                     Samplef=T, antes=T,fecPand = "Dic_2019")$DM_comp
  ## After pandemi evaluation set vs Desagregates
  B3 <- Eval_RMSE.DM(AC1C.0, AC1C.0D,h, Obsh1, NomVar,
                     typeFores = "_D",
                     saveTab, path = pathplot, 
                     NomSale = Namestab1[3],
                     Samplef=T, antes=F,fecPand = "Dic_2019")$DM_comp
  
  A1[,which(names(A1) == "PYXIPP")] <- "-"
  A1[,which(names(A1) == "PYFXIPP")] <- "-"
  A2[,which(names(A2) == "PYXIPP")] <- "-"
  A2[,which(names(A2) == "PYFXIPP")] <- "-"
  A3[,which(names(A3) == "PYXIPP")] <- "-"
  A3[,which(names(A3) == "PYFXIPP")] <- "-"
  
  JointTabNew(A1, B1, path= pathplot, 
              NomSale = paste(Namestab2[1],LetraSG,sep ="_"))
  JointTabNew(A2, B2, path= pathplot, 
              NomSale = paste(Namestab2[2],LetraSG,sep ="_"))
  JointTabNew(A3, B3, path= pathplot, 
              NomSale = paste(Namestab2[3],LetraSG,sep ="_"))
  
  
}

NomVar[BestA1]
NomVar[BestA2]
NomVar[BestA3]

