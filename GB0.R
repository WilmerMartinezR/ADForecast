

GB_0 <- read_excel(paste(path_Datos,DataFile,sep=""),sheet = SheetFile)

# Define months names
months <- c("Ene","Feb","Mar","Abr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dic")
#months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dec")
first.year <- as.numeric(substr(FechaIniSample,5,8))
last.year <- as.numeric(substr(FechaFinSample,5,8)) 
aos <- rep(first.year:last.year,each=12)
Maos <- paste(months,aos,sep = "_")
lastObs <- which(substr(FechaFinSample,1,3) == months) 
names(GB_0) <- c("Codigo",names(GB_0)[2:6],Maos[1:(length(Maos)-lastObs)])

# Choose a subset or a sample period
fechaIni0 <- which(names(GB_0) == FechaIniSample) 

GB1 <- data.frame(GB_0[,c(1:6,fechaIni0:ncol(GB_0))])


#############################################################################
##               Calculate year-on-year inflation
#############################################################################
h <- 12
VAtransf <- TRUE 
if(VAtransf){
  orden0 <- 12
  GB1.VA <- GB1
  
  for(i in 1:nrow(GB1)){
    #j=h+1;i=1
    GB2 <- t(GB1[i,-c(1:(fechaIni0-1)) ])
    iNA <- which(is.na(GB2) )
    if(length(iNA) > 0) GB2[iNA] <- mean(GB2,na.rm=T) 
    GB1.VA[i,(fechaIni0+h):ncol(GB1)] <- t(vAnual(GB2, orden0))
  }
  GB1.VA <- GB1.VA[,-c(fechaIni0:(fechaIni0+h-1))]
}

if(fechaIni0 > 6) GB1.VA <- GB1.VA[,-c(6:(fechaIni0-1))]
GB1.VA <- data.frame(GB1.VA)

