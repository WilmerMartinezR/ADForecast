
# VERSION: 1.0

# ########################################################################
#   In this code defines the object path_permaNent, line 47
#   
# ########################################################################


# ########################################################################
#                           Libraries
# ########################################################################

# List of libraries to use
load.lib <- c("readxl", "multDM", "gtools", "plyr", "mvnTest", "rlist",
              "tsDyn", "MTS", "vars", "forecast", "portes", "urca",
              "imputeTS", "dplyr", "bvarsv", "tidyverse", "aTSA", 
              "data.table", "tseries", "mFilter", "DescTools", "svars",
              "ggplot2", "moments", "mlVAR", "jcolors", "tsDyn", "sarima",
              "sysid","lawstat","FinTS", "xtable", "glmnet","randomForest",
              "astsa")

# Identify packages to install
install.lib <- load.lib[!load.lib %in% installed.packages()]

# Imprimir los nombres de los paquetes no instalados o confirmar que todos están instalados
if (length(install.lib) > 0) {
  
  cat("These packages need to be installed:\n")
  cat(install.lib, sep = "\n")
  
  # Instalar los paquetes no instalados
  install.packages(install.lib)
} else {
  cat("All packages are installed.\n")
}

# Load packages
sapply(load.lib, require, character = TRUE)

# ########################################################################
#                           Files path
# ########################################################################

# Name of the permanent path
# Define the path of the main folder
path_permaNent <- "/Users/..."

# Name of directory to be created
Name_Newdir <- "Results_Pronos/"
# To create folders to store results
dir.create(paste(path_permaNent,Name_Newdir,sep = ""))

### Save results
Name_Newdir2 <- "Results_Results/"
dir.create(paste(path_permaNent,Name_Newdir2,sep = ""))

path_Datos <- paste(path_permaNent,"Data/",sep = "")
path_Code <- paste(path_permaNent,"Code/",sep = "")

## The directory with the forecast output
path.out <- paste(path_permaNent,Name_Newdir,sep = "")
## The directory with the results (summary table)
path.Results <- paste(path_permaNent,Name_Newdir2,sep = "")
pathplot <- path.Results

## The directory with the plot results

Name_Newdir3 <- "PlotsDesag/"
dir.create(paste(path_permaNent,Name_Newdir3,sep = ""))
pathplotStore <- paste(path_permaNent,Name_Newdir3,sep = "")


# ########################################################################
#                           Data files
# ########################################################################

# Define
Country <- "UK" # UK for the United Kingdom, US for the United States and CO for Colombia

DataFile <- "BASEUKFINALDEF_VW.xlsx"
SheetFile <- "Hoja1"


FechaIniSample <- "Ene_2010"
FechaFinSample <- "Jun_2023"
source(paste(path_Code, "Funciones_insumo.R",sep=""))
source(paste(path_Code, "GB0.R",sep=""))
source(paste(path_Code, "Models.R",sep=""))

# ########################################################################
# Define criterion PCA
# ########################################################################

CritACP <- 0.7
source(paste(path_Code, "PCA_Yt.R",sep=""))

# ########################################################################
# Define forecast horizon and initial sample size
# ########################################################################

h <- 12 # horizon forecast
W.Wide <- 60 # initial sample size

# ########################################################################
# Run models and forecast.
# Next lines try to run in separate R sessions or in parallel.
# ########################################################################

source(paste(path_Code, "M0.R",sep=""))
source(paste(path_Code, "ARpDirect.R",sep=""))
source(paste(path_Code, "RidgeYt.R",sep=""))
source(paste(path_Code, "LassoYt.R",sep=""))
source(paste(path_Code, "RFYt.R",sep=""))
