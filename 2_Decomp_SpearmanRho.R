library(glmmTMB)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#######################################################################################################################################################################################################################
##################################################################################################### read in data & custom functions ##################################################################################
#######################################################################################################################################################################################################################

# set working directory to master folder...
setwd("C:/...")

# and paths to folders to data, scripts, etc.
InDir <- "Papers/ClimateDecomposition/Archiving/"
InDirData <- "MountainBirds/ModelData/"

# Path to final models
CYJulyPath <- "Papers/ClimateDecomposition/Archiving/CurrentMayJuly"
PYJulyPath <- "Papers/ClimateDecomposition/Archiving/PreviousMayJuly"
CYJunePath <- "Papers/ClimateDecomposition/Archiving/CurrentMayJune"

############################## Custom functions
source(paste0(InDir, "Funs_DecompStudy.R"))

############################## Species count and covariate data Current: May-July
load(paste0(InDirData, "ModelData_A_CY_3WayDecomp.RData"))

############################## Repeat with: Species count and covariate data Previous: May-July
# load(paste0(InDirData, "ModelData_A_PY_3WayDecomp.RData"))

############################## Repeat with: Species count and covariate data Current: May-June
# load(paste0(InDirData, "ModelData_A_MayJune_3WayDecomp.RData"))

# Means and SDs for standardization
Means_SDs_MayJuly <- read.csv(paste0(InDirData, "Means_SDs_CovarDataSummer_SC.csv"))
Means_SDs_MayJune <- read.csv(paste0(InDirData, "Means_SDs_CovarDataSummer_SC_MayJune.csv"))

# Cross-validatin results
CV <- read.csv(paste0(InDir, "Results/CVCYJuly.csv"))

#######################################################################################################################################################################################################################
###################################################################### select the focal species ##################################################################################
#######################################################################################################################################################################################################################

CVSel <- CV[CV$Pears_Ind_All >= 0.5 & CV$Pears_Ind_No >= 0.4 & CV$Pears_Ind_Swe >= 0.4 & CV$Pears_Ind_Fi >= 0.4, ]

###############################################################################################################################################################################
################################################################ create Covar data for all routes #######################################################################
###############################################################################################################################################################################

SketchCovarsCYJuly <- FunCreateSketchCovars(InCovarData = CovarDataSummer$SC)
SketchCovarsPYJuly <- FunCreateSketchCovars(InCovarData = CovarDataSummer$SP)
SketchCovarsCYJune <- FunCreateSketchCovars(InCovarData = CovarDataSummer$SC)

###############################################################################################################################################################################
############################################################## calculate Spearman rho between climate components ##############################################################
###############################################################################################################################################################################

# CovarData are the Sketch data created above. This is used to create predicted values.
# Range data is the full set of covariate data. This is used to obtain the covariate ranges

CYJulyModels <- list.files(CYJulyPath, pattern = "MAM")
RhosCYJuly <- FunCalcSpearmanRhos(ModelPath = CYJulyPath, ModelList = CYJulyModels, CovarData = SketchCovarsCYJuly, RangeData = CovarDataSummer$SC, StandMeanSD = Means_SDs_MayJuly)
write.csv(RhosCYJuly, "SpearmanRhos_CYMayJuly.csv", row.names = F)

PYJulyModels <- list.files(PYJulyPath, pattern = "MAM")
RhosPYJuly <- FunCalcSpearmanRhos(ModelPath = PYJulyPath, ModelList = PYJulyModels, CovarData = SketchCovarsPYJuly, RangeData = CovarDataSummer$SP, StandMeanSD = Means_SDs_MayJuly)
write.csv(RhosPYJuly, "SpearmanRhos_PYMayJuly.csv", row.names = F)

CYJuneModels <- list.files(CYJunePath, pattern = "MAM")
RhosCYJune <- FunCalcSpearmanRhos(ModelPath = CYJunePath, ModelList = CYJuneModels, CovarData = SketchCovarsCYJune, RangeData = CovarDataSummer$SC, StandMeanSD = Means_SDs_MayJune)
write.csv(RhosCYJune, "SpearmanRhos_CYMayJune.csv", row.names = F)

##############################################################################################################################################################################
######################################################################## Raster for predicted abundance ######################################################################
##############################################################################################################################################################################

CYJulyModels <- list.files(CYJulyPath, pattern = "MAM")
FunSketchModelsSpatial(ModelPath = CYJulyPath, ModelList = CYJulyModels, CovarData = SketchCovarsCYJuly, RangeData = CovarDataSummer$SC, StandMeanSD = Means_SDs_MayJuly)
FunSketchModelsTemporal(ModelPath = CYJulyPath, ModelList = CYJulyModels, CovarData = SketchCovarsCYJuly, RangeData = CovarDataSummer$SC, StandMeanSD = Means_SDs_MayJuly)
FunSketchModelsResidual(ModelPath = CYJulyPath, ModelList = CYJulyModels, CovarData = SketchCovarsCYJuly, RangeData = CovarDataSummer$SC, StandMeanSD = Means_SDs_MayJuly)





