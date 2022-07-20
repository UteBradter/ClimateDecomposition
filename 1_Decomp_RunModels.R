library(glmmTMB)

#######################################################################################################################################################################################################################
##################################################################################################### read in data, custom functions, set paths #######################################################################
#######################################################################################################################################################################################################################

# set working directory to master folder...
setwd("C:/...")

# and paths to folders to data, scripts, etc.
InDir <- "Papers/ClimateDecomposition/Archiving/"
InDirData <- "MountainBirds/ModelData/"

# set paths to store model results
CurPathCYJuly <- "Papers/ClimateDecomposition/Archiving/CurrentMayJuly"
CurPathPYJuly <- "Papers/ClimateDecomposition/Archiving/PreviousMayJuly"
CurPathCYJune <- "Papers/ClimateDecomposition/Archiving/CurrentMayJune"

############################## Custom functions
source(paste0(InDir, "Funs_DecompStudy.R"))

############################## Species count and covariate data Current: May-July
load(paste0(InDirData, "ModelData_A_CY_3WayDecomp.RData"))

############################## Repeat with: Species count and covariate data Previous: May-July
# load(paste0(InDirData, "ModelData_A_PY_3WayDecomp.RData"))

############################## Repeat with: Species count and covariate data Current: May-June
# load(paste0(InDirData, "ModelData_A_MayJune_3WayDecomp.RData"))

############################## Species names, Euring code, Covariate groups, etc.
SpeciesInfo <- read.csv(paste0(InDir, "SpeciesInfo.csv"), stringsAsFactors = F)

############################## model family and zero-inflation
FamilyInfo <- read.csv(paste0(InDir, "FamilyInfo.csv"))

#######################################################################################################################################################################################################################
################################################################################################### set up data #######################################################################################################
#######################################################################################################################################################################################################################

# sort by Euring code
SpeciesInfo <- SpeciesInfo[order(SpeciesInfo$Euring), ]
SpeciesInfo <- SpeciesInfo[!duplicated(SpeciesInfo$Euring), ]

# counts of the 39 species
SFNInds <- SFNInds[SFNInds$euring %in% SpeciesInfo$Euring, ]	

# Dipper separately (convergence problems in some models with neg bin; repeat with Poisson for Dipper only)
SpInfoDipper <- SpeciesInfo[SpeciesInfo$englishname == "Dipper", ]
FamInfoDipper <- FamilyInfo
FamInfoDipper[FamInfoDipper$Species == "Dipper", ]$Family <- "poisson"

# the 33 species
SpInfo33 <- SpeciesInfo[SpeciesInfo$englishname != "Mallard" & SpeciesInfo$englishname != "Common_gull" & SpeciesInfo$englishname != "Raven" & 
SpeciesInfo$englishname != "Common_redpoll" & SpeciesInfo$englishname != "Fieldfare" & SpeciesInfo$englishname != "Willow_warbler" , ]

#######################################################################################################################################################################################################################
##################################################################################################### Model formulas ##################################################################################################
#######################################################################################################################################################################################################################

CountForm <- list(

	S_M_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MA_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MDA_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Decid.for + Agri + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MvDA_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Decid.for + Agri + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MvFA_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Decid.for + Other.for + Agri + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MvD_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Decid.for + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MAA2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + I(Agri^2) + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MARad_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MARadWe2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + Wetland + I(Wetland^2) + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MAWa2A2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + I(Agri^2) + Wetland + Water + I(Water^2) + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MARadWa2A2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + I(Agri^2) + Wetland + Water + I(Water^2) + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MAWa2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + Wetland + Water + I(Water^2) + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MAWe2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + Wetland + I(Wetland^2) + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MAWa2We2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + Wetland + I(Wetland^2) + Water + I(Water^2) + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MAWa2We2A2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Agri + I(Agri^2) + Wetland + I(Wetland^2) + Water + I(Water^2) + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MDARadA2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Decid.for + Agri + I(Agri^2) + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MDARadWe2A2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Fjall.open + Decid.for + Agri + I(Agri^2) + Wetland + I(Wetland^2) + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MvDARad_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Decid.for + Agri + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MvDARadWe2A2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + Rad + yrStand + I(yrStand^2) + Fjall.veg + Decid.for + Agri + I(Agri^2) + Wetland + I(Wetland^2) + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route),

	S_MvFAA2_Form = Count ~ ST_XY * SR_XY + ST_T * SR_T + ST_RXYT * SR_RXYT + Sl + yrStand + I(yrStand^2) + Fjall.veg + Decid.for + Other.for + Agri + I(Agri^2) + Wetland + Water + 
	offset(log(Effort)) + Survey + Unit + (1|Route)
)

#######################################################################################################################################################################################################################
######################################################################################################### run the models ##############################################################################################
#######################################################################################################################################################################################################################

###################################### Current May-July: Output saved to CurPathCYJuly

OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpeciesInfo, CovarData = CovarDataSummer, Season = "Summer", Year = "Current", CovariateForm = CountForm, 
FamilyBySpecies = FamilyInfo, VarSel = T, PdfName = "CurMayJuly")

# Dipper with Poisson
OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpInfoDipper, CovarData = CovarDataSummer, Season = "Summer", Year = "Current", CovariateForm = CountForm, 
FamilyBySpecies = FamInfoDipper, VarSel = T, PdfName = "CurMayJuly")

####################################### Previous May-July: Output saved to CurPathPYJuly

OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpInfo33, CovarData = CovarDataSummer, Season = "Summer", Year = "Previous", CovariateForm = CountForm, 
FamilyBySpecies = FamilyInfo, VarSel = T, PdfName = "PrevMayJuly")

# Dipper with Poisson
OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpInfoDipper, CovarData = CovarDataSummer, Season = "Summer", Year = "Previous", CovariateForm = CountForm, 
FamilyBySpecies = FamInfoDipper, VarSel = T, PdfName = "PrevMayJuly")

# Long-tailed skua without zero-inflation
SpInfoSkua <- SpeciesInfo[SpeciesInfo$englishname == "Long-tailed_skua", ]
FamInfoSkua <- FamilyInfo
FamInfoSkua[FamInfoDipper$Species == "Long-tailed_skua", ]$ZeroInfl <- "No"
OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpInfoSkua, CovarData = CovarDataSummer, Season = "Summer", Year = "Previous", CovariateForm = CountForm, 
FamilyBySpecies = FamInfoSkua, VarSel = T, PdfName = "PrevMayJuly")

######################################## Current May-June: Output saved to CurPathCYJune 

OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpInfo33, CovarData = CovarDataSummer, Season = "Summer", Year = "Current", CovariateForm = CountForm, 
FamilyBySpecies = FamilyInfo, VarSel = T, PdfName = "CurMayJune")

# Dipper with Poisson:
OutSummary <- FunGLMMs(SpeciesData = SFNInds, SpeciesInfo = SpInfoDipper, CovarData = CovarDataSummer, Season = "Summer", Year = "Current", CovariateForm = CountForm, 
FamilyBySpecies = FamInfoDipper, VarSel = T, PdfName = "CurMayJune")

#######################################################################################################################################################################################################################
################################################################################### extract regression coefficients and standard errors ###############################################################################
#######################################################################################################################################################################################################################

# Current May-July

MAMModelsCYJuly <- list.files(CurPathCYJuly, pattern = "MAM")
BetasMAMCYJuly <- FunExtractBetas(Path = CurPath, MyModelList = MAMModelsCYJuly, Type = "MAM")

write.csv(BetasMAMCYJuly$Betas, "BetasMAMCYJuly.csv", row.names = F)
write.csv(BetasMAMCYJuly$StdErrors, "SEsMAMCYJuly.csv", row.names = F)

FullModelsCYJuly <- list.files(CurPathCYJuly, pattern = "Full")
BetasFullCYJuly <- FunExtractBetas(Path = CurPath, MyModelList = FullModelsCYJuly, Type = "Full")

write.csv(BetasFullCYJuly$Betas, "BetasFullCYJuly.csv", row.names = F)
write.csv(BetasFullCYJuly$StdErrors, "SEsFullCYJuly.csv", row.names = F)

# Previous May-July

MAMModelsPYJuly <- list.files(CurPathPYJuly, pattern = "MAM")
BetasMAMPYJuly <- FunExtractBetas(Path = CurPath, MyModelList = MAMModelsPYJuly, Type = "MAM")

write.csv(BetasMAMPYJuly$Betas, "BetasMAMPYJuly.csv", row.names = F)
write.csv(BetasMAMPYJuly$StdErrors, "SEsMAMPYJuly.csv", row.names = F)

FullModelsPYJuly <- list.files(CurPathPYJuly, pattern = "Full")
BetasFullPYJuly <- FunExtractBetas(Path = CurPath, MyModelList = FullModelsPYJuly, Type = "Full")

write.csv(BetasFullPYJuly$Betas, "BetasFullPYJuly.csv", row.names = F)
write.csv(BetasFullPYJuly$StdErrors, "SEsFullPYJuly.csv", row.names = F)

# Current May-June

MAMModelsCYJune <- list.files(CurPathCYJune, pattern = "MAM")
BetasMAMCYJune <- FunExtractBetas(Path = CurPath, MyModelList = MAMModelsCYJune, Type = "MAM")

write.csv(BetasMAMCYJune$Betas, "BetasMAMCYJune.csv", row.names = F)
write.csv(BetasMAMCYJune$StdErrors, "SEsMAMCYJune.csv", row.names = F)

FullModelsCYJune <- list.files(CurPathCYJune, pattern = "Full")
BetasFullCYJune <- FunExtractBetas(Path = CurPath, MyModelList = FullModelsCYJune, Type = "Full")

write.csv(BetasFullCYJune$Betas, "BetasFullCYJune.csv", row.names = F)
write.csv(BetasFullCYJune$StdErrors, "SEsFullCYJune.csv", row.names = F)

#######################################################################################################################################################################################################################
############################################################################################### cross-validation ##################################################################################################
#######################################################################################################################################################################################################################

CYJulyModels <- list.files(CurPathCYJuly, pattern = "MAM")
CVResults <- FunCV(ModelPath = CurPathCYJuly, ModelList = CYJulyModels, SpeciesData = SFNInds, SpeciesInfo = SpeciesInfo, CovarData = CovarDataSummer, Season = "Summer", Year = "Current")
write.csv(CVResults, "CVCYJuly.csv", row.names = F)



