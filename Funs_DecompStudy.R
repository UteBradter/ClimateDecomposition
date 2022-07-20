library(glmmTMB)

##########################################################################################################################################################################################################
#################################################################################### Standardize covariates ##################################################################################
##########################################################################################################################################################################################################

FunStand <- function(Covars){
# divide out the coordinates, year and effort to retain for residual diagnostic plots, etc.
	SpaceTimeEffort <- Covars[, c("x", "y", "yr", "Effort")]
# identify columns of covariates that will not be standardized (Route, Survey, Unit
	SetFalse <- which(colnames(Covars) %in% c("Route", "Survey", "Unit") == T)
# standardize 
	Stand <- apply(Covars[, -c(SetFalse)], 2, FUN = function(xx){(xx - mean(xx)) / sd(xx)})
# find indices for x, y, yr, Effort
	Inds <- c(grep("yr", colnames(Stand)), grep("x", colnames(Stand)), grep("y", colnames(Stand)), grep("Effort", colnames(Stand)))
	Inds <- unique(Inds)
	colnames(Stand)[Inds] <- paste0(colnames(Stand)[Inds], "Stand")
# add everything together
	Covars <- cbind(Covars[, SetFalse], Stand, SpaceTimeEffort)
	return(Covars)
}

##########################################################################################################################################################################################################
#################################################################################### Run models with tryCatch ##################################################################################
##########################################################################################################################################################################################################

FunRunModel <- function(MF, ZMF, CMD, Fam){
	m1 <- NULL
	m1 <- tryCatch({glmmTMB(MF, ziformula = ZMF, data = CMD, family = Fam)},
	error = function(cond){
		message(MF)
		message(ZMF)
		message(cond)
		return(NULL)
	},
	warning = function(cond) {
		message(MF)
		message(ZMF)
		message(cond)
		return(m1)
	})
	return(m1)
}

##########################################################################################################################################################################################################
#################################################################################### At least one element is not NA ##################################################################################
##########################################################################################################################################################################################################

FunNotAllIsNA <- function(x){
	any(!is.na(x)) == T
}

##########################################################################################################################################################################################################
############################################################################################## GLMMs, variable selection #################################################################################
##########################################################################################################################################################################################################

# uses package glmmTMB
# error distribution and possibility of zero-inflation is supplied via function arguments
# writes out the full model
# if the argument for variable selection is set to TRUE, writes out a table with the AIC values for each submodel
# SpeciesInfo contains the English and latin names of the focal species, the EURING code and the habitat formulae used to select the covariate form
# FamilyBySpecies contains information on the error structure used and on zero-inflation
# Family and ZeroInfl = NULL were used for the initial comparison between different error structures and whether or not to include zero-inflation. 
# If FamilyBySpecies supplied, these arguments are no longer of consequence

FunGLMMs <- function(SpeciesData, SpeciesInfo, CovarData, Year, Season = Null, CovariateForm, Family, ZeroInfl = NULL, FamilyBySpecies = NULL, VarSel = F, PdfName){

# set up the species info data			
	CurSpInfo <- SpeciesInfo[!duplicated(SpeciesInfo$Euring), ]
# for each species...		
	for (i in 1 : length(CurSpInfo$Euring)){
		print(CurSpInfo$englishname[i])
		ModelListFull <- NULL
		ModelListMAM <- NULL
# set up an empty data frame for summaries (can be enlarged to hold results from model validation)
		TestSummary <- as.data.frame(matrix(nrow = 1, ncol = 2))
		colnames(TestSummary) <- c("Species", "AIC")
# extract species counts
		CurSpecData <- SpeciesData[SpeciesData$euring == CurSpInfo$Euring[i], ]
# extract species information
		CurSpecInfo <- CurSpInfo[CurSpInfo$Euring == CurSpInfo$Euring[i], ]	
		TestSummary$Species <- as.character(CurSpecInfo$Vetenskapligt.namn)
# identify the covariate data
		if(is.null(Season)){
			SeasonShort <- ifelse(CurSpecInfo$Fjall_skog_area == "Summer", "S", "Y")
		}else{
			if(!Season %in% c("Summer", "Year"))stop("Argument to Season must be either Summer or Year")
			SeasonShort <- substr(Season, 1, 1)
		}
		CurCovData <- paste0(SeasonShort, ifelse(Year == "Current", "C", "P"))
		CurCovData <- CovarData[[which(names(CovarData) == CurCovData)]]
# standardize
		CurCovData <- FunStand(CurCovData)
# combine species and covariate data
		CurModelData <- merge(CurSpecData[, c("Route", "yr", "Count", "Survey")], CurCovData, by = c("Route", "yr", "Survey"), all.y = T)
# zero-fill
		CurModelData[is.na(CurModelData$Count), ]$Count <- 0
# find the relevant model formulae
		ModelFormula <- CovariateForm[[which(names(CovariateForm) == paste(SeasonShort, CurSpecInfo$HabFinal, "Form", sep = "_"))]]
# identify the family, if not supplied as argument
		if(!is.null(FamilyBySpecies)){
			CurFamily <- FamilyBySpecies[FamilyBySpecies$Euring == CurSpecInfo$Euring, ]
			CurFamily <- CurFamily[!duplicated(CurFamily$Euring), ]
			if("Family" %in% colnames(CurFamily)){
				Family <- CurFamily$Family
				if(CurFamily$ZeroInfl == "Yes"){ZeroInfl <- "Zero-inflation"}else{ZeroInfl <- NULL}
			}else{
				Family <- CurFamily$FinalMin
				if(substr(Family, 1, 2) == "Zi"){ZeroInfl <- "Zero-inflation"}else{ZeroInfl <- NULL}
				if(Family == "Pois" | Family == "ZiPois"){Family <- "poisson"}
				if(Family == "NB2" | Family == "ZiNB2"){Family <- "nbinom2"}
			}
		}
# create the model formula for the zero-inflation part (if there is one)
		if(!is.null(ZeroInfl)){
# remove random effect and some fixed effects from binomial part of model (I remove Survey and Unit because neither of them should have an effect on the probability of a false zero; Effort may)
			ZiModelFormula <- update(ModelFormula, ~. - Unit -Survey -(1|Route))
# removing offset is not supported by update, so have to do it the long way...
# ...remove the response variable...			
			ZiModelFormula <- as.character(ZiModelFormula)[3]
# ...remove the offset...
			ZiModelFormula <- unlist(strsplit(ZiModelFormula, split = "\\+"))
			ZiModelFormula <- ZiModelFormula[- grep("offset", ZiModelFormula)]
# ... back to formula...
			ZiModelFormula <- as.formula(paste("~", paste(ZiModelFormula, collapse = "+"), sep = " "))
# ... add effort as covariate (here not as offset, which would equal the number of trials)
			ZiModelFormula <- update(ZiModelFormula, ~. + EffortStand)
		}else{
			ZiModelFormula <- as.formula("~0")
		}
# run the models, catching errors and warnings, if this argument was chosen
		m1 <- NULL
		m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
# save model
		ModelOutName <- gsub(" ", "_", CurSpecInfo$englishname)
		ModelListFull <- append(ModelListFull, list(NextModel = m1))
		names(ModelListFull)[length(names(ModelListFull))] <- paste("M", ModelOutName, sep = "_")
		save(list = "ModelListFull", file = paste0(PdfName, "_Full_", CurSpecInfo$englishname, ".RData"))
# optionally, variable selection
		if(VarSel == T){
# set up empty data frame for VarSel AICs; this will be read out for each species for safety, as VarSel takes long for each species
			AICs <- as.data.frame(matrix(nrow = 1, ncol = 33))
			colnames(AICs) <- c("Species", "Habitat", "Family", "ZeroInfl", "Full", "yrStand2", "ST_RXYTxSR_RXYT", "ST_TxSR_T", "ST_XYxSR_XY", "yrStand", "ST_RXYT", "SR_RXYT", "ST_RXYT_SR_RXYT", "ST_T", "SR_T", "ST_T_SR_T",
 "ST_XY", "SR_XY", "ST_XY_SR_XY", "Agri2", "Water2", "Wetland2", "Rad", "Sl", "Agri", "Water", "Wetland", "Other.for", "Decid.for", "Other_Decid.for", "Fjall.open", "Fjall.veg", "Fjall.open_veg")	
			AICs$Species <- as.character(CurSpecInfo$Vetenskapligt.namn)
			AICs$Habitat <- as.character(CurSpecInfo$HabFinal)
			AICs$Family <- as.character(m1$modelInfo$family)[1]
# VarSel: backwards; in Zi-model first; throughout the process the current best model gets assigned as mMin, and the minimum AIC is extracted after each step
# full model		
			AICs$Full <- try(AIC(m1))
			try(mMin <- m1)
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			if(!is.null(ZeroInfl)){
				ZiAICs <- AICs[, -(1:5)]
				colnames(ZiAICs) <- paste("Zi", colnames(ZiAICs), sep = "_")
				AICs <- cbind(AICs, ZiAICs)
			}

# remove quadratics if present...
# ... for Agri...
			if(!is.null(ZeroInfl) & "I(Agri^2)" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - I(Agri^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Agri2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Agri2 > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + I(Agri^2))
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("I(Agri^2)" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. - I(Agri^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Agri2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Agri2 > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + I(Agri^2))
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# ... for Water...
			if(!is.null(ZeroInfl) & "I(Water^2)" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - I(Water^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Water2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Water2 > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + I(Water^2))
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("I(Water^2)" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. - I(Water^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Water2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Water2 > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + I(Water^2))
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# ... for Wetland...
			if(!is.null(ZeroInfl) & "I(Wetland^2)" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - I(Wetland^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family) 
				if(!is.null(m1)){
					AICs$Zi_Wetland2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Wetland2 > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + I(Wetland^2))
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("I(Wetland^2)" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. - I(Wetland^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Wetland2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Wetland2 > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + I(Wetland^2))
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# remove the quadratic for year...	
			if(!is.null(ZeroInfl)){
				ZiModelFormula <- update(ZiModelFormula, ~. -I(yrStand^2))
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_yrStand2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_yrStand2 > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +I(yrStand^2))
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
		
			ModelFormula <- update(ModelFormula, ~. -I(yrStand^2))
			m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
			if(!is.null(m1)){
				AICs$yrStand2 <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(AICs$yrStand2 > MinAIC){
					ModelFormula <- update(ModelFormula, ~. +I(yrStand^2))
				}else{
					mMin <- m1
				}
			}			
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)

# remove the residual climate interaction...
			if(!is.null(ZeroInfl)){
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_RXYT:SR_RXYT)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_ST_RXYTxSR_RXYT <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced model becomes the best model.
					if(AICs$Zi_ST_RXYTxSR_RXYT > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +ST_RXYT:SR_RXYT)
					}else{
						mMin <- m1
					}			
					MinAIC <- min(AICs[, -(1:4)], na.rm = T)
				}
			}

			ModelFormula <- update(ModelFormula, ~. -ST_RXYT:SR_RXYT)
			m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
			if(!is.null(m1)){
				AICs$ST_RXYTxSR_RXYT <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(AICs$ST_RXYTxSR_RXYT > MinAIC){
					ModelFormula <- update(ModelFormula, ~. +ST_RXYT:SR_RXYT)
				}else{
					mMin <- m1
				}
			}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)

# remove the temporal average interaction...
			if(!is.null(ZeroInfl)){
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_T:SR_T)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_ST_TxSR_T <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_ST_TxSR_T > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +ST_T:SR_T)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			ModelFormula <- update(ModelFormula, ~. -ST_T:SR_T)
			m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
			if(!is.null(m1)){
				AICs$ST_TxSR_T <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(AICs$ST_TxSR_T > MinAIC){
					ModelFormula <- update(ModelFormula, ~. +ST_T:SR_T)
				}else{
					mMin <- m1
				}
			}			
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)

# remove the spatial average interaction...
			if(!is.null(ZeroInfl)){
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_XY:SR_XY)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_ST_XYxSR_XY <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_ST_XYxSR_XY > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +ST_XY:SR_XY)
					}else{
						mMin <- m1
					}	
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			ModelFormula <- update(ModelFormula, ~. -ST_XY:SR_XY)
			m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
			if(!is.null(m1)){
				AICs$ST_XYxSR_XY <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(AICs$ST_XYxSR_XY > MinAIC){
					ModelFormula <- update(ModelFormula, ~. +ST_XY:SR_XY)
				}else{
					mMin <- m1
				}
			}			
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)

# if the quadratic year term has been removed, remove the year term...
			if(!is.null(ZeroInfl) & !"I(yrStand^2)" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. -yrStand)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_yrStand <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_yrStand > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +yrStand)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if(!"I(yrStand^2)" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. -yrStand)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$yrStand <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$yrStand > MinAIC){
						ModelFormula <- update(ModelFormula, ~. +yrStand)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# if the residual interaction has been removed, remove both of the residual terms separately and together...
			if(!is.null(ZeroInfl) & !"ST_RXYT:SR_RXYT" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_RXYT)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Zi_ST_RXYT <- try(AIC(m1))}
				ZiModelFormula <- update(ZiModelFormula, ~. +ST_RXYT -SR_RXYT)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Zi_SR_RXYT <- try(AIC(m2))}
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_RXYT -SR_RXYT)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Zi_ST_RXYT_SR_RXYT <- try(AIC(m3))}

				AICSub <- AICs[, c("Zi_ST_RXYT", "Zi_SR_RXYT", "Zi_ST_RXYT_SR_RXYT")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub)){
  					ZiModelFormula <- update(ZiModelFormula, ~. +ST_RXYT +SR_RXYT)
				}else{
					if(MinSub > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +ST_RXYT +SR_RXYT)
					}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
						MinTerm <- names(which.min(AICSub))
						if(MinTerm == "Zi_ST_RXYT"){
							ZiModelFormula <- update(ZiModelFormula, ~. +SR_RXYT)
							mMin <- m1
						}
							if(MinTerm == "Zi_SR_RXYT"){
							ZiModelFormula <- update(ZiModelFormula, ~. +ST_RXYT)
							mMin <- m2
						}
						if(MinTerm == "Zi_ST_RXYT_SR_RXYT"){mMin <- m3}
					}
				}
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if(!"ST_RXYT:SR_RXYT" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. -ST_RXYT)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$ST_RXYT <- try(AIC(m1))}
				ModelFormula <- update(ModelFormula, ~. +ST_RXYT -SR_RXYT)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$SR_RXYT <- try(AIC(m2))}
				ModelFormula <- update(ModelFormula, ~. -ST_RXYT -SR_RXYT)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$ST_RXYT_SR_RXYT <- try(AIC(m3))}

				AICSub <- AICs[, c("ST_RXYT", "SR_RXYT", "ST_RXYT_SR_RXYT")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub)){
  					ModelFormula <- update(ModelFormula, ~. +ST_RXYT +SR_RXYT)
				}else{
					if(MinSub > MinAIC){
						ModelFormula <- update(ModelFormula, ~. +ST_RXYT +SR_RXYT)
					}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
						MinTerm <- names(which.min(AICSub))
						if(MinTerm == "ST_RXYT"){
							ModelFormula <- update(ModelFormula, ~. +SR_RXYT)
							mMin <- m1
						}
							if(MinTerm == "SR_RXYT"){
							ModelFormula <- update(ModelFormula, ~. +ST_RXYT)
							mMin <- m2
						}
						if(MinTerm == "ST_RXYT_SR_RXYT"){mMin <- m3}
					}
				}
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}


# if the temporal interaction has been removed, remove both of the terms separately and together...
			if(!is.null(ZeroInfl) & !"ST_T:SR_T" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_T)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Zi_ST_T <- try(AIC(m1))}
				ZiModelFormula <- update(ZiModelFormula, ~. +ST_T -SR_T)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Zi_SR_T <- try(AIC(m2))}
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_T -SR_T)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Zi_ST_T_SR_T <- try(AIC(m3))}
				AICSub <- AICs[, c("Zi_ST_T", "Zi_SR_T", "Zi_ST_T_SR_T")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ZiModelFormula <- update(ZiModelFormula, ~. +ST_T +SR_T)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "Zi_ST_T"){
						ZiModelFormula <- update(ZiModelFormula, ~. +SR_T)
						mMin <- m1
					}
					if(MinTerm == "Zi_SR_T"){
						ZiModelFormula <- update(ZiModelFormula, ~. +ST_T)
						mMin <- m2
					}
					if(MinTerm == "Zi_ST_T_SR_T"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
			

			if(!"ST_T:SR_T" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. -ST_T)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family) 
				if(!is.null(m1)){AICs$ST_T <- try(AIC(m1))}
				ModelFormula <- update(ModelFormula, ~. +ST_T -SR_T)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$SR_T <- try(AIC(m2))}
				ModelFormula <- update(ModelFormula, ~. -ST_T -SR_T)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$ST_T_SR_T <- try(AIC(m3))}
				AICSub <- AICs[, c("ST_T", "SR_T", "ST_T_SR_T")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ModelFormula <- update(ModelFormula, ~. +ST_T +SR_T)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "ST_T"){
						ModelFormula <- update(ModelFormula, ~. +SR_T)
						mMin <- m1
					}
						if(MinTerm == "SR_T"){
						ModelFormula <- update(ModelFormula, ~. +ST_T)
						mMin <- m2
					}
					if(MinTerm == "ST_T_SR_T"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
			

# if the spatial interaction has been removed, remove both of the terms separately and together...
			if(!is.null(ZeroInfl) & !"ST_XY:SR_XY" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_XY)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Zi_ST_XY <- try(AIC(m1))}
				ZiModelFormula <- update(ZiModelFormula, ~. +ST_XY -SR_XY)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Zi_SR_XY <- try(AIC(m2))}
				ZiModelFormula <- update(ZiModelFormula, ~. -ST_XY -SR_XY)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Zi_ST_XY_SR_XY <- try(AIC(m3))}
				AICSub <- AICs[, c("Zi_ST_XY", "Zi_SR_XY", "Zi_ST_XY_SR_XY")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ZiModelFormula <- update(ZiModelFormula, ~. +ST_XY +SR_XY)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "Zi_ST_XY"){
						ZiModelFormula <- update(ZiModelFormula, ~. +SR_XY)
						mMin <- m1
					}
						if(MinTerm == "Zi_SR_XY"){
						ZiModelFormula <- update(ZiModelFormula, ~. +ST_XY)
						mMin <- m2
					}
					if(MinTerm == "Zi_ST_XY_SR_XY"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}


			if(!"ST_XY:SR_XY" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. -ST_XY)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$ST_XY <- try(AIC(m1))}
				ModelFormula <- update(ModelFormula, ~. +ST_XY -SR_XY)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$SR_XY <- try(AIC(m2))}
				ModelFormula <- update(ModelFormula, ~. -ST_XY -SR_XY)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$ST_XY_SR_XY <- try(AIC(m3))}

				AICSub <- AICs[, c("ST_XY", "SR_XY", "ST_XY_SR_XY")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ModelFormula <- update(ModelFormula, ~. +ST_XY +SR_XY)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "ST_XY"){
						ModelFormula <- update(ModelFormula, ~. +SR_XY)
						mMin <- m1
					}
						if(MinTerm == "SR_XY"){
						ModelFormula <- update(ModelFormula, ~. +ST_XY)
						mMin <- m2
					}
					if(MinTerm == "ST_XY_SR_XY"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
			
				
# remove rad if present...
			if(!is.null(ZeroInfl) & "Rad" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. -Rad)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Rad <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Rad > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Rad)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Rad" %in% as.character(attr(terms(ModelFormula), "term.labels"))){			
				ModelFormula <- update(ModelFormula, ~. -Rad)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Rad <- try(AIC(m1))	
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Rad > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Rad)
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# remove slope...
			if(!is.null(ZeroInfl)){
				ZiModelFormula <- update(ZiModelFormula, ~. -Sl)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Sl <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Sl > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +Sl)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			ModelFormula <- update(ModelFormula, ~. -Sl)
			m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
			if(!is.null(m1)){
				AICs$Sl <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(AICs$Sl > MinAIC){
					ModelFormula <- update(ModelFormula, ~. +Sl)
				}else{
					mMin <- m1
				}
			}		
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)

# remove Agri if present and no Agri^2 is present
			if(!is.null(ZeroInfl) & !"I(Agri^2)" %in% as.character(attr(terms(ZiModelFormula), "variables")) & "Agri" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. -Agri)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Agri <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Agri > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. +Agri)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Agri" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"I(Agri^2)" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. -Agri)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Agri <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Agri > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Agri)
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# remove Urban.Infrastr if present
			if(!is.null(ZeroInfl) & "Urban.Infrastr" %in% as.character(attr(terms(ZiModelFormula), "term.labels"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Urban.Infrastr)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Urban.Infrastr <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Urban.Infrastr > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Urban.Infrastr)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Urban.Infrastr" %in% as.character(attr(terms(ModelFormula), "term.labels"))){
				ModelFormula <- update(ModelFormula, ~. - Urban.Infrastr)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Urban.Infrastr <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Urban.Infrastr > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Urban.Infrastr)
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# remove Water if present and Water^2 is not
			if(!is.null(ZeroInfl) & "Water" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & !"I(Water^2)" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Water)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family) 
				if(!is.null(m1)){
					AICs$Zi_Water <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Water > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Water)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Water" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"I(Water^2)" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Water)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Water <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Water > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Water)
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# remove Wetland if present and Wetland^2 is not
			if(!is.null(ZeroInfl) & "Wetland" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & !"I(Wetland^2)" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Wetland)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_Wetland <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_Wetland > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Wetland)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Wetland" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"I(Wetland^2)" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Wetland)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Wetland <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Wetland > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Wetland)
					}else{
						mMin <- m1
					}
				}		
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# if only one forest term is present, remove it
			if(!is.null(ZeroInfl) & "Other.for" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & !"Decid.for" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Other.for)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_OtherFor <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_OtherFor > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Other.for)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if(!is.null(ZeroInfl) & "Decid.for" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & !"Other.for" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Decid.for)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_DecidFor <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_DecidFor > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Decid.for)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# if forest terms are present, remove both of the terms separately and together...
			if(!is.null(ZeroInfl) & "Other.for" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & "Decid.for" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Other.for)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Zi_Other.for <- try(AIC(m1))}
				ZiModelFormula <- update(ZiModelFormula, ~. + Other.for - Decid.for)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Zi_Decid.for <- try(AIC(m2))}
				ZiModelFormula <- update(ZiModelFormula, ~. - Other.for - Decid.for)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Zi_Other_Decid.for <- try(AIC(m3))}

				AICSub <- AICs[, c("Zi_Other.for", "Zi_Decid.for", "Zi_Other_Decid.for")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ZiModelFormula <- update(ZiModelFormula, ~. + Other.for + Decid.for)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "Zi_Other.for"){
						ZiModelFormula <- update(ZiModelFormula, ~. + Decid.for)
						mMin <- m1
					}
						if(MinTerm == "Zi_Decid.for"){
						ZiModelFormula <- update(ZiModelFormula, ~. + Other.for)
						mMin <- m2
					}
					if(MinTerm == "Zi_Other_Decid.for"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
			

# if only one forest term is present, remove it
			if("Other.for" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"Decid.for" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Other.for)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$OtherFor <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$OtherFor > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Other.for)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Decid.for" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"Other.for" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Decid.for)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$DecidFor <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$DecidFor > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Decid.for)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# if both forest terms are present, remove each separately and both together
			if("Other.for" %in% as.character(attr(terms(ModelFormula), "term.labels")) & "Decid.for" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Other.for)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Other.for <- try(AIC(m1))}
				ModelFormula <- update(ModelFormula, ~. + Other.for - Decid.for)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Decid.for <- try(AIC(m2))}
				ModelFormula <- update(ModelFormula, ~. - Other.for - Decid.for)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Other_Decid.for <- try(AIC(m3))}

				AICSub <- AICs[, c("Other.for", "Decid.for", "Other_Decid.for")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ModelFormula <- update(ModelFormula, ~. + Other.for + Decid.for)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "Other.for"){
						ModelFormula <- update(ModelFormula, ~. + Decid.for)
						mMin <- m1
					}
						if(MinTerm == "Decid.for"){
						ModelFormula <- update(ModelFormula, ~. + Other.for)
						mMin <- m2
					}
					if(MinTerm == "Other_Decid.for"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
			

# if one fjell term is present, remove it
			if(!is.null(ZeroInfl) & "Fjall.open" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & !"Fjall.veg" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Fjall.open)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_FjallOpen <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_FjallOpen > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Fjall.open)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if(!is.null(ZeroInfl) & "Fjall.veg" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & !"Fjall.open" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Fjall.veg)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){
					AICs$Zi_FjallVeg <- try(AIC(m1))
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
					if(AICs$Zi_FjallVeg > MinAIC){
						ZiModelFormula <- update(ZiModelFormula, ~. + Fjall.veg)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# if both fjell terms are present, remove both of the terms separately and together...
			if(!is.null(ZeroInfl) & "Fjall.open" %in% as.character(attr(terms(ZiModelFormula), "term.labels")) & "Fjall.veg" %in% as.character(attr(terms(ZiModelFormula), "variables"))){
				ZiModelFormula <- update(ZiModelFormula, ~. - Fjall.open)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Zi_Fjall.open <- try(AIC(m1))}
				ZiModelFormula <- update(ZiModelFormula, ~. + Fjall.open - Fjall.veg)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Zi_Fjall.veg <- try(AIC(m2))}
				ZiModelFormula <- update(ZiModelFormula, ~. - Fjall.open - Fjall.veg)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Zi_Fjall.open_veg <- try(AIC(m3))}

				AICSub <- AICs[, c("Zi_Fjall.open", "Zi_Fjall.veg", "Zi_Fjall.open_veg")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ZiModelFormula <- update(ZiModelFormula, ~. + Fjall.open + Fjall.veg)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "Zi_Fjall.open"){
						ZiModelFormula <- update(ZiModelFormula, ~. + Fjall.veg)
						mMin <- m1
					}
						if(MinTerm == "Zi_Fjall.veg"){
						ZiModelFormula <- update(ZiModelFormula, ~. + Fjall.open)
						mMin <- m2
					}
					if(MinTerm == "Zi_Fjall.open_veg"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}
			

# if only one forest term is present, remove it
			if("Fjall.open" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"Fjall.veg" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Fjall.open)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$FjallOpen <- try(AIC(m1))}
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(!is.na(AICs$FjallOpen)){
					if(AICs$FjallOpen > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Fjall.open)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			if("Fjall.veg" %in% as.character(attr(terms(ModelFormula), "term.labels")) & !"Fjall.open" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Fjall.veg)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$FjallVeg <- try(AIC(m1))}
# ...if the reduced model has a higher AIC, add the removed term back in, otherwise the reduced  model becomes the best model.
				if(!is.na(AICs$FjallVeg)){
					if(AICs$FjallVeg > MinAIC){
						ModelFormula <- update(ModelFormula, ~. + Fjall.veg)
					}else{
						mMin <- m1
					}
				}			
				MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

# if both fjell terms are present, remove each separately and both together
			if("Fjall.open" %in% as.character(attr(terms(ModelFormula), "term.labels")) & "Fjall.veg" %in% as.character(attr(terms(ModelFormula), "variables"))){
				ModelFormula <- update(ModelFormula, ~. - Fjall.open)
				m1 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m1)){AICs$Fjall.open <- try(AIC(m1))}
				ModelFormula <- update(ModelFormula, ~. + Fjall.open - Fjall.veg)
				m2 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m2)){AICs$Fjall.veg <- try(AIC(m2))}
				ModelFormula <- update(ModelFormula, ~. - Fjall.open - Fjall.veg)
				m3 <- FunRunModel(MF = ModelFormula, ZMF = ZiModelFormula, CMD = CurModelData, Fam = Family)
				if(!is.null(m3)){AICs$Fjall.open_veg <- try(AIC(m3))}
				
				AICSub <- AICs[, c("Fjall.open", "Fjall.veg", "Fjall.open_veg")]
# if the model with the lowest AIC is the old one or all of the models m1-m3 had errors (MinSub will be infinite), add both terms back in and retain the previous mMin
				MinSub <- try(min(unlist(AICSub), na.rm = T))
				if(is.infinite(MinSub) | MinSub > MinAIC){
  					ModelFormula <- update(ModelFormula, ~. + Fjall.open + Fjall.veg)
				}else{
# if one of m1-m3 has the lowest AIC, update mMin and model formula
					MinTerm <- names(which.min(AICSub))
					if(MinTerm == "Fjall.open"){
						ModelFormula <- update(ModelFormula, ~. + Fjall.veg)
						mMin <- m1
					}
						if(MinTerm == "Fjall.veg"){
						ModelFormula <- update(ModelFormula, ~. + Fjall.open)
						mMin <- m2
					}
					if(MinTerm == "Fjall.open_veg"){mMin <- m3}
				}
			MinAIC <- min(AICs[, -(1:4)], na.rm = T)
			}

			try(m1 <- mMin)
			write.csv(AICs, paste0("AICs_", PdfName, "_", CurSpecInfo$englishname, ".csv"), row.names = F)
# save model
			ModelOutName <- gsub(" ", "_", CurSpecInfo$englishname)
			ModelListMAM <- append(ModelListMAM, list(NextModel = m1))
			names(ModelListMAM)[length(names(ModelListMAM))] <- paste("M", ModelOutName, sep = "_")
			save(list = "ModelListMAM", file = paste0(PdfName, "_MAM_", CurSpecInfo$englishname, ".RData"))

		}
	}
}

##########################################################################################################################################################################################################
###################################################### Extract coefficient estimates and standard errors from a batch of models ##########################################################################
##########################################################################################################################################################################################################

FunExtractBetas <- function(Path, MyModelList, Type){
# set up empty data frames..
	Covars <- c("(Intercept)", "ST_XY", "SR_XY", "ST_XY:SR_XY", "ST_T", "SR_T", "ST_T:SR_T", "ST_RXYT", "SR_RXYT", "ST_RXYT:SR_RXYT", "I(yrStand^2)", "yrStand", "Sl", "Rad", "I(Water^2)", 
"Water", "I(Wetland^2)", "Wetland", "Fjall.veg", "Fjall.open", "Decid.for", "Other.for", "I(Agri^2)", "Agri", "Urban.Infrastr", "SurveyPoint", "UnitPair")
# add zi-covars
	Covars <- c(Covars, paste0("zi.", Covars), "zi.EffortStand")
# ... one data frame for the coefficients...
	OutBetas <- as.data.frame(matrix(nrow = length(MyModelList), ncol = length(Covars) + 1))
	colnames(OutBetas) <- c("Species", Covars)
# ... and one data frame for the standard errors
	OutSEs <- OutBetas
	for(i in 1 : length(MyModelList)){
		CurModel <- load(paste0(CurPath, "/", MyModelList[i]))
		if(Type == "MAM"){
			CurModel <- ModelListMAM[[1]]
		}else{
			CurModel <- ModelListFull[[1]]
		}
		Ind <- regexpr(Type, MyModelList[i])
# extract species name from file name
		CurSpec <- substr(MyModelList[i], Ind, nchar(MyModelList[i]))
		CurSpec <- gsub(".RData", "", CurSpec)
		CurSpec <- gsub(paste0(Type, "_"), "", CurSpec)
		OutBetas$Species[i] <- CurSpec
		OutSEs$Species[i] <- CurSpec
# extract coefficient estimates and standard errors for the betas in the count and zero-inflation parts
		SEs <- summary(CurModel$sdr, "fixed")
		SEs <- SEs[dimnames(SEs)[[1]] == "beta" | dimnames(SEs)[[1]] == "betazi", ]
		SEs <- as.data.frame(SEs)
# extract the covariate names
		CoefList <- fixef(CurModel)	
		Coefs <- CoefList$cond
		if(length(CoefList$zi) > 0){
			names(CoefList$zi) <- paste0("zi.", names(CoefList$zi))
			Coefs <- c(Coefs, CoefList$zi)
		}
		rownames(SEs) <- names(Coefs)
		colnames(SEs) <- c("Estimate", "StdError")
# check that the covariate names are matched correctly
		if(all((SEs$Estimate - Coefs) == 0) == F){print("Covariate names not matched correctly")}
# add the coefficient estimates and standard errors to the data frames
		OutBetas[i, match(rownames(SEs), colnames(OutBetas))] <- SEs$Estimate
		OutSEs[i, match(rownames(SEs), colnames(OutSEs))] <- SEs$StdError
	}
	Out <- list(Betas = OutBetas, StdErrors = OutSEs)
	return(Out)
}	

##########################################################################################################################################################################################################
###################################################### Standardize based on the means and sds from the data used to fit models ##########################################################################
##########################################################################################################################################################################################################

FunStandWithModelData <- function(InData, StandMeanSD){
	InData$ST_XY <- (InData$ST_XY - StandMeanSD$M_ST_XY) / StandMeanSD$S_ST_XY
	InData$SR_XY <- (InData$SR_XY - StandMeanSD$M_SR_XY) / StandMeanSD$S_SR_XY
	InData$ST_T <- (InData$ST_T - StandMeanSD$M_ST_T) / StandMeanSD$S_ST_T
	InData$SR_T <- (InData$SR_T - StandMeanSD$M_SR_T) / StandMeanSD$S_SR_T
	InData$ST_RXYT <- (InData$ST_RXYT - StandMeanSD$M_ST_RXYT) / StandMeanSD$S_ST_RXYT
	InData$SR_RXYT <- (InData$SR_RXYT - StandMeanSD$M_SR_RXYT) / StandMeanSD$S_SR_RXYT

	InData$DEM <- (InData$DEM - StandMeanSD$M_DEM) / StandMeanSD$S_DEM
	InData$Sl <- (InData$Sl - StandMeanSD$M_Sl) / StandMeanSD$S_Sl
	InData$Rad <- (InData$Rad - StandMeanSD$M_Rad) / StandMeanSD$S_Rad
	InData$Urban.Infrastr <- (InData$Urban.Infrastr - StandMeanSD$M_Urban.Infrastr) / StandMeanSD$S_Urban.Infrastr
	InData$Agri <- (InData$Agri - StandMeanSD$M_Agri) / StandMeanSD$S_Agri
	InData$Decid.for <- (InData$Decid.for - StandMeanSD$M_Decid.for) / StandMeanSD$S_Decid.for
	InData$Other.for <- (InData$Other.for - StandMeanSD$M_Other.for) / StandMeanSD$S_Other.for
	InData$Fjall.veg <- (InData$Fjall.veg - StandMeanSD$M_Fjall.veg) / StandMeanSD$S_Fjall.veg
	InData$Fjall.open <- (InData$Fjall.open - StandMeanSD$M_Fjall.open) / StandMeanSD$S_Fjall.open
	InData$Wetland <- (InData$Wetland - StandMeanSD$M_Wetland) / StandMeanSD$S_Wetland
	InData$Water <- (InData$Water - StandMeanSD$M_Water) / StandMeanSD$S_Water
	return(InData)
}

##########################################################################################################################################################################################################
########################################################################## Create sketch data for model predictions #############################################################################
##########################################################################################################################################################################################################

FunCreateSketchCovars <- function(InCovarData){
	SketchCovars <- InCovarData
	SketchCovars <- SketchCovars[SketchCovars$Survey == "Line", ]					# 1749 survey routes with lines; including also points, there are 1756 survey routes
# sort by year and select unique routes
	SketchCovars <- SketchCovars[order(SketchCovars$yr),  ]
	SketchCovars <- SketchCovars[!duplicated(SketchCovars$Route), ]
	SketchCovars$Unit <- "Ind"
	SketchCovars$Unit <- factor(SketchCovars$Unit, levels = c("Ind", "Pair"))
	return(SketchCovars)
}

##########################################################################################################################################################################################################
########################################################################## Calculate Spearman rho between climate components #############################################################################
##########################################################################################################################################################################################################

FunCalcSpearmanRhos <- function(ModelPath, ModelList, CovarData, RangeData, StandMeanSD){
	OutCors <- as.data.frame(matrix(ncol = 10, nrow = length(ModelList)))
	colnames(OutCors) <- c("Species", "Sp_T", "Sp_Res", "T_Res", "Temp_Sp_T", "Temp_Sp_Res", "Temp_T_Res", "Prec_Sp_T", "Prec_Sp_Res", "Prec_T_Res")
	for(i in 1 : length(ModelList)){
# get current mdoel
		CurModel <- load(paste(ModelPath, ModelList[i], sep = "/"))
		CurModel <- ModelListMAM[[1]]
		CurSpecies <- ModelList[i]
		CurSpecies <- gsub(".RData", "", CurSpecies)
		CurSpecies <- gsub("CurMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("CurMayJune_MAM_", "", CurSpecies)
		CurSpecies <- gsub("PrevMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("_", " ", CurSpecies)
		OutCors$Species[i] <- CurSpecies
# create the sketch data...
# ...for spatial...
		Min_ST_XY <- min(RangeData$ST_XY)
		Max_ST_XY <- max(RangeData$ST_XY)
		Min_SR_XY <- min(RangeData$SR_XY)
		Max_SR_XY <- max(RangeData$SR_XY)
		AddDataSp <- expand.grid(ST_XY = seq(Min_ST_XY, Max_ST_XY, length.out = 10), SR_XY = seq(Min_SR_XY, Max_SR_XY, length.out = 10))
# ...for temporal...
		Min_ST_T <- min(RangeData$ST_T)
		Max_ST_T <- max(RangeData$ST_T)
		Min_SR_T <- min(RangeData$SR_T)
		Max_SR_T <- max(RangeData$SR_T)
		AddDataT <- expand.grid(ST_T = seq(Min_ST_T, Max_ST_T, length.out = 10), SR_T = seq(Min_SR_T, Max_SR_T, length.out = 10))
# ...for residual
		Min_ST_RXYT <- min(RangeData$ST_RXYT)
		Max_ST_RXYT <- max(RangeData$ST_RXYT)
		Min_SR_RXYT <- min(RangeData$SR_RXYT)
		Max_SR_RXYT <- max(RangeData$SR_RXYT)
		AddDataRes <- expand.grid(ST_RXYT = seq(Min_ST_RXYT, Max_ST_RXYT, length.out = 10), SR_RXYT = seq(Min_SR_RXYT, Max_SR_RXYT, length.out = 10))
# repeat the sketch data for as many times as there are AddData points
		CurCovarDataSp <- CovarData[rep(seq_len(nrow(CovarData)), nrow(AddDataSp)), ]
		CurCovarDataT <- CovarData[rep(seq_len(nrow(CovarData)), nrow(AddDataT)), ]
		CurCovarDataRes <- CovarData[rep(seq_len(nrow(CovarData)), nrow(AddDataRes)), ]
# repeat the AddData as many times as there are sketch data points
		AddDataSp <- AddDataSp[rep(seq_len(nrow(AddDataSp)), each = nrow(CovarData)), ]
		AddDataT <- AddDataT[rep(seq_len(nrow(AddDataT)), each = nrow(CovarData)), ]
		AddDataRes <- AddDataRes[rep(seq_len(nrow(AddDataRes)), each = nrow(CovarData)), ]
# add together
		CurCovarDataSp$ST_XY <- NULL
		CurCovarDataSp$SR_XY <- NULL
		CurCovarDataSp <- cbind(CurCovarDataSp, AddDataSp)

		CurCovarDataT$ST_T <- NULL
		CurCovarDataT$SR_T <- NULL
		CurCovarDataT <- cbind(CurCovarDataT, AddDataT)

		CurCovarDataRes$ST_RXYT <- NULL
		CurCovarDataRes$SR_RXYT <- NULL
		CurCovarDataRes <- cbind(CurCovarDataRes, AddDataRes)

# standardize covariate data
		CurPredDataSp <- FunStandWithModelData(CurCovarDataSp, StandMeanSD)
		CurPredDataT <- FunStandWithModelData(CurCovarDataT, StandMeanSD)
		CurPredDataRes <- FunStandWithModelData(CurCovarDataRes, StandMeanSD)
# set year and Effort to 0 / 8
		CurPredDataSp$yrStand <- 0
		CurPredDataSp$EffortStand <- 0
		CurPredDataSp$Effort <- 8
		CurPredDataT$yrStand <- 0
		CurPredDataT$EffortStand <- 0
		CurPredDataT$Effort <- 8
		CurPredDataRes$yrStand <- 0
		CurPredDataRes$EffortStand <- 0
		CurPredDataRes$Effort <- 8
# set T and Res to 0 for spatial component
		CurPredDataSp$ST_T <- 0
		CurPredDataSp$SR_T <- 0
		CurPredDataSp$ST_RXYT <- 0
		CurPredDataSp$SR_RXYT <- 0
# set Res to 0 for temporal component
		CurPredDataT$ST_RXYT <- 0
		CurPredDataT$SR_RXYT <- 0
# set T to 0 for residual component
		CurPredDataRes$ST_T <- 0
		CurPredDataRes$SR_T <- 0
# predict on standardized data and attach to unstandardized data
		CurCovarDataSp$Preds <-  predict(CurModel, CurPredDataSp, type = "response")
		CurCovarDataT$Preds <-  predict(CurModel, CurPredDataT, type = "response")
		CurCovarDataRes$Preds <-  predict(CurModel, CurPredDataRes, type = "response")

# aggregate per Temperature and precipitation
		SpPredsPerTempPrec <- aggregate(CurCovarDataSp$Preds, by = list(CurCovarDataSp$ST_XY, CurCovarDataSp$SR_XY), sum)		# elements of by are coerced to factors
		TPredsPerTempPrec <- aggregate(CurCovarDataT$Preds, by = list(CurCovarDataT$ST_T, CurCovarDataT$SR_T), sum)	
		ResPredsPerTempPrec <- aggregate(CurCovarDataRes$Preds, by = list(CurCovarDataRes$ST_RXYT, CurCovarDataRes$SR_RXYT), sum)

# calculate Spearman rhos over the temperature and precipitation gradient
		OutCors$Sp_T[i] <- cor(SpPredsPerTempPrec$x, TPredsPerTempPrec$x, method = "spearman")
		OutCors$Sp_Res[i] <- cor(SpPredsPerTempPrec$x, ResPredsPerTempPrec$x, method = "spearman")	
		OutCors$T_Res[i] <- cor(TPredsPerTempPrec$x, ResPredsPerTempPrec$x, method = "spearman")

# calculate Spearman rhos over the temperature gradient...
# make sure they are ordered first by precipitation, then by temperature
		SpPredsPerTempPrec <- SpPredsPerTempPrec[order(SpPredsPerTempPrec$Group.2, SpPredsPerTempPrec$Group.1), ]	
		TPredsPerTempPrec <- TPredsPerTempPrec[order(TPredsPerTempPrec$Group.2, TPredsPerTempPrec$Group.1), ]	
		ResPredsPerTempPrec <- ResPredsPerTempPrec[order(ResPredsPerTempPrec$Group.2, ResPredsPerTempPrec$Group.1), ]
	
# the first ten rows have a temp gradient and prec is fixed to the lowest amount - the last ten rows have a temp gradient and prec is fixed to the highest amount
# I calculate Spearman rho for each temp gradient and fixed prec
		CorsSp_T <- numeric(10)
		CorsSp_Res <- numeric(10)
		CorsT_Res <- numeric(10)
		for(j in 1 : 10){
			CorsSp_T[j] <- cor(SpPredsPerTempPrec$x[(j*10-9) : (j*10)], TPredsPerTempPrec$x[(j*10-9) : (j*10)], method = "spearman")
			CorsSp_Res[j] <- cor(SpPredsPerTempPrec$x[(j*10-9) : (j*10)], ResPredsPerTempPrec$x[(j*10-9) : (j*10)], method = "spearman")
			CorsT_Res[j] <- cor(TPredsPerTempPrec$x[(j*10-9) : (j*10)], ResPredsPerTempPrec$x[(j*10-9) : (j*10)], method = "spearman")
		}
# I set NA values to 0 (NA are the result of no variation in either of the vectors)
		CorsSp_T[is.na(CorsSp_T)] <- 0
		CorsSp_Res[is.na(CorsSp_Res)] <- 0
		CorsT_Res[is.na(CorsT_Res)] <- 0
# I calculate and save the mean of the 10 Spearman rhos
		OutCors$Temp_Sp_T[i] <- mean(CorsSp_T)
		OutCors$Temp_Sp_Res[i] <- mean(CorsSp_Res)
		OutCors$Temp_T_Res[i] <- mean(CorsT_Res)
	
# calculate Spearman rhos over the precipitation gradient...
# make sure they are ordered first by temperature, then by precipitation
		SpPredsPerTempPrec <- SpPredsPerTempPrec[order(SpPredsPerTempPrec$Group.1, SpPredsPerTempPrec$Group.2), ]	
		TPredsPerTempPrec <- TPredsPerTempPrec[order(TPredsPerTempPrec$Group.1, TPredsPerTempPrec$Group.2), ]	
		ResPredsPerTempPrec <- ResPredsPerTempPrec[order(ResPredsPerTempPrec$Group.1, ResPredsPerTempPrec$Group.2), ]
	
# the first ten rows have a prec gradient and temp is fixed to the lowest amount - the last ten rows have a prec gradient and temp is fixed to the highest amount
# I calculate Spearman rho for each prec gradient and fixed temp
		CorsSp_T <- numeric(10)
		CorsSp_Res <- numeric(10)
		CorsT_Res <- numeric(10)
		for(j in 1 : 10){
			CorsSp_T[j] <- cor(SpPredsPerTempPrec$x[(j*10-9) : (j*10)], TPredsPerTempPrec$x[(j*10-9) : (j*10)], method = "spearman")
			CorsSp_Res[j] <- cor(SpPredsPerTempPrec$x[(j*10-9) : (j*10)], ResPredsPerTempPrec$x[(j*10-9) : (j*10)], method = "spearman")
			CorsT_Res[j] <- cor(TPredsPerTempPrec$x[(j*10-9) : (j*10)], ResPredsPerTempPrec$x[(j*10-9) : (j*10)], method = "spearman")
		}
# I set NA values to 0 (NA are the result of no variation in either of the vectors)
		CorsSp_T[is.na(CorsSp_T)] <- 0
		CorsSp_Res[is.na(CorsSp_Res)] <- 0
		CorsT_Res[is.na(CorsT_Res)] <- 0
# I calculate and save the mean of the 10 Spearman rhos
		OutCors$Prec_Sp_T[i] <- mean(CorsSp_T)
		OutCors$Prec_Sp_Res[i] <- mean(CorsSp_Res)
		OutCors$Prec_T_Res[i] <- mean(CorsT_Res)	
	}
# replace NA values with 0 for the Spearman rho over both the temp and prec gradient
	OutCors$Sp_T[is.na(OutCors$Sp_T)] <- 0
	OutCors$Sp_Res[is.na(OutCors$Sp_Res)] <- 0
	OutCors$T_Res[is.na(OutCors$T_Res)] <- 0
	return(OutCors)
}

###########################################################################################################################################################################################
############################################################################### create spatial raster #####################################################################################
###########################################################################################################################################################################################

FunSketchModelsSpatial <- function(ModelPath, ModelList, CovarData, RangeData, StandMeanSD){
	for(i in 1 : length(ModelList)){
# get current mdoel
	# get current mdoel
		CurModel <- load(paste(ModelPath, ModelList[i], sep = "/"))
		CurModel <- ModelListMAM[[1]]
		CurSpecies <- ModelList[i]
		CurSpecies <- gsub(".RData", "", CurSpecies)
		CurSpecies <- gsub("CurMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("CurMayJune_MAM_", "", CurSpecies)
		CurSpecies <- gsub("PrevMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("_", " ", CurSpecies)
		Min_ST_XY <- min(RangeData$ST_XY)
		Max_ST_XY <- max(RangeData$ST_XY)
		Min_SR_XY <- min(RangeData$SR_XY)
		Max_SR_XY <- max(RangeData$SR_XY)
		AddData <- expand.grid(ST_XY = seq(Min_ST_XY, Max_ST_XY, length.out = 10), SR_XY = seq(Min_SR_XY, Max_SR_XY, length.out = 10))
# repeat the sketch data for as many times as there are AddData points
		CurCovarData <- CovarData[rep(seq_len(nrow(CovarData)), nrow(AddData)), ]
# repeat the AddData as many times as there are sketch data points
		AddData <- AddData[rep(seq_len(nrow(AddData)), each = nrow(CovarData)), ]
# add together
		CurCovarData$ST_XY <- NULL
		CurCovarData$SR_XY <- NULL
		CurCovarData <- cbind(CurCovarData, AddData)
# standardize covariate data
		CurPredData <- FunStandWithModelData(CurCovarData, StandMeanSD)
# set year and Effort to 0 / 8
		CurPredData$yrStand <- 0
		CurPredData$EffortStand <- 0
		CurPredData$Effort <- 8
		CurPredData$ST_T <- 0
		CurPredData$SR_T <- 0
		CurPredData$ST_RXYT <- 0
		CurPredData$SR_RXYT <- 0
# predict on standardized data and attach to unstandardized data
		CurCovarData$Preds <-  predict(CurModel, CurPredData, type = "response")
		CurCovarData$ST_XY <- factor(CurCovarData$ST_XY)
		CurCovarData$SR_XY <- factor(CurCovarData$SR_XY)
# For each combination of the AddData, I sum the predictions over all 1749 routes.
		CurPlotData <- aggregate(CurCovarData$Preds, by = list(CurCovarData$ST_XY, CurCovarData$SR_XY), sum)
		colnames(CurPlotData) <- c("ST_XY", "SR_XY", "Abundance")
		CurPlotData$ST_XY <- as.numeric(as.character(CurPlotData$ST_XY))
		CurPlotData$SR_XY <- as.numeric(as.character(CurPlotData$SR_XY))
		X_Axis <- round(unique(CurPlotData$ST_XY), 1)
		Y_Axis <- round(unique(CurPlotData$SR_XY), 1)

		P1 <- ggplot(CurPlotData, aes(x = ST_XY, y = SR_XY)) + geom_raster(aes(fill = Abundance)) + theme_void() + 
		scale_fill_gradientn(colours = c("grey90", "grey20")) + theme(legend.position = "none")
		ggsave(paste0("Spatial_", CurSpecies, ".pdf"), P1)
	}
}

###########################################################################################################################################################################################
############################################################################### create temporal raster #####################################################################################
###########################################################################################################################################################################################

FunSketchModelsTemporal <- function(ModelPath, ModelList, CovarData, RangeData, StandMeanSD){
	for(i in 1 : length(ModelList)){
# get current mdoel
		CurModel <- load(paste(ModelPath, ModelList[i], sep = "/"))
		CurModel <- ModelListMAM[[1]]
		CurSpecies <- ModelList[i]
		CurSpecies <- gsub(".RData", "", CurSpecies)
		CurSpecies <- gsub("CurMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("CurMayJune_MAM_", "", CurSpecies)
		CurSpecies <- gsub("PrevMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("_", " ", CurSpecies)
		Min_ST_T <- min(RangeData$ST_T)
		Max_ST_T <- max(RangeData$ST_T)
		Min_SR_T <- min(RangeData$SR_T)
		Max_SR_T <- max(RangeData$SR_T)
		AddData <- expand.grid(ST_T = seq(Min_ST_T, Max_ST_T, length.out = 10), SR_T = seq(Min_SR_T, Max_SR_T, length.out = 10))
# repeat the sketch data for as many times as there are AddData points
		CurCovarData <- CovarData[rep(seq_len(nrow(CovarData)), nrow(AddData)), ]
# repeat the AddData as many times as there are sketch data points
		AddData <- AddData[rep(seq_len(nrow(AddData)), each = nrow(CovarData)), ]
# add together
		CurCovarData$ST_T <- NULL
		CurCovarData$SR_T <- NULL
		CurCovarData <- cbind(CurCovarData, AddData)
# standardize covariate data
		CurPredData <- FunStandWithModelData(CurCovarData, StandMeanSD)
# set year and Effort to 0 / 8
		CurPredData$yrStand <- 0
		CurPredData$EffortStand <- 0
		CurPredData$Effort <- 8
		CurPredData$ST_RXYT <- 0
		CurPredData$SR_RXYT <- 0
# predict on standardized data and attach to unstandardized data
		CurCovarData$Preds <-  predict(CurModel, CurPredData, type = "response")
		CurCovarData$ST_T <- factor(CurCovarData$ST_T)
		CurCovarData$SR_T <- factor(CurCovarData$SR_T)
# For each combination of the AddData, I sum the predictions over all 1749 routes.
		CurPlotData <- aggregate(CurCovarData$Preds, by = list(CurCovarData$ST_T, CurCovarData$SR_T), sum)
		colnames(CurPlotData) <- c("ST_T", "SR_T", "Abundance")
		CurPlotData$ST_T <- as.numeric(as.character(CurPlotData$ST_T))
		CurPlotData$SR_T <- as.numeric(as.character(CurPlotData$SR_T))
		X_Axis <- round(unique(CurPlotData$ST_T), 1)
		Y_Axis <- round(unique(CurPlotData$SR_T), 1)

		P1 <- ggplot(CurPlotData, aes(x = ST_T, y = SR_T)) + geom_raster(aes(fill = Abundance)) + theme_void() + scale_fill_gradientn(colours = c("grey90", "grey20")) + theme(legend.position = "none")
		ggsave(paste0("Temporal_", CurSpecies, ".pdf"), P1)
	}
}

###########################################################################################################################################################################################
############################################################################### create residual raster #####################################################################################
###########################################################################################################################################################################################

FunSketchModelsResidual <- function(ModelPath, ModelList, CovarData, RangeData, StandMeanSD){
	for(i in 1 : length(ModelList)){
# get current mdoel
		CurModel <- load(paste(ModelPath, ModelList[i], sep = "/"))
		CurModel <- ModelListMAM[[1]]
		CurSpecies <- ModelList[i]
		CurSpecies <- gsub(".RData", "", CurSpecies)
		CurSpecies <- gsub("CurMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("CurMayJune_MAM_", "", CurSpecies)
		CurSpecies <- gsub("PrevMayJuly_MAM_", "", CurSpecies)
		CurSpecies <- gsub("_", " ", CurSpecies)
		Min_ST_RXYT <- min(RangeData$ST_RXYT)
		Max_ST_RXYT <- max(RangeData$ST_RXYT)
		Min_SR_RXYT <- min(RangeData$SR_RXYT)
		Max_SR_RXYT <- max(RangeData$SR_RXYT)
		AddData <- expand.grid(ST_RXYT = seq(Min_ST_RXYT, Max_ST_RXYT, length.out = 10), SR_RXYT = seq(Min_SR_RXYT, Max_SR_RXYT, length.out = 10))
# repeat the sketch data for as many times as there are AddData points
		CurCovarData <- CovarData[rep(seq_len(nrow(CovarData)), nrow(AddData)), ]
# repeat the AddData as many times as there are sketch data points
		AddData <- AddData[rep(seq_len(nrow(AddData)), each = nrow(CovarData)), ]
# add together
		CurCovarData$ST_RXYT <- NULL
		CurCovarData$SR_RXYT <- NULL
		CurCovarData <- cbind(CurCovarData, AddData)
# standardize covariate data
		CurPredData <- FunStandWithModelData(CurCovarData, StandMeanSD)
# set year and Effort to 0 / 8
		CurPredData$yrStand <- 0
		CurPredData$EffortStand <- 0
		CurPredData$Effort <- 8
		CurPredData$ST_T <- 0
		CurPredData$SR_T <- 0
# predict on standardized data and attach to unstandardized data
		CurCovarData$Preds <-  predict(CurModel, CurPredData, type = "response")
		CurCovarData$ST_RXYT <- factor(CurCovarData$ST_RXYT)
		CurCovarData$SR_RXYT <- factor(CurCovarData$SR_RXYT)
# For each combination of the AddData, I sum the predictions over all 1749 routes.
		CurPlotData <- aggregate(CurCovarData$Preds, by = list(CurCovarData$ST_RXYT, CurCovarData$SR_RXYT), sum)
		colnames(CurPlotData) <- c("ST_RXYT", "SR_RXYT", "Abundance")
		CurPlotData$ST_RXYT <- as.numeric(as.character(CurPlotData$ST_RXYT))
		CurPlotData$SR_RXYT <- as.numeric(as.character(CurPlotData$SR_RXYT))
		X_Axis <- round(unique(CurPlotData$ST_RXYT), 1)
		Y_Axis <- round(unique(CurPlotData$SR_RXYT), 1)

		P1 <- ggplot(CurPlotData, aes(x = ST_RXYT, y = SR_RXYT)) + geom_raster(aes(fill = Abundance)) + theme_void() + scale_fill_gradientn(colours = c("grey90", "grey20")) + theme(legend.position = "none")
		ggsave(paste0("Residual_", CurSpecies, ".pdf"), P1)
	}
}

###########################################################################################################################################################################################
################################################################################# Cross-validation ########################################################################################
###########################################################################################################################################################################################

FunCV <- function(ModelPath, ModelList, SpeciesData, SpeciesInfo, CovarData, Season, Year, PdfName){
	CVSummary <- as.data.frame(matrix(nrow = length(ModelList), ncol = 6))
	colnames(CVSummary) <- c("Species", "Pears_Ind_All", "Pears_Ind_No", "Pears_Ind_Swe", "Pears_Ind_Fi", "CV_Pears")
	for(i in 1 : length(ModelList)){
		CurModel <- load(paste(ModelPath, ModelList[i], sep = "/"))
		CurModel <- ModelListMAM[[1]]
		CurSpecies <- ModelList[i]
		CurSpecies <- gsub(".RData", "", CurSpecies)
		CurSpecies <- gsub("CurMayJuly_MAM_", "", CurSpecies)
		CVSummary$Species[i] <- CurSpecies		
# get current species info
		CurSpecInfo <- SpeciesInfo[SpeciesInfo$englishname == CVSummary$Species[i], ]
# extract species counts
		CurSpecData <- SpeciesData[SpeciesData$euring == CurSpecInfo$Euring, ]
# identify the covariate data
		if(is.null(Season)){
			SeasonShort <- ifelse(CurSpecInfo$Fjall_skog_area == "Summer", "S", "Y")
		}else{
			if(!Season %in% c("Summer", "Year"))stop("Argument to Season must be either Summer or Year")
			SeasonShort <- substr(Season, 1, 1)
		}
		CurCovData <- paste0(SeasonShort, ifelse(Year == "Current", "C", "P"))
		CurCovData <- CovarData[[which(names(CovarData) == CurCovData)]]
# standardize
		CurCovData <- FunStand(CurCovData)
# combine species and covariate data
		CurModelData <- merge(CurSpecData[, c("Route", "yr", "Count", "Survey")], CurCovData, by = c("Route", "yr", "Survey"), all.y = T)
# zero-fill
		CurModelData[is.na(CurModelData$Count), ]$Count <- 0
# compare fitted with observed...
# calculate predicted values at individual level
		CurModelData$PredsInd <- predict(CurModel, CurModelData, type = "response")		
# ... globally...		
		CVSummary$Pears_Ind_All[i] <- cor(CurModelData$Count, CurModelData$PredsInd, method = "pearson")
# ... per country...
		Norway <- CurModelData[substr(CurModelData$Route, 1, 1) == "N", ]
		CVSummary$Pears_Ind_No[i] <- cor(Norway$Count, Norway$PredsInd, method = "pearson")
		Finland <- CurModelData[substr(CurModelData$Route, 1, 1) == "F", ]
		CVSummary$Pears_Ind_Fi[i] <- cor(Finland$Count, Finland$PredsInd, method = "pearson")
		Sweden <- CurModelData[substr(CurModelData$Route, 1, 1) != "F" & substr(CurModelData$Route, 1, 1) != "N", ]
		CVSummary$Pears_Ind_Swe[i] <- cor(Sweden$Count, Sweden$PredsInd, method = "pearson")
# Cross-validation for years 2010 - 2018 (because early years have no/few routes in Norway and Finland...
# ...extract model formulas
		CallText <- unlist(strsplit(as.character(CurModel$call), split = ","))
		ModelFormula <- formula(CallText[2])
		ZiModelFormula <- formula(CallText[5])
# ...extract Family
		Fam <- summary(CurModel)$family
		TestOut <- NULL
		for(yr in 1996 : 2018){
			Train <- CurModelData[CurModelData$yr != yr, ]
			Test <- CurModelData[CurModelData$yr == yr, ]
			m1 <- NULL
			try(m1 <- glmmTMB(ModelFormula, ziformula = ZiModelFormula, data = Train, family = Fam))
			if(class(m1) == "glmmTMB" & !is.na(try(AIC(m1)))){
		  		try(Test$Preds <- predict(m1, Test, type = "response", allow.new.levels = T))
		  		TestOut <- rbind(TestOut, Test)
		  	}
		}		
		CVSummary$CV_Pears[i] <- try(cor(TestOut$Count, TestOut$Preds, method = "pearson"))
# safety write out
		write.csv(CVSummary, "CVSummary_SafetyWriteOut.csv", row.names = F)
	}		# end for each mdoel
	return(CVSummary)
}
