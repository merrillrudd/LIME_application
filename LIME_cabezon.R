rm(list=ls())

## load LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

## load LBSPR package
library(LBSPR)

## load Rfishbase
library(rfishbase)

###################################
## Directories
###################################
main_dir <- "C:\\Git_Projects\\LIME_application"

R_dir <- file.path(main_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

###################################
## Identify species
###################################

stock_name <- "cabezon_NCS"
proper_name <- "Scorpaenichthys marmoratus"

###################################
## Species directories
###################################

sp_dir <- file.path(main_dir, stock_name)
dir.create(sp_dir, showWarnings=FALSE)

data_dir <- file.path(sp_dir, "data")
dir.create(data_dir, showWarnings=FALSE)

res_dir <- file.path(sp_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

figs_dir <- file.path(sp_dir, "figures")
dir.create(figs_dir, showWarnings=FALSE)

###################################
## Life history info
###################################

## check fishbase
lw <- length_weight(proper_name)
lwa_med <- median(lw$a, na.rm=TRUE)
lwb_med <- median(lw$b, na.rm=TRUE)

growth <- popgrowth(proper_name)
linf_med <- median(growth$TLinfinity, na.rm=TRUE)
if(is.na(linf_med)) linf_med <- median(growth$Loo, na.rm=TRUE)
vbk_med <- median(growth$K, na.rm=TRUE)
winf_med <- median(growth$Winfinity, na.rm=TRUE)
temp_med <- median(growth$Temperature)

## choose values (from stock assessment)
linf_toUse <- 68.8
vbk_toUse <- 0.2
t0_toUse <- -1.19
lwa_toUse <- 9.2e-6
lwb_toUse <- 3.187

## use natural mortality tool
M_toUse <- 0.25

## add life history information for species chosen for assessment
## siganus sutor
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, lwa=lwa_toUse, lwb=lwb_toUse, M=M_toUse, S50=34.6, selex_input="length", M50=34.6, maturity_input="length", binwidth=2, nseasons=1)

par(mfrow=c(2,2))
plot(plist$L_a, type="l", lwd=3)
plot(plist$W_a, type="l", lwd=3)
plot(plist$Mat_l, type="l", lwd=3)
plot(plist$S_l, type="l", lwd=3)

###################################
## Load data
###################################
## read raw data
lc <- as.matrix(read.csv(file.path(data_dir, "Cabezon_LIME_LC.csv"), row=1))
LC <- lc[which(rowSums(lc)>0),]
bins <- seq(2, ncol(LC)*2, by=2)
colnames(LC) <- bins

plot_LC(Inputs=list("LF"=LC))


## abundance index
It_raw <- read.csv(file.path(data_dir, "Cabezon_LIME_It.csv"))
I_t <- as.vector(It_raw[,2])
names(I_t) <- It_raw[,1]
I_t <- I_t[which(is.na(I_t)==FALSE)]

###################################
## Run assessments
###################################

####################################
## LC only
####################################

lconly_dir <- file.path(res_dir, "lc_only")
dir.create(lconly_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LC))
years_t <- min(years_o):max(years_o)

input_data <- list("years"=years_t, "LF"=LC, "I_t"=I_t)

		####################################
		## LC only
		####################################
		lconly_lbspr <- file.path(lconly_dir, "LBSPR")
		dir.create(lconly_lbspr, showWarnings=FALSE)

		run <- run_LBSPR(modpath=lconly_lbspr, lh=plist, species=stock_name, input_data=input_data, rewrite=TRUE, simulation=FALSE, write=TRUE)	
		lbspr_res <- readRDS(file.path(lconly_lbspr, "LBSPR_results.rds"))

		####################################
		## LC only, LOGISTIC SELECTIVITY
		####################################
		slogis_lconly <- file.path(lconly_dir, "selex_logistic")
		dir.create(slogis_lconly, showWarnings=FALSE)

			##############################################
			## LC only, LOGISTIC SELECTIVITY - DEFAULT
			##############################################
			slogis_lconly_default <- file.path(slogis_lconly, "default")
			dir.create(slogis_lconly_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=slogis_lconly_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2))
				
				Report <- readRDS(file.path(slogis_lconly_default, "Report.rds"))
				Sdreport <- readRDS(file.path(slogis_lconly_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(slogis_lconly_default, "Inputs.rds"))

				png(file.path(slogis_lconly_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(slogis_lconly_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()


		####################################
		## LC only, DOME SELECTIVITY
		####################################
		sdome_lconly <- file.path(lconly_dir, "selex_dome")
		dir.create(sdome_lconly, showWarnings=FALSE)

		plist_lowdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, lwa=lwa_toUse, lwb=lwb_toUse, M=M_toUse, S50=34.6, selex_input="length", M50=34.6, maturity_input="length", binwidth=2, nseasons=1, selex_type="dome", dome_sd=50)
		plist_highdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, lwa=lwa_toUse, lwb=lwb_toUse, M=M_toUse, S50=34.6, selex_input="length", M50=34.6, maturity_input="length", binwidth=2, nseasons=1, selex_type="dome", dome_sd=25)

			plot(plist$S_l, type="l", lwd=3)
			lines(plist_lowdome$S_l, col="goldenrod", lwd=3, lty=2)
			lines(plist_highdome$S_l, col="purple", lwd=3, lty=3)

			##############################################
			## LC only, DOME SELECTIVITY - LOW DOME DEFAULT
			##############################################
			lowdome_lconly_default <- file.path(sdome_lconly, "low_dome_default")
			dir.create(lowdome_lconly_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=lowdome_lconly_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_lowdome$S_l)
				
				Report <- readRDS(file.path(lowdome_lconly_default, "Report.rds"))
				Sdreport <- readRDS(file.path(lowdome_lconly_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(lowdome_lconly_default, "Inputs.rds"))

				png(file.path(lowdome_lconly_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(lowdome_lconly_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()

			##############################################
			## LC only, DOME SELECTIVITY - HIGH DOME DEFAULT
			##############################################
			highdome_lconly_default <- file.path(sdome_lconly, "high_dome_default")
			dir.create(highdome_lconly_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=highdome_lconly_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_highdome$S_l)
				
				Report <- readRDS(file.path(highdome_lconly_default, "Report.rds"))
				Sdreport <- readRDS(file.path(highdome_lconly_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(highdome_lconly_default, "Inputs.rds"))

				png(file.path(highdome_lconly_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(highdome_lconly_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()

####################################
## LC +Index
####################################

lcindex_dir <- file.path(res_dir, "lc_index")
dir.create(lcindex_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LC))
years_t <- min(years_o):max(years_o)

input_data <- list("years"=years_t, "LF"=LC, "I_t"=I_t)

		####################################
		## LC + Index, LOGISTIC SELECTIVITY
		####################################
		slogis_lcindex <- file.path(lcindex_dir, "selex_logistic")
		dir.create(slogis_lcindex, showWarnings=FALSE)

			##############################################
			## LC + Index, LOGISTIC SELECTIVITY - DEFAULT
			##############################################
			slogis_lcindex_default <- file.path(slogis_lcindex, "default")
			dir.create(slogis_lcindex_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=slogis_lcindex_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="Index_LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2))
				
				Report <- readRDS(file.path(slogis_lcindex_default, "Report.rds"))
				Sdreport <- readRDS(file.path(slogis_lcindex_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(slogis_lcindex_default, "Inputs.rds"))

				png(file.path(slogis_lcindex_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(slogis_lcindex_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()


		####################################
		## LC + Index, DOME SELECTIVITY
		####################################
		sdome_lcindex <- file.path(lcindex_dir, "selex_dome")
		dir.create(sdome_lcindex, showWarnings=FALSE)

		plist_lowdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, lwa=lwa_toUse, lwb=lwb_toUse, M=M_toUse, S50=34.6, selex_input="length", M50=34.6, maturity_input="length", binwidth=2, nseasons=1, selex_type="dome", dome_sd=50)
		plist_highdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, lwa=lwa_toUse, lwb=lwb_toUse, M=M_toUse, S50=34.6, selex_input="length", M50=34.6, maturity_input="length", binwidth=2, nseasons=1, selex_type="dome", dome_sd=25)

			plot(plist$S_l, type="l", lwd=3)
			lines(plist_lowdome$S_l, col="goldenrod", lwd=3, lty=2)
			lines(plist_highdome$S_l, col="purple", lwd=3, lty=3)

			##############################################
			## LC + Index, DOME SELECTIVITY - LOW DOME DEFAULT
			##############################################
			lowdome_lcindex_default <- file.path(sdome_lcindex, "low_dome_default")
			dir.create(lowdome_lcindex_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=lowdome_lcindex_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="Index_LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_lowdome$S_l)
				
				Report <- readRDS(file.path(lowdome_lcindex_default, "Report.rds"))
				Sdreport <- readRDS(file.path(lowdome_lcindex_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(lowdome_lcindex_default, "Inputs.rds"))

				png(file.path(lowdome_lcindex_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(lowdome_lcindex_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()

			##############################################
			## LC + Index, DOME SELECTIVITY - HIGH DOME DEFAULT
			##############################################
			highdome_lcindex_default <- file.path(sdome_lcindex, "high_dome_default")
			dir.create(highdome_lcindex_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=highdome_lcindex_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="Index_LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_highdome$S_l)
				
				Report <- readRDS(file.path(highdome_lcindex_default, "Report.rds"))
				Sdreport <- readRDS(file.path(highdome_lcindex_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(highdome_lcindex_default, "Inputs.rds"))

				png(file.path(highdome_lcindex_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(highdome_lcindex_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()
