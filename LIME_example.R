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

stock_name <- "example"
proper_name <- "Lutjanus guttatus"

###################################
## Species directories
###################################

sp_dir <- file.path(main_dir, "example")
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

## choose values (from a study or from FishBase)
linf_toUse <- linf_med
vbk_toUse <- vbk_med
t0_toUse <- -0.01

## use natural mortality tool
M_dist <- c(0.15, 0.29, 0.31, 0.50, 0.41)
M_toUse <- median(M_dist)

## add life history information for species chosen for assessment
## siganus sutor
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", binwidth=1, nseasons=1)

par(mfrow=c(2,2))
plot(plist$L_a, type="l", lwd=3)
plot(plist$W_a, type="l", lwd=3)
plot(plist$Mat_l, type="l", lwd=3)
plot(plist$S_l, type="l", lwd=3)

###################################
## Load data
###################################
example_data <- generate_data(modpath=data_dir, data_avail="LC10", Fdynamics="Ramp", Rdynamics="AR", lh=plist, write=FALSE, Nyears=10, comp_sample=200, rewrite=TRUE, init_depl=0.6, itervec=1)
saveRDS(example_data, file.path(data_dir, "example_data.rds"))

example_data <- readRDS(file.path(data_dir, "example_data.rds"))

LC <- example_data$LF
plot_LC(Inputs=list("LF"=LC))

###################################
## Run assessments
###################################

####################################
## ALL YEARS
####################################

allyears_dir <- file.path(res_dir, "all_years")
dir.create(allyears_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LC))
years_t <- min(years_o):max(years_o)

input_data <- list("years"=years_t, "LF"=LC)

		####################################
		## ALL YEARS, LBSPR
		####################################
		allyears_lbspr <- file.path(allyears_dir, "LBSPR")
		dir.create(allyears_lbspr, showWarnings=FALSE)

		run <- run_LBSPR(modpath=allyears_lbspr, lh=plist, species=stock_name, input_data=input_data, rewrite=TRUE, simulation=FALSE, write=TRUE)	
		lbspr_res <- readRDS(file.path(allyears_lbspr, "LBSPR_results.rds"))

		####################################
		## ALL YEARS, LOGISTIC SELECTIVITY
		####################################
		slogis_allyears <- file.path(allyears_dir, "selex_logistic")
		dir.create(slogis_allyears, showWarnings=FALSE)

			##############################################
			## ALL YEARS, LOGISTIC SELECTIVITY - DEFAULT
			##############################################
			slogis_allyears_default <- file.path(slogis_allyears, "default")
			dir.create(slogis_allyears_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=slogis_allyears_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2))
				
				Report <- readRDS(file.path(slogis_allyears_default, "Report.rds"))
				Sdreport <- readRDS(file.path(slogis_allyears_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(slogis_allyears_default, "Inputs.rds"))

				png(file.path(slogis_allyears_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res, True=example_data)
				dev.off()				

				png(file.path(slogis_allyears_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()


		####################################
		## ALL YEARS, DOME SELECTIVITY
		####################################
		sdome_allyears <- file.path(allyears_dir, "selex_dome")
		dir.create(sdome_allyears, showWarnings=FALSE)

			plist_lowdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", binwidth=1, nseasons=1, selex_type="dome", dome_sd=70)
			plist_highdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", binwidth=1, nseasons=1, selex_type="dome", dome_sd=40)

			plot(plist$S_l, type="l", lwd=3)
			lines(plist_lowdome$S_l, col="goldenrod", lwd=3, lty=2)
			lines(plist_highdome$S_l, col="purple", lwd=3, lty=3)

			##############################################
			## ALL YEARS, DOME SELECTIVITY - LOW DOME DEFAULT
			##############################################
			lowdome_allyears_default <- file.path(sdome_allyears, "low_dome_default")
			dir.create(lowdome_allyears_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=lowdome_allyears_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_lowdome$S_l)
				
				Report <- readRDS(file.path(lowdome_allyears_default, "Report.rds"))
				Sdreport <- readRDS(file.path(lowdome_allyears_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(lowdome_allyears_default, "Inputs.rds"))

				png(file.path(lowdome_allyears_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res, True=example_data)
				dev.off()				

				png(file.path(lowdome_allyears_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()

			##############################################
			## ALL YEARS, DOME SELECTIVITY - HIGH DOME DEFAULT
			##############################################
			highdome_allyears_default <- file.path(sdome_allyears, "high_dome_default")
			dir.create(highdome_allyears_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=highdome_allyears_default, write=TRUE, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_highdome$S_l)
				
				Report <- readRDS(file.path(highdome_allyears_default, "Report.rds"))
				Sdreport <- readRDS(file.path(highdome_allyears_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(highdome_allyears_default, "Inputs.rds"))

				png(file.path(highdome_allyears_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res, True=example_data)
				dev.off()				

				png(file.path(highdome_allyears_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()
