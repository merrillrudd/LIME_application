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

stock_name <- "costa_rican_snapper"
proper_name <- "Lutjanus guttatus"

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
## Load data
###################################
## read weighted data
LCprop <- readRDS(file.path(data_dir, "CR_Snapper_LFprop_weighted.rds"))
LC <- readRDS(file.path(data_dir, "CR_Snapper_LF_weighted.rds"))
bins <- as.numeric(colnames(LC))
years <- rownames(LC)


png(file.path(figs_dir, "LF_weighted.png"), height=10, width=14, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LC), lc_years=years, gears=NULL)
dev.off()

## weighted data by gear
LCprop_bygear <- readRDS(file.path(data_dir, "CR_Snapper_LFprop_gear_year.rds"))
gearmat <- readRDS(file.path(data_dir, "Samples_per_gear.rds"))
gears <- colnames(gearmat)

png(file.path(figs_dir, "LF_weighted_byGear.png"), height=10, width=14, res=200, units="in")
plot_LCdata(lc_list=LCprop_bygear, lc_years=years, gears=gears)
dev.off()

## raw unweighted length comp
LCprop_raw <- readRDS(file.path(data_dir, "CR_Snapper_LFprop_raw.rds"))
LC_raw <- readRDS(file.path(data_dir, "CR_Snapper_LF_raw.rds"))

png(file.path(figs_dir, "LF_compare_weighted_unweighted.png"), height=10, width=14, res=200, units="in")
plot_LCfits(Inputs2=list("LF"=LC), Inputs=list("LF"=LC_raw))
dev.off()

## LC separate by gears
LC_bline <- readRDS(file.path(data_dir, "CR_Snapper_LF_bottom_longline.rds"))
LCprop_bline <- t(sapply(1:length(years), function(x) LC_bline[x,]/sum(LC_bline[x,])))
rownames(LCprop_bline) <- rownames(LC_bline)

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

## choose values (from Bystrom thesis)
linf_toUse <- 64.58
vbk_toUse <- 0.21
t0_toUse <- -0.01
lwa_toUse <- 0.0245
lwb_toUse <- 2.790

## use natural mortality tool
M_toUse <- 0.43

## add life history information for species chosen for assessment
## siganus sutor
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwb_toUse, S50=32, S95=39, selex_input="length", M50=34, maturity_input="length", M=M_toUse, SigmaR=0.7, SigmaF=0.2, F1=0.34, CVlen=0.1, nseasons=1, binwidth=bins[1], AgeMax=16)

par(mfrow=c(2,2))
plot(plist$L_a, type="l", lwd=3)
plot(plist$W_a, type="l", lwd=3)
plot(plist$Mat_l, type="l", lwd=3)
plot(plist$S_l, type="l", lwd=3)

png(file.path(figs_dir, "LF_weighted_S50.png"), height=10, width=14, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LC), lc_years=years, gears=NULL, S50=plist$SL50)
dev.off()

###################################
## Run assessments
###################################

####################################
## LC only
####################################

############### weighted length comps ##########################
### weighted length comp - counts
out <- file.path(res_dir, "LC_counts")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist)
Report <- readRDS(file.path(res_dir, "LC_counts", "Report.rds"))

## weighted length comps - counts - fixed selectivity
out <- file.path(res_dir, "LC_counts_Sfixed")
dir.create(out, showWarnings=FALSE)

# run_model_set(dir=out, LC=LC, lh=plist, fix_param="logS50")
run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l)
Report <- readRDS(file.path(res_dir, "LC_counts_Sfixed", "Report.rds"))

## weighted length comps - counts - fixed selectivity -- lower sigmaR
out <- file.path(res_dir, "LC_counts_Sfixed_lowerSigmaR")
dir.create(out, showWarnings=FALSE)

plist_lowerSigR <- with(plist, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, S95=SL95, M95=ML95, selex_input="length", maturity_input="length", AgeMax=AgeMax, F1=F1, nseasons=1, SigmaR=0.3, SigmaF=SigmaF))

run_model_set(dir=out, LC=LC, lh=plist_lowerSigR, S_l_input=plist_lowerSigR$S_l, est_sigma=NULL)

### weighted length comps - proportions
out <- file.path(res_dir, "LC_prop")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist)

### weighted length comps - proportions - fixed selectivity
out <- file.path(res_dir, "LC_prop_Sfixed")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist, S_l_input=plist$S_l)



## weighted length comps - dome-shaped selex
plist_lowdome <- with(plist, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, S95=SL95, M95=ML95, selex_input="length", maturity_input="length", selex_type="dome", dome_sd=150, AgeMax=AgeMax, F1=F1, nseasons=1, SigmaR=SigmaR, SigmaF=SigmaF))
plist_highdome <- with(plist, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, S95=SL95, M95=ML95, selex_input="length", maturity_input="length", selex_type="dome", dome_sd=40, AgeMax=AgeMax, F1=F1, nseasons=1, SigmaR=SigmaR, SigmaF=SigmaF))
par(mfrow=c(1,1))
plot(plist$S_l, type="l", lwd=3)
lines(plist_lowdome$S_l, col="goldenrod", lwd=3, lty=2)
lines(plist_highdome$S_l, col="purple", lwd=3, lty=3)



## weighted length comps - counts - low dome
out <- file.path(res_dir, "LC_counts_lowdome")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist_lowdome$S_l)

## weighted length comps - counts - high dome
out <- file.path(res_dir, "LC_counts_highdome")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist_highdome$S_l)


############### unweighted length comps ##########################
### unweighted length comp - counts
out <- file.path(res_dir, "LC_counts_unweighted")
dir.create(out, showWarnings=FALSE)

# run_model_set(dir=out, LC=LC_raw, lh=plist)
Report <- readRDS(file.path(out, "Report.rds"))
Report$jnll

### unweighted length comps - proportions
out <- file.path(res_dir, "LC_prop_unweighted")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop_raw, lh=plist)


############### bottom longline only ##########################
### bottom longline only length comp - counts
out <- file.path(res_dir, "LC_counts_BL")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC_bline, lh=plist)

### bottom longline only length comps - proportions
out <- file.path(res_dir, "LC_prop_BL")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop_bline, lh=plist)








		####################################
		## LC only, DOME SELECTIVITY
		####################################
		sdome_lconly <- file.path(lconly_dir, "selex_dome")
		dir.create(sdome_lconly, showWarnings=FALSE)

		plist_lowdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwb_toUse, S50=lbspr_res$SL50[length(lbspr_res$SL50)], S95=lbspr_res$SL95[length(lbspr_res$SL95)], selex_input="length", M50=34, maturity_input="length", SigmaR=0.3, M=M_toUse, F1=0.34, CVlen=0.1, nseasons=1, binwidth=bins[1], selex_type="dome", dome_sd=60)
		plist_highdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwb_toUse, S50=lbspr_res$SL50[length(lbspr_res$SL50)], S95=lbspr_res$SL95[length(lbspr_res$SL95)],  selex_input="length", M50=34, maturity_input="length", SigmaR=0.3, M=M_toUse, F1=0.34, CVlen=0.1, nseasons=1, binwidth=bins[1], selex_type="dome", dome_sd=30)

			plot(plist$S_l, type="l", lwd=3)
			lines(plist_lowdome$S_l, col="goldenrod", lwd=3, lty=2)
			lines(plist_highdome$S_l, col="purple", lwd=3, lty=3)

			##############################################
			## LC only, DOME SELECTIVITY - LOW DOME DEFAULT
			##############################################
			lowdome_lconly_default <- file.path(sdome_lconly, "low_dome_default")
			dir.create(lowdome_lconly_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=lowdome_lconly_default, write=TRUE, lh=plist_lowdome, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_lowdome$S_l)
				
				Report <- readRDS(file.path(lowdome_lconly_default, "Report.rds"))
				Sdreport <- readRDS(file.path(lowdome_lconly_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(lowdome_lconly_default, "Inputs.rds"))

				png(file.path(lowdome_lconly_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(lowdome_lconly_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()

			##############################################
			## LC only, DOME SELECTIVITY - HIGH DOME DEFAULT
			##############################################
			highdome_lconly_default <- file.path(sdome_lconly, "high_dome_default")
			dir.create(highdome_lconly_default, showWarnings=FALSE)	

				## run models
				res <- run_LIME(modpath=highdome_lconly_default, write=TRUE, lh=plist_highdome, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2), S_l_input=plist_highdome$S_l)
				
				Report <- readRDS(file.path(highdome_lconly_default, "Report.rds"))
				Sdreport <- readRDS(file.path(highdome_lconly_default, "Sdreport.rds"))
				Inputs <- readRDS(file.path(highdome_lconly_default, "Inputs.rds"))

				png(file.path(highdome_lconly_default, "output.png"), width=16, height=10, res=200, units="in")
				plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
				dev.off()				

				png(file.path(highdome_lconly_default, "LCfits.png"), width=16, height=10, res=200, units="in")
				plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
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
