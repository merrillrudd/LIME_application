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

## read nominal CPUE Index
I_t <- readRDS(file.path(data_dir, "Nominal_CPUE.rds"))

## read weighted data
LCprop <- readRDS(file.path(data_dir, "CR_Snapper_LFprop_weighted.rds"))
LC <- readRDS(file.path(data_dir, "CR_Snapper_LF_weighted.rds"))
LC2 <- readRDS(file.path(data_dir, "CR_Snapper_LF_weighted_2cm.rds"))
LC15 <- readRDS(file.path(data_dir, "CR_Snapper_LF_weighted_15bins.rds"))
bins <- as.numeric(colnames(LC))
bins2 <- as.numeric(colnames(LC2))
bins15 <- as.numeric(colnames(LC15))
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
LC2_bline <- readRDS(file.path(data_dir, "CR_Snapper_LF_bottom_longline_2cm.rds"))
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
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwb_toUse, S50=32, S95=39, selex_input="length", M50=34, maturity_input="length", M=M_toUse, SigmaR=0.7, SigmaF=0.2, F1=0.34, CVlen=0.1, nseasons=1, binwidth=1, AgeMax=22)

plist2 <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwb_toUse, S50=32, S95=39, selex_input="length", M50=34, maturity_input="length", M=M_toUse, SigmaR=0.7, SigmaF=0.2, F1=0.34, CVlen=0.1, nseasons=1, binwidth=2, AgeMax=22)


par(mfrow=c(2,2))
plot(plist$L_a, type="l", lwd=3)
plot(plist$W_a, type="l", lwd=3)
plot(plist$Mat_l, type="l", lwd=3)
plot(plist$S_l, type="l", lwd=3)

png(file.path(figs_dir, "LF_weighted_S50.png"), height=10, width=14, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LC), lc_years=years, gears=NULL, S50=plist$SL50)
dev.off()

dev.new()
plot_LCdata(lc_list=list("LF"=LC2), lc_years=years, gears=NULL, S50=plist$SL50, ylim=c(0,0.2))


###################################
## Run assessments
###################################

############### weighted length comps ##########################
### weighted length comp - counts only
out <- file.path(res_dir, "LC_counts")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, data_avail="LC")
base_rep <- readRDS(file.path(res_dir, "LC_counts", "Report.rds"))
df <- readRDS(file.path(res_dir, "LC_counts", "check_convergence.rds"))
base_sdrep <- readRDS(file.path(res_dir, "LC_counts", "Sdreport.rds"))
base_inputs <- readRDS(file.path(res_dir, "LC_counts", "Inputs.rds"))
	### base model doesn't converge with abundance index
read_sdreport(summary(base_sdrep)[which(rownames(summary(base_sdrep))=="SPR_t"),], log=FALSE)
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=base_inputs$Data$n_a), Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=base_inputs$Data$n_a), Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, ref=0.3)$root, error=function(e) NA)

### weighted length comp - counts only - 2cm bins
out <- file.path(res_dir, "LC_counts_2cm")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC2, lh=plist2, data_avail="LC", lc_ylim=c(0, 0.2))
base_rep2 <- readRDS(file.path(res_dir, "LC_counts_2cm", "Report.rds"))
df2 <- readRDS(file.path(res_dir, "LC_counts_2cm", "check_convergence.rds"))
base_sdrep2 <- readRDS(file.path(res_dir, "LC_counts_2cm", "Sdreport.rds"))
base_inputs2 <- readRDS(file.path(res_dir, "LC_counts_2cm", "Inputs.rds"))

### weighted length comp - counts with index
out <- file.path(res_dir, "Index_LC_counts")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, I_t=I_t, data_avail="Index_LC")


############### weighted length comps w fixed selectivity ##########################
## weighted length comps - counts - fixed selectivity
out <- file.path(res_dir, "LC_counts_Sfixed")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, data_avail="LC")
sfixed_rep <- readRDS(file.path(res_dir, "LC_counts_Sfixed", "Report.rds"))
df <- readRDS(file.path(res_dir, "LC_counts_Sfixed", "check_convergence.rds"))
sfixed_sdrep <- readRDS(file.path(res_dir, "LC_counts_Sfixed", "Sdreport.rds"))
sfixed_inputs <- readRDS(file.path(res_dir, "LC_counts_Sfixed", "Inputs.rds"))
	### sfixed model doesn't converge with abundance index
read_sdreport(summary(sfixed_sdrep)[which(rownames(summary(sfixed_sdrep))=="SPR_t"),], log=FALSE)
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=sfixed_inputs$Data$n_a), Mat_a=sfixed_rep$Mat_a, W_a=sfixed_rep$W_a, M=sfixed_rep$M, S_a=sfixed_rep$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=sfixed_inputs$Data$n_a), Mat_a=sfixed_rep$Mat_a, W_a=sfixed_rep$W_a, M=sfixed_rep$M, S_a=sfixed_rep$S_a, ref=0.3)$root, error=function(e) NA)

## weighted length comps - counts - fixed selectivity - 2cm bins
out <- file.path(res_dir, "LC_counts_Sfixed_2cm")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC2, lh=plist2, S_l_input=plist2$S_l, data_avail="LC")

############### weighted length comps w fixed selectivity, recruitment deviations off ##########################
  ## weighted length comps - counts - fixed selectivity
out <- file.path(res_dir, "LC_counts_Sfixed_RecOff")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, data_avail="LC", param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.01,0.05), est_sigma=FALSE)
recoff_rep <- readRDS(file.path(res_dir, "LC_counts_Sfixed_RecOff", "Report.rds"))
df <- readRDS(file.path(res_dir, "LC_counts_Sfixed_RecOff", "check_convergence.rds"))
recoff_sdrep <- readRDS(file.path(res_dir, "LC_counts_Sfixed_RecOff", "Sdreport.rds"))
recoff_inputs <- readRDS(file.path(res_dir, "LC_counts_Sfixed_RecOff", "Inputs.rds"))

read_sdreport(summary(recoff_sdrep)[which(rownames(summary(recoff_sdrep))=="SPR_t"),], log=FALSE)
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=recoff_inputs$Data$n_a), Mat_a=recoff_rep$Mat_a, W_a=recoff_rep$W_a, M=recoff_rep$M, S_a=recoff_rep$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=recoff_inputs$Data$n_a), Mat_a=recoff_rep$Mat_a, W_a=recoff_rep$W_a, M=recoff_rep$M, S_a=recoff_rep$S_a, ref=0.3)$root, error=function(e) NA)


############### weighted length comps w fixed selectivity & fixed SigmaR ##########################
## weighted length comps - counts - fixed selectivity -- fix SigmaR
out <- file.path(res_dir, "LC_counts_Sfixed_SigmaRfixed")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, est_sigma=NULL, param_adjust="SigmaR", val_adjust=0.7)

## weighted length comps - counts - fixed selectivity -- fix SigmaR low
out <- file.path(res_dir, "LC_counts_Sfixed_SigmaRfixedlow")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, est_sigma=NULL, param_adjust="SigmaR", val_adjust=0.3)

############### weighted length comps  - proportions ##########################
### weighted length comps - proportions only
out <- file.path(res_dir, "LC_prop")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist, data_avail="LC")
bp_rep <- readRDS(file.path(res_dir, "LC_prop", "Report.rds"))



### weighted length comps - proportions + Index
out <- file.path(res_dir, "Index_LC_prop")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist, I_t=I_t, data_avail="Index_LC")


############### weighted length comps proportions w fixed selectivity ##########################
### weighted length comps - proportions - fixed selectivity
out <- file.path(res_dir, "LC_prop_Sfixed")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist, S_l_input=plist$S_l)

############### dome-shaped selectivity ##########################
## weighted length comps - dome-shaped selex
plist_lowdome <- with(plist, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, S95=SL95, M95=ML95, selex_input="length", maturity_input="length", selex_type="dome", dome_sd=75, AgeMax=AgeMax, F1=F1, nseasons=1, SigmaR=SigmaR, SigmaF=SigmaF))
plist_highdome <- with(plist, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, S95=SL95, M95=ML95, selex_input="length", maturity_input="length", selex_type="dome", dome_sd=35, AgeMax=AgeMax, F1=F1, nseasons=1, SigmaR=SigmaR, SigmaF=SigmaF))
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


################### unweighted length comps ##########################
### unweighted length comp - counts
out <- file.path(res_dir, "LC_counts_unweighted")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC_raw, lh=plist)
Report <- readRDS(file.path(res_dir, "LC_counts_unweighted", "Report.rds"))
# Report$jnll

### unweighted length comps - proportions
out <- file.path(res_dir, "LC_prop_unweighted")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop_raw, lh=plist)


######################## bottom longline only ##########################
### bottom longline only length comp - counts
out <- file.path(res_dir, "LC_counts_BL")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC_bline, lh=plist, rerun=TRUE)
bl_rep <- readRDS(file.path(res_dir, "LC_counts_BL", "Report.rds"))
df <- readRDS(file.path(res_dir, "LC_counts_BL", "check_convergence.rds"))
bl_sdrep <- readRDS(file.path(res_dir, "LC_counts_BL", "Sdreport.rds"))
bl_inputs <- readRDS(file.path(res_dir, "LC_counts_BL", "Inputs.rds"))
bl_lbspr <- readRDS(file.path(res_dir, "LC_counts_BL", "LBSPR_results.rds"))

read_sdreport(summary(bl_sdrep)[which(rownames(summary(bl_sdrep))=="SPR_t"),], log=FALSE)
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=bl_inputs$Data$n_a), Mat_a=bl_rep$Mat_a, W_a=bl_rep$W_a, M=bl_rep$M, S_a=bl_rep$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=bl_inputs$Data$n_a), Mat_a=bl_rep$Mat_a, W_a=bl_rep$W_a, M=bl_rep$M, S_a=bl_rep$S_a, ref=0.3)$root, error=function(e) NA)

### bottom longline only length comp - counts - 2cm
out <- file.path(res_dir, "LC_counts_BL_2cm")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC2_bline, lh=plist2, rerun=TRUE)

### bottom longline only length comps - proportions
out <- file.path(res_dir, "LC_prop_BL")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop_bline, lh=plist)


############### sensitivity analyses ##########################


### natural mortality
out <- file.path(res_dir, "M_prof")
dir.create(out, showWarnings=FALSE)

M_tool <- c(0.41, 0.34, 0.36, 0.36, 0.27, 0.16, 0.32, 0.32, 0.34)
M_vec <- seq(from=min(min(M_tool), (M_toUse-0.25*M_toUse)), to=max((M_toUse+0.25*M_toUse), max(M_tool)), length.out=10)

run_prof(dir=out, vec=M_vec, param="M", LC=LC, lh=plist)

### length at 50% maturity
out <- file.path(res_dir, "ML50_prof")
dir.create(out, showWarnings=FALSE)

ML50_vec <- seq((34-0.25*34), (34+0.25*34), length.out=10)

run_prof(dir=out, vec=ML50_vec, param="ML50", LC=LC, lh=plist)

### asymptotic length
out <- file.path(res_dir, "linf_prof")
dir.create(out, showWarnings=FALSE)

linf_vec <- seq((linf_toUse-0.25*linf_toUse), (linf_toUse+0.25*linf_toUse), length.out=10)

run_prof(dir=out, vec=linf_vec, param="linf", LC=LC, lh=plist)

### growth coefficient
out <- file.path(res_dir, "vbk_prof")
dir.create(out, showWarnings=FALSE)

vbk_vec <- seq((vbk_toUse-0.25*vbk_toUse), (vbk_toUse+0.25*vbk_toUse), length.out=10)

run_prof(dir=out, vec=vbk_vec, param="vbk", LC=LC, lh=plist)

## sigma F
out <- file.path(res_dir, "SigmaF_prof")
dir.create(out, showWarnings=FALSE)

SigmaF_vec <- seq(0.01, 0.4, length.out=10)

run_prof(dir=out, vec=SigmaF_vec, param="SigmaF", LC=LC, lh=plist)

## sigma F
out <- file.path(res_dir, "SigmaF_prof_2cm")
dir.create(out, showWarnings=FALSE)

SigmaF_vec <- seq(0.01, 0.4, length.out=10)

run_prof(dir=out, vec=SigmaF_vec, param="SigmaF", LC=LC2, lh=plist2)


## sigma R
out <- file.path(res_dir, "SigmaR_prof")
dir.create(out, showWarnings=FALSE)

SigmaR_vec <- seq(0.01, 1, length.out=10)

run_prof(dir=out, vec=SigmaR_vec, param="SigmaR", LC=LC, lh=plist)


## both sigmaR and sigmaF at the same time
SigmaRF_dir <- file.path(res_dir, "SigmaRF_prof")
dir.create(SigmaRF_dir, showWarnings=FALSE)

dir <- SigmaRF_dir

for(i in 1:length(SigmaR_vec)){
	out <- file.path(dir, paste0("SigmaR_SigmaF_", i))
	dir.create(out, showWarnings=FALSE)

	param_adjust <- c("SigmaR", "SigmaF")
	val_adjust <- c(SigmaR_vec[i], SigmaF_vec[i])

	  lh_new <- plist
	  for(p in 1:length(param_adjust)){
	  	lh_new[[param_adjust[p]]] <- val_adjust[p]
	  }

	years_o <- as.numeric(rownames(LC))
	years_m <- (min(years_o)-0):max(years_o)

	input_data <- list("years"=years_m, "LF"=LC, "C_t"=NULL, "I_t"=NULL)


				run <- run_LBSPR(modpath=out, lh=lh_new, species="x", input_data=input_data, rewrite=TRUE, simulation=FALSE)	
				lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

			res <- run_LIME(modpath=out, lh=lh_new, input_data=input_data, est_sigma=FALSE, data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, S_l_input=-1, F_up=5, fix_param=FALSE, param_adjust=param_adjust, val_adjust=val_adjust)

			if(file.exists(file.path(out, "LBSPR_results.rds"))) lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))
			if(file.exists(file.path(out, "LBSPR_results.rds"))==FALSE) lbspr_res <- NULL
				
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
			dev.off()				

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
			dev.off()
	}

  ## plot likelihood profile
	nll_vec <- rep(NA, length(SigmaR_vec))
	col_vec <- rep("black", length(SigmaR_vec))
	for(i in 1:length(SigmaR_vec)){
	out <- file.path(dir, paste0("SigmaR_SigmaF_", i))
		Report <- readRDS(file.path(out, "Report.rds"))
		if(round(max(Report$F_t))==10) col_vec[i] <- "gray"
		nll_vec[i] <- Report$jnll
	}
	png(file.path(dir, "likelihood_profile.png"), width=16, height=10, res=200, units="in")
	plot(x=SigmaR_vec, y=nll_vec, xlab="SigmaR SigmaF combo", ylab="NLL", xaxs="i", yaxs="i", pch=19, cex=2, xpd=NA, ylim=c(0.95*min(nll_vec), 1.05*max(nll_vec)), col=col_vec)
	dev.off()


## both sigmaR and sigmaF at the same time
SigmaRF_dir <- file.path(res_dir, "SigmaRF_prof_2cm")
dir.create(SigmaRF_dir, showWarnings=FALSE)

dir <- SigmaRF_dir

for(i in 1:length(SigmaR_vec)){
	out <- file.path(dir, paste0("SigmaR_SigmaF_", i))
	dir.create(out, showWarnings=FALSE)

	param_adjust <- c("SigmaR", "SigmaF")
	val_adjust <- c(SigmaR_vec[i], SigmaF_vec[i])

	  lh_new <- plist2
	  for(p in 1:length(param_adjust)){
	  	lh_new[[param_adjust[p]]] <- val_adjust[p]
	  }

	years_o <- as.numeric(rownames(LC))
	years_m <- (min(years_o)-0):max(years_o)

	input_data <- list("years"=years_m, "LF"=LC2, "C_t"=NULL, "I_t"=NULL)


				run <- run_LBSPR(modpath=out, lh=lh_new, species="x", input_data=input_data, rewrite=TRUE, simulation=FALSE)	
				lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

			res <- run_LIME(modpath=out, lh=lh_new, input_data=input_data, est_sigma=FALSE, data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, S_l_input=-1, F_up=5, fix_param=FALSE, param_adjust=param_adjust, val_adjust=val_adjust)

			if(file.exists(file.path(out, "LBSPR_results.rds"))) lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))
			if(file.exists(file.path(out, "LBSPR_results.rds"))==FALSE) lbspr_res <- NULL
				
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
			dev.off()				

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res, ylim=c(0,0.2))
			dev.off()
	}

  ## plot likelihood profile
	nll_vec <- rep(NA, length(SigmaR_vec))
	col_vec <- rep("black", length(SigmaR_vec))
	for(i in 1:length(SigmaR_vec)){
	out <- file.path(dir, paste0("SigmaR_SigmaF_", i))
		Report <- readRDS(file.path(out, "Report.rds"))
		if(round(max(Report$F_t))==10) col_vec[i] <- "gray"
		nll_vec[i] <- Report$jnll
	}
	png(file.path(dir, "likelihood_profile.png"), width=16, height=10, res=200, units="in")
	plot(x=SigmaR_vec, y=nll_vec, xlab="SigmaR SigmaF combo", ylab="NLL", xaxs="i", yaxs="i", pch=19, cex=2, xpd=NA, ylim=c(0.95*min(nll_vec), 1.05*max(nll_vec)), col=col_vec)
	dev.off()







## summarize sensitivity analyses
sens_param <- c("M", "linf", "vbk", "ML50")
sens_vec <- list()
sens_vec[["M"]] <- M_vec
sens_vec[["linf"]] <- linf_vec
sens_vec[["vbk"]] <- vbk_vec
sens_vec[["ML50"]] <- ML50_vec

all_years <- min(years):max(years)
lc_years <- years

png(file.path(figs_dir, "Compare_sensitivities.png"), height=16, width=24, res=200, units="in")
ramp <- colorRamp(c("red", "orange"))
col_vec <- rgb( ramp(seq(0, 1, length = 5)), max = 255)
nf <- layout(matrix(c(1,6,11,16,
				1,6,11,16,
				2,7,12,17,
				3,8,13,18,
				3,8,13,18,
				3,8,13,18,
				4,9,14,19,
				4,9,14,19,
				4,9,14,19,
				5,10,15,20,
				5,10,15,20,
				5,10,15,20), nrow=12, ncol=4, byrow=TRUE))
layout.show(nf)
par(mar=c(0,3,0,0), omi=c(1,1,1,1))
for(p in 1:length(sens_param)){
	base <- file.path(res_dir, "LC_counts")
	base_rep <- readRDS(file.path(base, "Report.rds"))
	base_sdrep <- readRDS(file.path(base, "Sdreport.rds"))

	out <- file.path(res_dir, paste0(sens_param[p], "_prof"))
	nll_vec <- rep(NA, length(sens_vec[[p]]))
	for(i in 1:length(sens_vec[[p]])){
		dir <- file.path(out, paste0(sens_param[p], "_", sens_vec[[p]][i]))
		Report <- readRDS(file.path(dir, "Report.rds"))
		nll_vec[i] <- Report$jnll
	}
	plot(x=sens_vec[[p]], y=nll_vec, pch=19, cex=2, xlab="", ylab=sens_param[p], ylim=c(0.9*min(nll_vec), 1.1*max(nll_vec)), col=c(col_vec[1], rep("black", length(sens_vec[[p]])-2,), col_vec[5]))
	points(x=plist[[sens_param[p]]], y=base_rep$jnll, pch=19, cex=2, col="blue")
	mtext(side=3, sens_param[p], font=2, line=2)
	if(p==1) mtext(side=2, "NLL", font=2, line=3)

	plot(x=1,y=1,type="n",axes=F,ann=F)

	Report1 <- readRDS(file.path(out, paste0(sens_param[p], "_", sens_vec[[p]][1]), "Report.rds"))
	Report2 <- readRDS(file.path(out, paste0(sens_param[p], "_", sens_vec[[p]][length(sens_vec[[p]])]), "Report.rds"))
	Sdreport1 <- readRDS(file.path(out, paste0(sens_param[p], "_", sens_vec[[p]][1]), "Sdreport.rds"))
	Sdreport2 <- readRDS(file.path(out, paste0(sens_param[p], "_", sens_vec[[p]][length(sens_vec[[p]])]),"Sdreport.rds"))

	ylim <- c(0, 1.5)
	plot(x=seq_along(all_years), y=Report1$F_t, lwd=2, col=col_vec[1], ylim=ylim, type="l", xaxt="n", ylab="", xlab="Year", xaxs="i", yaxs="i", cex.lab=2)
	points(x=which(all_years %in% lc_years), y=Report1$F_t[which(all_years %in% lc_years)], col=col_vec[1], pch=19, xpd=NA)
	sd <- summary(Sdreport1)[which(rownames(summary(Sdreport1))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[1],50), border=NA)

  	lines(x=seq_along(all_years), y=Report2$F_t, lwd=2, col=col_vec[5])
  	points(x=which(all_years %in% lc_years), y=Report2$F_t[which(all_years %in% lc_years)], col=col_vec[5], pch=19, xpd=NA)
 	sd <- summary(Sdreport2)[which(rownames(summary(Sdreport2))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[5],50), border=NA)


  	lines(x=seq_along(all_years), y=base_rep$F_t, lwd=2, col="#0000AA")
  	points(x=which(all_years %in% lc_years), y=base_rep$F_t[which(all_years %in% lc_years)], col="#0000AA", pch=19, xpd=NA)
 	sd <- summary(base_sdrep)[which(rownames(summary(base_sdrep))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0("#0000AA",50), border=NA)
  	if(p==1) mtext(side=2, "Fishing mortality", font=2, line=3)

	ylim <- c(0, 2.2)
	plot(x=seq_along(all_years), y=Report1$R_t, lwd=2, col=col_vec[1], ylim=ylim, type="l", xaxt="n", ylab="", xlab="Year", xaxs="i", yaxs="i", cex.lab=2)
	points(x=which(all_years %in% lc_years), y=Report1$R_t[which(all_years %in% lc_years)], col=col_vec[1], pch=19, xpd=NA)
	sd <- summary(Sdreport1)[which(rownames(summary(Sdreport1))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[1],50), border=NA)

  	lines(x=seq_along(all_years), y=Report2$R_t, lwd=2, col=col_vec[5])
  	points(x=which(all_years %in% lc_years), y=Report2$R_t[which(all_years %in% lc_years)], col=col_vec[5], pch=19, xpd=NA)
 	sd <- summary(Sdreport2)[which(rownames(summary(Sdreport2))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[5],50), border=NA)


  	lines(x=seq_along(all_years), y=base_rep$R_t, lwd=2, col="#0000AA")
  	points(x=which(all_years %in% lc_years), y=base_rep$R_t[which(all_years %in% lc_years)], col="#0000AA", pch=19, xpd=NA)
 	sd <- summary(base_sdrep)[which(rownames(summary(base_sdrep))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0("#0000AA",50), border=NA)
  	if(p==1) mtext(side=2, "Recruitment", font=2, line=3)

	ylim <- c(0, 1.2)
	plot(x=seq_along(all_years), y=Report1$SPR_t, lwd=2, col=col_vec[1], ylim=ylim, type="l", xaxt="n", ylab="", xlab="Year", xaxs="i", yaxs="i", cex.lab=2)
	points(x=which(all_years %in% lc_years), y=Report1$SPR_t[which(all_years %in% lc_years)], col=col_vec[1], pch=19, xpd=NA)
	sd <- summary(Sdreport1)[which(rownames(summary(Sdreport1))=="lSPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[1],50), border=NA)
  	axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))

  	lines(x=seq_along(all_years), y=Report2$SPR_t, lwd=2, col=col_vec[5])
  	points(x=which(all_years %in% lc_years), y=Report2$SPR_t[which(all_years %in% lc_years)], col=col_vec[5], pch=19, xpd=NA)
 	sd <- summary(Sdreport2)[which(rownames(summary(Sdreport2))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[5],50), border=NA)


  	lines(x=seq_along(all_years), y=base_rep$SPR_t, lwd=2, col="#0000AA")
  	points(x=which(all_years %in% lc_years), y=base_rep$SPR_t[which(all_years %in% lc_years)], col="#0000AA", pch=19, xpd=NA)
 	sd <- summary(base_sdrep)[which(rownames(summary(base_sdrep))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0("#0000AA",50), border=NA)
  	if(p==1) mtext(side=2, "SPR", font=2, line=3)
}
mtext(side=1, outer=TRUE, "Year", font=2, line=4)
dev.off()




sens_var <- c("SigmaF", "SigmaR", "SigmaR_SigmaF")
var_vec <- list()
var_vec[["SigmaF"]] <- SigmaF_vec
var_vec[["SigmaR"]] <- SigmaR_vec

png(file.path(figs_dir, "Compare_variance_sensitivities.png"), height=16, width=24, res=200, units="in")
ramp <- colorRamp(c("red", "orange"))
col_vec <- rgb( ramp(seq(0, 1, length = 5)), max = 255)
nf <- layout(matrix(c(1,6,11,
				1,6,11,
				2,7,12,
				3,8,13,
				3,8,13,
				3,8,13,
				4,9,14,
				4,9,14,
				4,9,14,
				5,10,15,
				5,10,15,
				5,10,15), nrow=12, ncol=3, byrow=TRUE))
layout.show(nf)
par(mar=c(0,3,0,0), omi=c(1,1,1,1))
for(p in 1:length(sens_var)){
	base <- file.path(res_dir, "LC_counts")
	base_rep <- readRDS(file.path(base, "Report.rds"))
	base_sdrep <- readRDS(file.path(base, "Sdreport.rds"))

	if(p<3) out <- file.path(res_dir, paste0(sens_var[p], "_prof"))
	if(p==3) out <- file.path(res_dir, paste0("SigmaRF_prof"))
	if(p<3) nll_vec <- rep(NA, length(var_vec[[p]]))
	if(p==3) nll_vec <- rep(NA, length(var_vec[[p-1]]))
	for(i in 1:length(nll_vec)){
		if(p<3) dir <- file.path(out, paste0(sens_var[p], "_", var_vec[[p]][i]))
		if(p==3) dir <- file.path(out, paste0(sens_var[p], "_", i))
		Report <- readRDS(file.path(dir, "Report.rds"))
		nll_vec[i] <- Report$jnll
	}
	if(p<3) plot(x=var_vec[[p]], y=nll_vec, pch=19, cex=2, xlab="", ylab=sens_var[p], ylim=c(0.9*min(nll_vec), 1.1*max(nll_vec)), col=c(col_vec[1], rep("black", length(nll_vec)-2,), col_vec[5]))
	if(p==3){
		plot(x=1:length(nll_vec), y=nll_vec, pch=19, cex=2, xlab="", ylab=sens_var[p], ylim=c(0.9*min(nll_vec), 1.1*max(nll_vec)), col=c(col_vec[1], rep("black", length(nll_vec)-2,), col_vec[5]), xaxt="n")
		axis(1, at=c(1, 10), labels=c("Low variability", "High variability"))
	}
	if(p<3) points(x=plist[[sens_var[p]]], y=base_rep$jnll, pch=19, cex=2, col="blue")
	if(p<3) mtext(side=3, sens_var[p], font=2, line=2)
	if(p==3) mtext(side=3, "Both", font=2, line=2)
	if(p==1) mtext(side=2, "NLL", font=2, line=3)

	plot(x=1,y=1,type="n",axes=F,ann=F)

	if(p<3){
		Report1 <- readRDS(file.path(out, paste0(sens_var[p], "_", var_vec[[p]][1]), "Report.rds"))
		Report2 <- readRDS(file.path(out, paste0(sens_var[p], "_", var_vec[[p]][length(var_vec[[p]])]), "Report.rds"))
		Sdreport1 <- readRDS(file.path(out, paste0(sens_var[p], "_", var_vec[[p]][1]), "Sdreport.rds"))
		Sdreport2 <- readRDS(file.path(out, paste0(sens_var[p], "_", var_vec[[p]][length(var_vec[[p]])]),"Sdreport.rds"))
	}
	if(p==3){
		Report1 <- readRDS(file.path(out, paste0("SigmaR_SigmaF_", 1), "Report.rds"))
		Report2 <- readRDS(file.path(out, paste0("SigmaR_SigmaF_", length(nll_vec)), "Report.rds"))
		Sdreport1 <- readRDS(file.path(out, paste0("SigmaR_SigmaF_", 1), "Sdreport.rds"))
		Sdreport2 <- readRDS(file.path(out, paste0("SigmaR_SigmaF_", length(nll_vec)),"Sdreport.rds"))
	}

	ylim <- c(0, 1.5)
	plot(x=seq_along(all_years), y=Report1$F_t, lwd=2, col=col_vec[1], ylim=ylim, type="l", xaxt="n", ylab="", xlab="Year", xaxs="i", yaxs="i", cex.lab=2)
	points(x=which(all_years %in% lc_years), y=Report1$F_t[which(all_years %in% lc_years)], col=col_vec[1], pch=19, xpd=NA)
	sd <- summary(Sdreport1)[which(rownames(summary(Sdreport1))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[1],50), border=NA)

  	lines(x=seq_along(all_years), y=Report2$F_t, lwd=2, col=col_vec[5])
  	points(x=which(all_years %in% lc_years), y=Report2$F_t[which(all_years %in% lc_years)], col=col_vec[5], pch=19, xpd=NA)
 	sd <- summary(Sdreport2)[which(rownames(summary(Sdreport2))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[5],50), border=NA)


  	lines(x=seq_along(all_years), y=base_rep$F_t, lwd=2, col="#0000AA")
  	points(x=which(all_years %in% lc_years), y=base_rep$F_t[which(all_years %in% lc_years)], col="#0000AA", pch=19, xpd=NA)
 	sd <- summary(base_sdrep)[which(rownames(summary(base_sdrep))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0("#0000AA",50), border=NA)
  	if(p==1) mtext(side=2, "Fishing mortality", font=2, line=3)

	ylim <- c(0, 2.2)
	plot(x=seq_along(all_years), y=Report1$R_t, lwd=2, col=col_vec[1], ylim=ylim, type="l", xaxt="n", ylab="", xlab="Year", xaxs="i", yaxs="i", cex.lab=2)
	points(x=which(all_years %in% lc_years), y=Report1$R_t[which(all_years %in% lc_years)], col=col_vec[1], pch=19, xpd=NA)
	sd <- summary(Sdreport1)[which(rownames(summary(Sdreport1))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[1],50), border=NA)

  	lines(x=seq_along(all_years), y=Report2$R_t, lwd=2, col=col_vec[5])
  	points(x=which(all_years %in% lc_years), y=Report2$R_t[which(all_years %in% lc_years)], col=col_vec[5], pch=19, xpd=NA)
 	sd <- summary(Sdreport2)[which(rownames(summary(Sdreport2))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[5],50), border=NA)


  	lines(x=seq_along(all_years), y=base_rep$R_t, lwd=2, col="#0000AA")
  	points(x=which(all_years %in% lc_years), y=base_rep$R_t[which(all_years %in% lc_years)], col="#0000AA", pch=19, xpd=NA)
 	sd <- summary(base_sdrep)[which(rownames(summary(base_sdrep))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0("#0000AA",50), border=NA)
  	if(p==1) mtext(side=2, "Recruitment", font=2, line=3)

	ylim <- c(0, 1.2)
	plot(x=seq_along(all_years), y=Report1$SPR_t, lwd=2, col=col_vec[1], ylim=ylim, type="l", xaxt="n", ylab="", xlab="Year", xaxs="i", yaxs="i", cex.lab=2)
	points(x=which(all_years %in% lc_years), y=Report1$SPR_t[which(all_years %in% lc_years)], col=col_vec[1], pch=19, xpd=NA)
	sd <- summary(Sdreport1)[which(rownames(summary(Sdreport1))=="lSPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[1],50), border=NA)
  	axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))

  	lines(x=seq_along(all_years), y=Report2$SPR_t, lwd=2, col=col_vec[5])
  	points(x=which(all_years %in% lc_years), y=Report2$SPR_t[which(all_years %in% lc_years)], col=col_vec[5], pch=19, xpd=NA)
 	sd <- summary(Sdreport2)[which(rownames(summary(Sdreport2))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_vec[5],50), border=NA)


  	lines(x=seq_along(all_years), y=base_rep$SPR_t, lwd=2, col="#0000AA")
  	points(x=which(all_years %in% lc_years), y=base_rep$SPR_t[which(all_years %in% lc_years)], col="#0000AA", pch=19, xpd=NA)
 	sd <- summary(base_sdrep)[which(rownames(summary(base_sdrep))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0("#0000AA",50), border=NA)
  	if(p==1) mtext(side=2, "SPR", font=2, line=3)
}
mtext(side=1, outer=TRUE, "Year", font=2, line=4)
dev.off()


############### retrospective analysis ##########################

out <- file.path(res_dir, "retrospective")
dir.create(out, showWarnings=FALSE)

years <- as.numeric(rownames(LC))

for(y in 1:length(years)){
	if(y==1){
		LCsub <- t(as.matrix(LC[1,]))
		rownames(LCsub) <- years[1]
	}
	if(y>1){
		LCsub <- LC[1:y,]
	}

	out2 <- file.path(out, y)
	dir.create(out2, showWarnings=FALSE)

	run_model_set(dir=out2, LC=LCsub, lh=plist, data_avail="LC", plot=FALSE)
}


png(file.path(figs_dir, "Retrospective.png"), height=7, width=16, res=200, units="in")
col_vec <- brewer.pal((length(years)), "Blues")
par(mfrow=c(1,3))
plot(x=years,y=rep(1,length(years)),type="n", ylim=c(0,2), xlim=c(min(years), max(years)), xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality")
for(y in 2:length(years)){
	out2 <- file.path(out, y)
	Report <- readRDS(file.path(out2, "Report.rds"))
	Sdreport <- readRDS(file.path(out2, "Sdreport.rds"))
	lbspr <- readRDS(file.path(out2, "LBSPR_results.rds"))

	lines(x=years[1:y], y=Report$F_t, col=col_vec[y], pch=19, xpd=NA, lwd=2)
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(years[1:y], rev(years[1:y])), col=paste0(col_vec[y],50), border=NA)
}

plot(x=years,y=rep(1,length(years)),type="n", ylim=c(0,3), xlim=c(min(years), max(years)),  xaxs="i", yaxs="i", xlab="Year", ylab="Recruitment")
for(y in 2:length(years)){
	out2 <- file.path(out, y)
	Report <- readRDS(file.path(out2, "Report.rds"))
	Sdreport <- readRDS(file.path(out2, "Sdreport.rds"))
	lbspr <- readRDS(file.path(out2, "LBSPR_results.rds"))

	lines(x=years[1:y], y=Report$R_t, col=col_vec[y], pch=19, xpd=NA, lwd=2)
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=TRUE), x=c(years[1:y], rev(years[1:y])), col=paste0(col_vec[y],50), border=NA)
}


plot(x=years,y=rep(1,length(years)),type="n", ylim=c(0,1.1), xlim=c(min(years), max(years)),  xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
for(y in 2:length(years)){
	out2 <- file.path(out, y)
	Report <- readRDS(file.path(out2, "Report.rds"))
	Sdreport <- readRDS(file.path(out2, "Sdreport.rds"))
	lbspr <- readRDS(file.path(out2, "LBSPR_results.rds"))

	lines(x=years[1:y], y=Report$SPR_t, col=col_vec[y], pch=19, xpd=NA, lwd=2)
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
  	polygon( y=read_sdreport(sd, log=FALSE), x=c(years[1:y], rev(years[1:y])), col=paste0(col_vec[y],50), border=NA)
}
dev.off()




