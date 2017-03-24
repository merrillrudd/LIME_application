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

############### weighted length comps ##########################
### weighted length comp - counts with index
out <- file.path(res_dir, "Index_LC_counts")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, I_t=I_t, data_avail="Index_LC", rerun=FALSE)

### weighted length comp - counts only
out <- file.path(res_dir, "LC_counts")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, data_avail="LC", rerun=FALSE)
base_rep <- readRDS(file.path(res_dir, "LC_counts", "Report.rds"))
df <- readRDS(file.path(res_dir, "LC_counts", "check_convergence.rds"))
base_sdrep <- readRDS(file.path(res_dir, "LC_counts", "Sdreport.rds"))
base_inputs <- readRDS(file.path(res_dir, "LC_counts", "Inputs.rds"))
	### base model doesn't converge with abundance index
read_sdreport(summary(base_sdrep)[which(rownames(summary(base_sdrep))=="SPR_t"),])
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=base_inputs$Data$n_a), Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=base_inputs$Data$n_a), Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, ref=0.3)$root, error=function(e) NA)


############### weighted length comps w fixed selectivity ##########################
## weighted length comps - counts with index - fixed selectivity
out <- file.path(res_dir, "Index_LC_counts_Sfixed")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, I_t=I_t, data_avail="Index_LC")
df <- readRDS(file.path(res_dir, "Index_LC_counts_Sfixed", "check_convergence.rds"))
Report <- readRDS(file.path(res_dir, "Index_LC_counts_Sfixed", "Report.rds"))
	### base model doesn't converge with abundance index - doesn't add any information

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



############### weighted length comps w fixed selectivity & fixed SigmaR ##########################
## weighted length comps - counts - fixed selectivity -- fix SigmaR
out <- file.path(res_dir, "LC_counts_Sfixed_SigmaRfixed")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, est_sigma=NULL, param_adjust="SigmaR", val_adjust=0.7)

## weighted length comps - counts - fixed selectivity -- fix SigmaR low
out <- file.path(res_dir, "LC_counts_Sfixed_SigmaRfixedlow")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LC, lh=plist, S_l_input=plist$S_l, est_sigma=NULL, param_adjust="SigmaR", val_adjust=0.3)

############### weighted length comps ##########################
### weighted length comps - proportions + Index
out <- file.path(res_dir, "Index_LC_prop")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist, I_t=I_t, data_avail="Index_LC")

### weighted length comps - proportions only
out <- file.path(res_dir, "LC_prop")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop, lh=plist, data_avail="LC")
bp_rep <- readRDS(file.path(res_dir, "LC_prop", "Report.rds"))


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

run_model_set(dir=out, LC=LC_raw, lh=plist)
Report <- readRDS(file.path(res_dir, "LC_counts_unweighted", "Report.rds"))
# Report$jnll

### unweighted length comps - proportions
out <- file.path(res_dir, "LC_prop_unweighted")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop_raw, lh=plist)


############### bottom longline only ##########################
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

### bottom longline only length comps - proportions
out <- file.path(res_dir, "LC_prop_BL")
dir.create(out, showWarnings=FALSE)

run_model_set(dir=out, LC=LCprop_bline, lh=plist)







