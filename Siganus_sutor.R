rm(list=ls())

## load LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

## load LBSPR package
library(LBSPR)

## load Rfishbase
library(rfishbase)
library(gplots)


###################################
## Directories
###################################
main_dir <- "C:\\Git_Projects\\LIME_application"

R_dir <- file.path(main_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

region_dir <- file.path(main_dir, "kenyan_reef_fish")

data_dir <- file.path(region_dir, "data")
dir.create(data_dir, showWarnings=FALSE)


######################################
## Species life history
######################################

plist <- readRDS(file.path(data_dir, "Siganus_sutor_life_history_annual.rds"))
plist_m <- readRDS(file.path(data_dir, "Siganus_sutor_life_history_monthly.rds"))

## M/K
plist$M/plist$vbk

######################################
## Load data
######################################

LCraw <- readRDS(file.path(data_dir, "Siganus_sutor_LCraw.rds"))
LCm_raw <- readRDS(file.path(data_dir, "Siganus_sutor_LCraw_monthly.rds"))

LCweight <- readRDS(file.path(data_dir, "Siganus_sutor_LC_weighted.rds"))
LCm_weight <- readRDS(file.path(data_dir, "Siganus_sutor_LCmonthly_weighted.rds"))

LCm_bt <- readRDS(file.path(data_dir, "Siganus_sutor_LCmonthly_baskettrap.rds"))
LCm_highN <- readRDS(file.path(data_dir, "Siganus_sutor_LCmonthly_weighted_Nup.rds"))

####################################
## run assessments
####################################

res_dir <- file.path(stock_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LCraw))
years_t <- min(years_o):max(years_o)

months_o <- as.numeric(rownames(LCm_raw))
months_t <- min(months_o):max(months_o)

input_data_y <- list("years"=years_t, "LF"=lc_weight)
input_data_m <- list("years"=months_t, "LF"=lc_weight_mo)
input_data_bt <- list("years"=months_t, "LF"=LCm_bt)
input_data_up <- list("years"=months_t, "LF"=lc_weight_mo[which(rownames(lc_weight_mo) %in% months_up),])

############### weighted length comps ##########################
## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))


## annual LIME
out <- file.path(res_dir, "LCy")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_y$years, lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, ylim=c(0,0.3))
			dev.off()


### weighted length comp - counts only - all months
out <- file.path(res_dir, "LCm")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))


			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - all months - basket trap only
out <- file.path(res_dir, "LCm_baskettrap")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_bt, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_bt$years, lc_years=rownames(input_data_bt$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_bt$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCm_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_up, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


## dome-shaped selex
plist_m_lowdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=41)
plist_m_highdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=22)


par(mfrow=c(1,1), mar=c(5,5,4,4))
plot(plist_m$S_l, type="o", lwd=2, pch=19, cex.axis=2, xlab="Length (cm)", ylab="Proportion", cex.lab=2)
lines(plist_m_lowdome$S_l, type="o", lwd=2, lty=2, col="steelblue", pch=19)
lines(plist_m_highdome$S_l, type="o", lwd=2, lty=3, col="tomato", pch=19)

## low-dome
out <- file.path(res_dir, "LCm_lowdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=plist_m_lowdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


## low-dome
out <- file.path(res_dir, "LCm_highdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=plist_m_highdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()



### likelihood profiles
## check fishbase
growth <- popgrowth(proper_name)
# linf_vec <- seq(min(min(growth$Loo),0.75*linf_toUse),max(max(growth$Loo),1.25*linf_toUse),length.out=10)
linf_vec <- seq(0.75*linf_toUse, 1.25*linf_toUse, length.out=10)
# vbk_vec <- seq(min(min(growth$K),0.75*vbk_toUse), max(max(growth$K),1.25*vbk_toUse), length.out=10)
vbk_vec <- seq(0.75*vbk_toUse, 1.25*vbk_toUse, length.out=10)
# M_vec <- seq(min(min(M_dist),0.75*M_toUse), max(max(M_dist),1.25*M_toUse), length.out=10)
M_vec <- seq(0.75*M_toUse, 1.25*M_toUse, length.out=10)
ML50_vec <- seq(0.75*plist$ML50, 1.25*plist$ML50, length.out=10)


## linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(linf_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_toUse, linf=linf_vec[i], t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
linf_nll <- rep(NA, length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	linf_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "linf_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(linf_vec, linf_nll, pch=19, cex=2, xlab="Assumed Linf", ylab="NLL", ylim=c(0.8*min(linf_nll), max(linf_nll)*1.2))
points(linf_toUse, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(vbk_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_vec[i], linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
vbk_nll <- rep(NA, length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	vbk_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "vbk_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(vbk_vec, vbk_nll, pch=19, cex=2, xlab="Assumed vbk", ylab="NLL", ylim=c(0.8*min(vbk_nll), max(vbk_nll)*1.2))
points(vbk_toUse, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## M
out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(M_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_vec[i], lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
M_nll <- rep(NA, length(M_vec))
for(i in 1:length(M_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	M_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "M_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(M_vec, M_nll, pch=19, cex=2, xlab="Assumed M", ylab="NLL", ylim=c(0.8*min(M_nll), max(M_nll)*1.2))
points(M_toUse, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## ML50
out <- file.path(res_dir, "sens_ML50")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(ML50_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=ML50_vec[i], M95=NULL, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
ML50_nll <- rep(NA, length(ML50_vec))
for(i in 1:length(ML50_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	ML50_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "ML50_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(ML50_vec, ML50_nll, pch=19, cex=2, xlab="Assumed ML50", ylab="NLL", ylim=c(0.8*min(ML50_nll), max(ML50_nll)*1.2))
points(plist_m$ML50, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


#### output
## base model
out <- file.path(res_dir, "LCm")
dir.create(out, showWarnings=FALSE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))
			LBSPR <- readRDS(file.path(res_dir, "LBSPR_LCy", "LBSPR_results.rds"))

	F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
  	F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)
  	Report$F_y[length(Report$F_y)]
  	Report$F_y[length(Report$F_y)]/(F40*Inputs$Data$n_s)
  	Report$F_y[length(Report$F_y)]/(F30*Inputs$Data$n_s)

  	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	SPR_ci <- read_sdreport(sd, log=FALSE)
	SPR_ci[length(Report$SPR)]
	SPR_ci[length(Report$SPR)+1]
	Report$SPR[length(Report$SPR)]
	Report$sigma_R



# dev.new()
	all_years <- months_t
	lc_years <- months_o
	true_years <- years_t

png(file.path(fig_dir, "BaseModelOutput.png"), height=10, width=18, units="in", res=200)
plot_output(all_years=months_t, lc_years=months_o, true_years=years_t, Report=Report, Inputs=Inputs, Sdreport=Sdreport, LBSPR=LBSPR, lh=plist_m)
dev.off()

png(file.path(fig_dir, "BaseModelLCfits.png"), width=16, height=14, res=200, units="in")
plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=NULL, LBSPR=NULL, ylim=c(0,0.3))
dev.off()

## compare models
png(file.path(fig_dir, "CompareModels.png"), height=8, width=15, res=200, units="in")
dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "LCm_highN")
dir_list[[2]] <- file.path(res_dir, "LCm_baskettrap")
dir_list[[3]] <- file.path(res_dir, "LCm")

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}


col_list <- list()
col_list[[1]] <- col2hex("orange")
col_list[[2]] <- col2hex("orangered")
col_list[[3]] <- "#0000AA"

by <- round(length(true_years)/5)
lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
ilab <- which(true_years %in% lab)
ilab2 <- sapply(1:length(ilab), function(x){
  sub <- which(Inputs$Data$S_yrs %in% ilab[x])
  return(sub[length(sub)])
})

par(mfrow=c(1,3), mar=c(4,5,1,1), omi=c(0.2,0.2,0.2,0.2))
plot(x=1,y=1,type="n", xlim=c(1,length(true_years)), ylim=c(0,3.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Fishing mortality", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	if(Inputs$Data$n_s==1){
	  xY <- seq_along(all_years)
	  xLC <- which(all_years %in% lc_years)
	}
	if(Inputs$Data$n_s>1){
	  xY <- 1:Inputs$Data$n_y
	  xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
	}
	  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
	  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)	

	lines(x=xY, y=Report$F_y, lwd=2, col=col_list[[i]])
	points(x=xLC, y=Report$F_y[xLC], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
	  sd[,2][which(is.na(sd[,2]))] <- 0
	    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)  
	}
	if(i==1) axis(1, cex.axis=2, at=ilab, labels=lab)
	  if(i==length(dir_list)){
	  	abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2)
	  	abline(h=F30*Inputs$Data$n_s, lwd=2, lty=3)
	  }
if(i==1) legend("bottomleft", legend=c("All data", "Basket trap only", "High sampling months"), col=rev(unlist(col_list)), lty=1, lwd=4, cex=2, bty="n")
}

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,2.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Recruitment", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$R_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)	
}	

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,0.6), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="SPR", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
dev.off()

## compare selectivity models
png(file.path(fig_dir, "CompareSelex.png"), height=8, width=15, res=200, units="in")
dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "LCm_highdome")
dir_list[[2]] <- file.path(res_dir, "LCm_lowdome")
dir_list[[3]] <- file.path(res_dir, "LCm")

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}


col_list <- list()
col_list[[1]] <- col2hex("orange")
col_list[[2]] <- col2hex("orangered")
col_list[[3]] <- "#0000AA"

par(mfrow=c(2,2), mar=c(4,5,1,1), omi=c(0.2,0.2,0.2,0.2))
plot(x=1,y=1,type="n", xlim=c(1, length(plist_m$S_l)), ylim=c(0,1.1), xaxs="i", yaxs="i", cex.axis=2, xlab="Length bin (cm)", ylab="Selectivity", cex.lab=2, xaxt="n", yaxt="n")
axis(2, cex.axis=2, at=seq(0,1,by=0.25))
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=1:length(Report$S_l), y=Report$S_l, lwd=4, col=col_list[[i]])
	if(all(is.na(Sdreport))==FALSE){
	  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),]
	  sd[,2][which(is.na(sd[,2]))] <- 0
	    polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)  
	}
	if(i==1) axis(1, cex.axis=2, at=pretty(c(min(plist_m$highs),max(plist_m$highs))))

}

plot(x=1,y=1,type="n", xlim=c(1,length(true_years)), ylim=c(0,3.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Fishing mortality", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	if(Inputs$Data$n_s==1){
	  xY <- seq_along(all_years)
	  xLC <- which(all_years %in% lc_years)
	}
	if(Inputs$Data$n_s>1){
	  xY <- 1:Inputs$Data$n_y
	  xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
	}
	  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
	  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)	

	lines(x=xY, y=Report$F_y, lwd=2, col=col_list[[i]])
	points(x=xLC, y=Report$F_y[xLC], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
	  sd[,2][which(is.na(sd[,2]))] <- 0
	    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)  
	}
	if(i==1) axis(1, cex.axis=2, at=ilab, labels=lab)
	  if(i==length(dir_list)){
	  	abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2)
	  	abline(h=F30*Inputs$Data$n_s, lwd=2, lty=3)
	  }
if(i==1) legend("topright", legend=c("Logistic", "Low dome", "High dome"), col=rev(unlist(col_list)), lty=1, lwd=4, cex=2, bty="n")

}

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,2.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Recruitment", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$R_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)	
}	

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,0.6), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="SPR", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}

dev.off()


## sensitivity to biological inputs
png(file.path(fig_dir, "SensitivityBioInputs.png"), height=8, width=15, res=200, units="in")
dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_linf",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_linf",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}


col_list <- list()
col_list[[1]] <- col2hex("orange")
col_list[[2]] <- "#0000AA"
col_list[[3]] <- col2hex("orangered")

by <- round(length(true_years)/5)
lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
ilab <- which(true_years %in% lab)
ilab2 <- sapply(1:length(ilab), function(x){
  sub <- which(Inputs$Data$S_yrs %in% ilab[x])
  return(sub[length(sub)])
})

par(mfrow=c(2,2), mar=c(0,0,0,0), omi=c(1,1,1,1))
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	# if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
print.letter("Asymptotic length", xy=c(0.25,0.9), cex=2, font=2)
axis(2, seq(0,1,by=0.5), cex.axis=2)

dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_vbk",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_vbk",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	# if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
print.letter("von Bertalanffy k", xy=c(0.25,0.9), cex=2, font=2)

dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_M",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_M",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
print.letter("Natural mortality", xy=c(0.25,0.9), cex=2, font=2)
axis(2, seq(0,1,by=0.5), cex.axis=2)

dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_ML50",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_ML50",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
legend("topright", legend=c("-25%", "base", "+25%"), col=unlist(col_list), lwd=4, bty="n", cex=2, title="Input value")
print.letter("Length at 50% maturity", xy=c(0.25,0.9), cex=2, font=2)
mtext(side=1, "Year", cex=2, line=4, outer=TRUE)
mtext(side=2, "SPR", cex=2, line=4, outer=TRUE)
dev.off()