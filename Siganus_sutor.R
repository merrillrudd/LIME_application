rm(list=ls())

## load LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

## load LBSPR package
devtools::install_github("adrianhordyk/LBSPR", build.vignettes=TRUE, dependencies=TRUE)
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
stock_dir <- file.path(region_dir, "Siganus_sutor")
dir.create(stock_dir, showWarnings=FALSE)

data_dir <- file.path(region_dir, "data")
dir.create(data_dir, showWarnings=FALSE)

figs_dir <- file.path(stock_dir, "figures")
dir.create(figs_dir, showWarnings=FALSE)
######################################
## Species life history
######################################

plist <- readRDS(file.path(data_dir, "Siganus_sutor_life_history_annual.rds"))
plist_m <- readRDS(file.path(data_dir, "Siganus_sutor_life_history_monthly.rds"))

######################################
## Load data
######################################

LCraw <- readRDS(file.path(data_dir, "Siganus_sutor_LCraw.rds"))
LCm_raw <- readRDS(file.path(data_dir, "Siganus_sutor_LCraw_monthly.rds"))

LCweight <- readRDS(file.path(data_dir, "Siganus_sutor_LC_weighted.rds"))
LCm_weight <- readRDS(file.path(data_dir, "Siganus_sutor_LCmonthly_weighted.rds"))

LCm_bt <- readRDS(file.path(data_dir, "Siganus_sutor_LCmonthly_baskettrap.rds"))
LCm_highN <- readRDS(file.path(data_dir, "Siganus_sutor_LCmonthly_weighted_Nup.rds"))

LCy_bt <- readRDS(file.path(data_dir, "Siganus_sutor_LC_baskettrap.rds"))
LCy_up <- readRDS(file.path(data_dir, "Siganus_sutor_LC_weighted_Nup.rds"))
####################################
## run assessments
####################################

res_dir <- file.path(stock_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LCraw))
years_t <- min(years_o):max(years_o)

months_o <- as.numeric(rownames(LCm_raw))
months_t <- min(months_o):max(months_o)

input_data_y <- list("years"=years_t, "LF"=LCweight)
input_data_m <- list("years"=months_t, "LF"=LCm_weight)
input_data_bt <- list("years"=months_t, "LF"=LCm_bt)
input_data_up <- list("years"=months_t, "LF"=LCm_highN)
input_data_y_bt <- list("years"=years_t, "LF"=LCy_bt)
input_data_y_up <- list("years"=years_t, "LF"=LCy_up)

############### weighted length comps ##########################
## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))
  	lbspr_res$SPR[length(lbspr_res$SPR)] - 1.96*sqrt(lbspr_res$SPR_Var[length(lbspr_res$SPR_Var)])
  	lbspr_res$SPR[length(lbspr_res$SPR)] + 1.96*sqrt(lbspr_res$SPR_Var[length(lbspr_res$SPR_Var)])



	png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
	plot_LCfits(Inputs=input_data_y, Report=NULL, true_lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, ylim=c(0,0.2), ML50=plist$ML50, SL50=lbspr_res$SL50)
	dev.off()

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_baskettrap")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y_bt, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_bt <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_highN")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y_up, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_highN <- readRDS(file.path(out, "LBSPR_results.rds"))


## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_altM")
dir.create(out, showWarnings=FALSE)

plist_new <- plist
plist_new$M <- 1.49
	run <- run_LBSPR(modpath=out, lh=plist_new, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_altM <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LIME
out <- file.path(res_dir, "LCy")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

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

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			base_rep <- Report

			 F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
			 Report$F_t/F40
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=TRUE)[c(length(Report$F_y), length(Report$F_y)+1)]
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=FALSE)[c(length(Report$SPR_t), length(Report$SPR_t)+1)]
    		pthresh <- pnorm(0.3, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		ptarget <- pnorm(0.45, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]) - pnorm(0.35, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		pup <- 1-pnorm(0.40, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])


			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=lbspr_res, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()



### weighted length comp - counts only - all months - basket trap only
out <- file.path(res_dir, "LCm_baskettrap")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_bt, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_bt$years, lc_years=rownames(input_data_bt$LF), LBSPR=lbspr_res_bt, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_bt$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCm_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_up, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=lbspr_res_highN, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - all months
out <- file.path(res_dir, "LCm_altM")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI", "M"), val_adjust=c(0.2,0.737,0.2,0.2, 1.49), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))


			 F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
			 Report$F_t/F40
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=TRUE)[c(length(Report$F_y), length(Report$F_y)+1)]
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=FALSE)[c(length(Report$SPR_t), length(Report$SPR_t)+1)]
    		pthresh <- pnorm(0.3, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		ptarget <- pnorm(0.45, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]) - pnorm(0.35, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])



			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=lbspr_res_altM, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


### weighted length comp - counts only - all months
out <- file.path(res_dir, "LCm_sigRstartlow")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.3,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

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


### weighted length comp - counts only - all months
out <- file.path(res_dir, "LCm_selexstartalt")
dir.create(out, showWarnings=FALSE)

	plist_m_alt <- create_lh_list(vbk=plist$vbk, linf=plist$linf, lwa=plist$lwa, lwb=plist$lwb, S50=20, S95=25, selex_input="length", M50=plist$ML50, maturity_input="length", SigmaR=0.7, SigmaF=0.2, M=plist$M, F1=NULL, CVlen=0.1, nseasons=12, binwidth=1)

				res <- run_LIME(modpath=out, lh=plist_m_alt, input_data=input_data_m, est_sigma=c("log_sigma_R"), data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

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


## dome-shaped selex
plist_m_lowdome <- create_lh_list(vbk=plist$vbk, linf=plist$linf, t0=plist$t0, M=plist$M, lwa=plist$lwa, lwb=plist$lwb, S50=plist$SL50, S95=plist$SL95, M50=plist$ML50, M95=plist$ML95, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=41)
plist_m_highdome <- create_lh_list(vbk=plist$vbk, linf=plist$linf, t0=plist$t0, M=plist$M, lwa=plist$lwa, lwb=plist$lwb, S50=plist$SL50, S95=plist$SL95, M50=plist$ML50, M95=plist$ML95, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=22)


par(mfrow=c(1,1), mar=c(5,5,4,4))
plot(plist_m$S_l, type="o", lwd=2, pch=19, cex.axis=2, xlab="Length (cm)", ylab="Proportion", cex.lab=2)
lines(plist_m_lowdome$S_l, type="o", lwd=2, lty=2, col="steelblue", pch=19)
lines(plist_m_highdome$S_l, type="o", lwd=2, lty=3, col="tomato", pch=19)

## low-dome
out <- file.path(res_dir, "LCm_lowdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=plist_m_lowdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=lbspr_res, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


## high-dome
out <- file.path(res_dir, "LCm_highdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=plist_m_highdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

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



###########################
## likelihood profiles
###########################

## sigma R
out <- file.path(res_dir, "sens_SigmaR")
dir.create(out, showWarnings=FALSE)

sigR_vec <- seq(0.05, 0.95, length=10)

for(i in 1:length(sigR_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_m, input_data=input_data_m, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR"), val_adjust=c(0.2,sigR_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
sigR_like <- rep(NA, length(sigR_vec))
col_conv <- rep("black", length(sigR_vec))
for(i in 1:length(sigR_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	sigR_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "sigR_likeprof.png"), height=10, width=14, res=200, units="in")
plot(sigR_vec, sigR_like, pch=19, cex=2, xlab="Fixed SigmaR", ylab="NLL", col=col_conv)
points(rep1$sigma_R, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


## check fishbase and natural mortality toolfor life history distributions
proper_name <- "Siganus sutor"
genus <- strsplit(proper_name, " ")[[1]][1]

# growth2 <- popgrowth(species_list(Genus=genus))
# saveRDS(growth2, file.path(data_dir, "genus_info.rds"))

## linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

linf_genus <- growth2$Loo
linf_vec <- seq(min(linf_genus), max(linf_genus), length=10)

for(i in 1:length(linf_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_m, input_data=input_data_m, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI", "linf"), val_adjust=c(0.2,0.737,0.2,0.2, linf_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
linf_like <- rep(NA, length(linf_vec))
col_conv <- rep("black", length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	linf_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "linf_likeprof.png"), height=10, width=14, res=200, units="in")
plot(linf_vec, linf_like, pch=19, cex=2, xlab="Assumed Linf", ylab="NLL", col=col_conv)
points(plist_m$linf, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

vbk_genus <- growth2$K
vbk_vec <- seq(min(vbk_genus), max(vbk_genus), length=10)

for(i in 1:length(vbk_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_m, input_data=input_data_m, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","vbk"), val_adjust=c(0.2,0.737,vbk_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
vbk_like <- rep(NA, length(vbk_vec))
col_conv <- rep("black", length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	vbk_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "vbk_likeprof.png"), height=10, width=14, res=200, units="in")
plot(vbk_vec, vbk_like, pch=19, cex=2, xlab="Assumed vbk", ylab="NLL", col=col_conv)
points(plist_m$vbk, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## M
out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

M_tool <- c(1.37, 1.28, 1.37, 1.35, 0.95, 0.53, 1.30, 1.39, 1.88, 1.65)
M_vec <- seq(min(M_tool), max(M_tool), length=10)

for(i in 1:length(M_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

			res <- run_LIME(modpath=out2, lh=plist_m, input_data=input_data_m, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI", "M"), val_adjust=c(0.2,0.737,0.2,0.2, M_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
M_like <- rep(NA, length(M_vec))
col_conv <- rep("black", length(M_vec))
for(i in 1:length(M_vec)){
	inp <- readRDS(file.path(out, i, "Inputs.rds"))
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	M_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "M_likeprof.png"), height=10, width=14, res=200, units="in")
plot(M_vec, M_like, pch=19, cex=2, xlab="Assumed M", ylab="NLL", col=col_conv)
points(plist_m$M*12, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


png(file.path(figs_dir, "like_prof_compare.png"), height=10, width=15, units="in", res=200)

mat <- matrix(c(1,2,3,4,
				1,2,3,4,
				5,6,7,8,
				5,6,7,8,
				9,9,9,9,
				10,11,12,13,
				10,11,12,13), nrow=7, ncol=4, byrow=TRUE)
nf <- layout(mat)
layout.show(nf)

par(mar=c(0,4,0,0), omi=c(0.5,0.5,0.5,0.5))
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

## using input AgeMax=4, Linf=36.2, vbk=0.87, Age at maturity(years)=1
M_tool <- c(1.37, 1.28, 1.37, 1.35, 0.95, 0.53, 1.30, 1.39, 1.88, 1.65)
M_toUse <- median(M_tool)
M_h <- hist(M_tool, plot=FALSE)
M_x <- seq(min(M_tool), max(M_tool), length=60)
M_y <- dnorm(M_x, mean=mean(M_tool), sd=sd(M_tool))
M_y2 <- (M_y*diff(M_h$mids[1:2])*length(M_tool))/sum((M_y*diff(M_h$mids[1:2])*length(M_tool)))
plot(M_x, M_y2, col="black", lwd=2, type="l", xaxt="n", xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
rug(M_tool, ticksize=0.15, lwd=2)
# M_d <- rtnorm(1, mean=mean(M_tool), sd=sd(M_tool), lower=min(M_tool), upper=max(M_tool))
# abline(v=M_d, lty=2, col="red", lwd=2)
abline(v=plist$M, lwd=2, col="blue")
mtext(side=3, "M", cex=2.5, line=1)
mtext("Probability", side=2, line=3, cex=2)
print.letter("(a)", xy=c(0.9,0.95), cex=2)

linf_genus <- growth2$Loo
linf_h <- hist(linf_genus, plot=FALSE)
linf_x <- seq(min(linf_genus), max(linf_genus), length=60)
linf_y <- dnorm(linf_x, mean=mean(linf_genus), sd=sd(linf_genus))
linf_y2 <- linf_y*diff(linf_h$mids[1:2])*length(linf_genus)/sum(linf_y*diff(linf_h$mids[1:2])*length(linf_genus))
plot(linf_x, linf_y2, col="black", lwd=2, type="l", xaxt="n",xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
rug(linf_genus, ticksize=0.15, lwd=2)
# linf_d <- rtnorm(1, mean=mean(linf_genus), sd=sd(linf_genus), lower=min(linf_genus), upper=max(linf_genus))
# abline(v=linf_d, lty=2, col="red", lwd=2)
abline(v=plist$linf, lwd=2, col="blue")
mtext(side=3, "Linf", cex=2.5, line=1)
print.letter("(b)", xy=c(0.9,0.95), cex=2)

vbk_genus <- growth2$K
vbk_h <- hist(vbk_genus, plot=FALSE)
vbk_x <- seq(min(vbk_genus), max(vbk_genus), length=60)
vbk_y <- dnorm(vbk_x, mean=mean(vbk_genus), sd=sd(vbk_genus))
vbk_y2 <- vbk_y*diff(vbk_h$mids[1:2])*length(vbk_genus)/sum(vbk_y*diff(vbk_h$mids[1:2])*length(vbk_genus))
plot(vbk_x, vbk_y2, col="black", lwd=2, type="l", xaxt="n",xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
rug(vbk_genus, ticksize=0.15, lwd=2)
# vbk_d <- rtnorm(1, mean=mean(vbk_genus), sd=sd(vbk_genus), lower=min(vbk_genus), upper=max(vbk_genus))
# abline(v=vbk_d, lty=2, col="red", lwd=2)
abline(v=plist$vbk, lwd=2, col="blue")
mtext(side=3, "k", cex=2.5, line=1)
print.letter("(c)", xy=c(0.9,0.95), cex=2)

sigR_x <- seq(0, 1.5, length=60)
sigR_y <- dnorm(sigR_x, mean=mean(0.737), sd=0.353)/sum(dnorm(sigR_x, mean=mean(0.737), sd=0.353))
plot(sigR_x, sigR_y, col="black", lwd=2, type="l", xaxt="n",xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
# sigR_d <- rtnorm(1, mean=mean(sigR_genus), sd=sd(sigR_genus), lower=min(sigR_genus), upper=max(sigR_genus))
# abline(v=sigR_d, lty=2, col="red", lwd=2)
abline(v=0.737, lwd=2, col="blue")
mtext(side=3, "SigmaR", cex=2.5, line=1)
print.letter("(d)", xy=c(0.9,0.95), cex=2)



out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

M_vec <- seq(min(M_tool), max(M_tool), length=10)

M_like <- rep(NA, length(M_vec))
col_conv <- rep("black", length(M_vec))
for(i in 1:length(M_vec)){
	inp <- readRDS(file.path(out, i, "Inputs.rds"))
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	M_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
i1 <- 1
i2 <- 8
col_conv[i1] <- "orange"
col_conv[i2] <- "orangered"
plot(x=M_vec, y=M_like, pch=19, cex=3, col=col_conv, xpd=NA, cex.axis=2, xlab="", ylab="", yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(M_like),max(M_like)))[-length(pretty(c(min(M_like),max(M_like))))])
points(plist_m$M*12, rep1$jnll, pch=19, cex=3, col="blue")
mtext(side=2, "NLL", cex=2, line=3)
print.letter("(e)", xy=c(0.9,0.95), cex=2)

### linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

linf_genus <- growth2$Loo
linf_vec <- seq(min(linf_genus), max(linf_genus), length=10)

linf_like <- rep(NA, length(linf_vec))
col_conv_linf <- rep("black", length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	linf_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv_linf[i] <- "gray"
}
i1 <- which(col_conv_linf=="black")[1]
i2 <- 8
col_conv_linf[i1] <- "orange"
col_conv_linf[i2] <- "orangered"
plot(x=linf_vec, y=linf_like, pch=19, cex=3, col=col_conv_linf, xpd=NA, cex.axis=2, xlab="", ylab="", yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(linf_like),max(linf_like)))[-length(pretty(c(min(linf_like),max(linf_like))))])
points(plist_m$linf, rep1$jnll, pch=19, cex=3, col="blue")
print.letter("(f)", xy=c(0.9,0.95), cex=2)

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

vbk_genus <- growth2$K
vbk_vec <- seq(min(vbk_genus), max(vbk_genus), length=10)

vbk_like <- rep(NA, length(vbk_vec))
col_conv_vbk <- rep("black", length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	vbk_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv_vbk[i] <- "gray"
}
i1 <- 2
i2 <- 5
col_conv_vbk[i1] <- "orange"
col_conv_vbk[i2] <- "orangered"
plot(x=vbk_vec, y=vbk_like, pch=19, cex=3, col=col_conv_vbk, xpd=NA, cex.axis=2, xlab="", ylab="", yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(vbk_like),max(vbk_like)))[-length(pretty(c(min(vbk_like),max(vbk_like))))])
points(plist_m$vbk, rep1$jnll, pch=19, cex=3, col="blue")
print.letter("(g)", xy=c(0.9,0.95), cex=2)


## sigma R
out <- file.path(res_dir, "sens_SigmaR")
dir.create(out, showWarnings=FALSE)

sigR_vec <- seq(0.05, 0.95, length=10)

sigR_like <- rep(NA, length(sigR_vec))
col_conv_sigR <- rep("black", length(sigR_vec))
for(i in 1:length(sigR_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	sigR_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv_sigR[i] <- "gray"
}
i1 <- 4
i2 <- 7
col_conv_sigR[i1] <- "orange"
col_conv_sigR[i2] <- "orangered"
plot(x=sigR_vec, y=sigR_like, pch=19, cex=3, col=col_conv_sigR, xpd=NA, cex.axis=2, xlab="", ylab="", xlim=c(min(sigR_x), max(sigR_x)), yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(sigR_like),max(sigR_like)))[-length(pretty(c(min(sigR_like),max(sigR_like))))])
points(rep1$sigma_R, rep1$jnll, pch=19, cex=3, col="blue")
print.letter("(h)", xy=c(0.9,0.95), cex=2)


plot(x=1,y=1,type="n", axes=F, ann=F)
### SPR
out <- file.path(res_dir, "sens_M")
i1 <- 1
i2 <- 8
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=2, "SPR", line=3, cex=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(i)", xy=c(0.9,0.95), cex=2)

out <- file.path(res_dir, "sens_linf")
i1 <- which(col_conv_linf=="black")[1]
i2 <- 8
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(j)", xy=c(0.9,0.95), cex=2)

out <- file.path(res_dir, "sens_vbk")
i1 <- 2
i2 <- 5
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(k)", xy=c(0.9,0.95), cex=2)

out <- file.path(res_dir, "sens_SigmaR")
i1 <- 4
i2 <- 7
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(l)", xy=c(0.9,0.95), cex=2)

dev.off()



### model fits figure
out <- file.path(res_dir, "LCm")
dir.create(out, showWarnings=FALSE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			base_rep <- Report

out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

	obs1 <- input_data_m$LF
	all_mos <- 1:(length(years_t)*12)
	mo_yrs <- as.vector(sapply(1:length(years_t), function(x) rep(years_t[x], 12)))
	choose_mo <- rep(NA, length(years_o))
	for(i in 1:length(years_o)){
		sub_mo <- months_o[which(months_o %in% all_mos[which(mo_yrs==years_o[i])])]
		choose_mo[i] <- sub_mo[1]
	}

	obs2 <- input_data_y$LF

png(file.path(figs_dir, "CompareModelFits.png"), height=15, width=8, res=200, units="in")
par(mfcol=c(length(years_o),2), mar=c(0,0,0,0), omi=c(1,1,1,1))
for(i in 1:length(years_o)){
	barplot(obs1[which(rownames(obs1)==choose_mo[i]),]/sum(obs1[which(rownames(obs1)==choose_mo[i]),]), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,0.2), col=gray(0.6), border=NA, space=0)
	lines(base_rep$plb[choose_mo[i],], lwd=4, col="blue")
	if(i==length(years_o)) axis(1, cex.axis=2)
	axis(2, cex.axis=2, at=c(0,0.15), las=2)
	abline(v=plist$ML50, lty=2)
	box()
}
for(i in 1:length(years_o)){
	barplot(obs2[which(rownames(obs2)==years_o[i]),]/sum(obs2[which(rownames(obs2)==years_o[i]),]), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,0.2), col=gray(0.6), border=NA, space=0)
	lines(lbspr_res$pLF[,i], lwd=4, col="red")
	if(i==length(years_o)) axis(1, cex.axis=2)
	print.letter(years_o[i], xy=c(0.8,0.8), cex=3)
	abline(v=plist$ML50, lty=2)
	box()
}
mtext(side=1, "Length bin (cm)", cex=2, line=4, outer=TRUE)
mtext(side=2, "Proportions", cex=2, line=5, outer=TRUE)
dev.off()




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