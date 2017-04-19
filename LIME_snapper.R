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
## Species directories
###################################

sp_dir <- file.path(main_dir, "costa_rican_snapper")
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

## add life history information for species chosen for assessment
## siganus sutor
plist <- readRDS(file.path(data_dir, "CRSNAP_life_history_annual.rds"))
plist_q <- readRDS(file.path(data_dir, "CRSNAP_life_history_quarterly.rds"))
plist_m <- readRDS(file.path(data_dir, "CRSNAP_life_history_monthly.rds"))

###################################
## Load data
###################################

## read nominal CPUE Index
I_t <- readRDS(file.path(data_dir, "Nominal_CPUE.rds"))

## read weighted data
LCprop <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LFprop_weighted.rds"))
LC <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_weighted.rds"))
bins <- as.numeric(colnames(LC))

LCq <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_weighted_quarterly.rds"))

## weighted data by gear
LCprop_bygear <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LFprop_gear_year.rds"))
gearmat <- readRDS(file.path(data_dir, "Samples_per_gear.rds"))
gears <- colnames(gearmat)

LCup <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_weighted_Nup.rds"))
LCq_up <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LFquarterly_weighted_Nup.rds"))

## LC separate by gears
LC_bline <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_g1.rds"))
LCq_bline <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_g1_quarterly.rds"))

#########################################
## run assessments -- life history fixed
#########################################

years_o <- as.numeric(rownames(LC))
years_t <- min(years_o):max(years_o)

quarters_o <- as.numeric(rownames(LCq))
quarters_t <- 1:max(quarters_o)

input_data_y <- list("years"=years_t, "LF"=LC)
input_data_q <- list("years"=quarters_t, "LF"=LCq)
input_data_bline <- list("years"=years_t, "LF"=LC_bline)
input_data_up <- list("years"=years_t, "LF"=LCup)
input_data_bline_q <- list("years"=quarters_t, "LF"=LCq_bline)
input_data_up_q <- list("years"=quarters_t, "LF"=LCq_up)


############### weighted length comps ##########################
## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_highN")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_up, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_highN <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_bottomlongline")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_bline, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_bl <- readRDS(file.path(out, "LBSPR_results.rds"))

## quarterly LBPSR
out <- file.path(res_dir, "LBSPR_LCq")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_q, rewrite=TRUE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

## quarterly LBPSR
out <- file.path(res_dir, "LBSPR_LCq_highN")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_up_q, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_highN <- readRDS(file.path(out, "LBSPR_results.rds"))

## quarterly LBPSR
out <- file.path(res_dir, "LBSPR_LCq_bottomlongline")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_bline_q, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_bl <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LIME
out <- file.path(res_dir, "LCy")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_y$years, lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, true_years=years_t, lh=plist)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, ylim=c(0,0.3))
			dev.off()


### weighted length comp - counts only - all years - gear1
out <- file.path(res_dir, "LCy_bottomlongline")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist, input_data=input_data_bline, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_bline$years, lc_years=rownames(input_data_bline$LF), LBSPR=lbspr_res_bl, lh=plist, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_bline$LF), LBSPR=lbspr_res_bl, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCy_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist, input_data=input_data_up, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=lbspr_res_highN, lh=plist, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=lbspr_res_highN, ylim=c(0,0.3))
			dev.off()

###########################
## likelihood profiles
###########################

linf_vec <- seq(0.75*plist$linf, 1.25*plist$linf, length.out=10)
vbk_vec <- seq(0.75*plist$vbk, 1.25*plist$vbk, length.out=10)
ML50_vec <- seq(0.75*plist$ML50, 1.25*plist$ML50, length.out=10)

M_tool <- c(0.29, 0.23, 0.25, 0.25, 0.13, 0.16, 0.32, 0.34, 0.48, 0.41,0.43)
M_vec <- seq(min(M_tool), max(M_tool), length.out=10)

## linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(linf_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","linf"), val_adjust=c(0.2,0.7,linf_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCy", "Report.rds"))
linf_like <- rep(NA, length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	linf_like[i] <- rep$jnll
}
png(file.path(figs_dir, "linf_likeprof.png"), height=10, width=14, res=200, units="in")
plot(linf_vec, linf_like, pch=19, cex=2, xlab="Assumed Linf", ylab="NLL")
points(plist$linf, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(vbk_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","vbk"), val_adjust=c(0.2,0.7,vbk_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCy", "Report.rds"))
vbk_like <- rep(NA, length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	vbk_like[i] <- rep$jnll
}
png(file.path(figs_dir, "vbk_likeprof.png"), height=10, width=14, res=200, units="in")
plot(vbk_vec, vbk_like, pch=19, cex=2, xlab="Assumed vbk", ylab="NLL")
points(plist$vbk, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## M
out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(M_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","M"), val_adjust=c(0.2,0.7,M_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCy", "Report.rds"))
M_like <- rep(NA, length(M_vec))
for(i in 1:length(M_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	M_like[i] <- rep$jnll
}
png(file.path(figs_dir, "M_likeprof.png"), height=10, width=14, res=200, units="in")
plot(M_vec, M_like, pch=19, cex=2, xlab="Assumed M", ylab="NLL", )
points(plist$M, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


####################################
### MONTE CARLO
####################################
n <- 1000
dist <- data.frame("linf"=rlnorm(n,log(plist$linf),0.25), "vbk"=rlnorm(n,log(plist$vbk),0.25), "M"=rlnorm(n,log(plist$M),0.25))


par(mfrow=c(1,3), omi=c(0.5,1,0.5,0.5))
hist(dist$linf, col="goldenrod", xlab="", ylab="", main="", cex.axis=2, cex.lab=2, border=NA, yaxs="i")
mtext(side=3, "Asymptotic length (Linf)", cex=2)
abline(v=plist$linf, lwd=2, col="blue")
box()

hist(dist$vbk, col="goldenrod", xlab="", ylab="", main="", cex.axis=2, cex.lab=2, border=NA, yaxs="i")
mtext(side=3, "von Bertalanffy (k)", cex=2)
abline(v=plist$vbk, lwd=2, col="blue")
box()

hist(dist$M, col="goldenrod", xlab="", ylab="", main="", cex.axis=2, cex.lab=2, border=NA, yaxs="i")
mtext(side=3, "Natural mortality (M)", cex=2)
abline(v=plist$M, lwd=2, col="blue")
box()
mtext("Frequency", side=2, outer=TRUE, cex=2)