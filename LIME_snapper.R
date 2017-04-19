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
input_data_bline_q <- list("years"=quarters_t, "LF"=LCq_bline[which(rowSums(LCq_bline)>0),])
input_data_up_q <- list("years"=quarters_t, "LF"=LCq_up)


############### weighted length comps ##########################
## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=FALSE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_highN")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_up, rewrite=FALSE, simulation=FALSE)	
	lbspr_res_highN <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_bottomlongline")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_bline, rewrite=FALSE, simulation=FALSE)	
	lbspr_res_bl <- readRDS(file.path(out, "LBSPR_results.rds"))

## quarterly LBPSR
out <- file.path(res_dir, "LBSPR_LCq")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_q, rewrite=FALSE, simulation=FALSE)	
	lbspr_res_q <- readRDS(file.path(out, "LBSPR_results.rds"))

## quarterly LBPSR
out <- file.path(res_dir, "LBSPR_LCq_highN")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_up_q, rewrite=FALSE, simulation=FALSE)	
	lbspr_res_highN_q <- readRDS(file.path(out, "LBSPR_results.rds"))

## quarterly LBPSR
out <- file.path(res_dir, "LBSPR_LCq_bottomlongline")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_bline_q, rewrite=FALSE, simulation=FALSE)	
	lbspr_res_bl_q <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LIME
out <- file.path(res_dir, "LCy")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

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

## quarterly LIME
out <- file.path(res_dir, "LCq")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_q$years, lc_years=rownames(input_data_q$LF), LBSPR=lbspr_res, true_years=years_t, lh=plist_q)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_q$LF), LBSPR=lbspr_res_q, ylim=c(0,0.3))
			dev.off()


### weighted length comp - counts only - all years - gear1
out <- file.path(res_dir, "LCq_bottomlongline")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_bline_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_bline_q$years, lc_years=rownames(input_data_bline_q$LF), LBSPR=lbspr_res_bl_q, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_bline_q$LF), LBSPR=lbspr_res_bl_q, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCy_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist, input_data=input_data_up, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

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

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCq_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_up_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up_q$years, lc_years=rownames(input_data_up_q$LF), LBSPR=lbspr_res_highN_q, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up_q$LF), LBSPR=lbspr_res_highN_q, ylim=c(0,0.3))
			dev.off()

###########################
## likelihood profiles
###########################

## initial values (from Bystrom thesis)
linf_toUse <- plist$linf
vbk_toUse <- plist$vbk
t0_toUse <- plist$t0
lwa_toUse <- plist$lwa
lwb_toUse <- plist$lwb
M_toUse <- plist$M
ML50_toUse <- plist$ML50

## check fishbase and natural mortality toolfor life history distributions
proper_name <- "Lutjanus guttatus"
genus <- strsplit(proper_name, " ")[[1]][1]

# growth2 <- popgrowth(species_list(Genus=genus))
# saveRDS(growth2, file.path(data_dir, "genus_info.rds"))
growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

linf_genus <- growth2$Loo
linf_vec <- seq(min(linf_genus), max(linf_genus), length=10)

vbk_genus <- growth2$K
vbk_vec <- seq(min(vbk_genus), max(vbk_genus), length=10)

M_tool <- c(0.29, 0.23, 0.25, 0.25, 0.13, 0.16, 0.32, 0.34, 0.48, 0.41)
M_vec <- seq(min(M_tool), max(M_tool), length=10)

## linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(linf_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","linf"), val_adjust=c(0.2,0.2,linf_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
linf_like <- rep(NA, length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	linf_like[i] <- rep$jnll
}
png(file.path(figs_dir, "linf_likeprof.png"), height=10, width=14, res=200, units="in")
plot(linf_vec, linf_like, pch=19, cex=2, xlab="Assumed Linf", ylab="NLL")
points(plist_q$linf, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(vbk_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","vbk"), val_adjust=c(0.2,0.2,vbk_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
vbk_like <- rep(NA, length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	vbk_like[i] <- rep$jnll
}
png(file.path(figs_dir, "vbk_likeprof.png"), height=10, width=14, res=200, units="in")
plot(vbk_vec, vbk_like, pch=19, cex=2, xlab="Assumed vbk", ylab="NLL")
points(plist_q$vbk, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## M
out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(M_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","M"), val_adjust=c(0.2,0.2,M_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
M_like <- rep(NA, length(M_vec))
for(i in 1:length(M_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	M_like[i] <- rep$jnll
}
png(file.path(figs_dir, "M_likeprof.png"), height=10, width=14, res=200, units="in")
plot(M_vec, M_like, pch=19, cex=2, xlab="Assumed M", ylab="NLL", )
points(plist_q$M, rep1$jnll, pch=19, cex=2, col="blue")
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


############################
### life history figures

linf_genus <- growth2$Loo
linf_h <- hist(linf_genus, plot=FALSE)
linf_x <- seq(min(linf_genus), max(linf_genus), length=60)
linf_y <- dnorm(linf_x, mean=mean(linf_genus), sd=sd(linf_genus))
linf_y2 <- linf_y*diff(linf_h$mids[1:2])*length(linf_genus)
plot(linf_x, linf_y2, col="blue", lwd=2, type="l")
rug(linf_genus)
linf_d <- rtnorm(1, mean=mean(linf_genus), sd=sd(linf_genus), lower=min(linf_genus), upper=max(linf_genus))
abline(v=linf_d, lty=2, col="red", lwd=2)
abline(v=linf_toUse, lwd=2)

vbk_genus <- growth2$K
vbk_h <- hist(vbk_genus, plot=FALSE)
vbk_x <- seq(min(vbk_genus), max(vbk_genus), length=60)
vbk_y <- dnorm(vbk_x, mean=mean(vbk_genus), sd=sd(vbk_genus))
vbk_y2 <- vbk_y*diff(vbk_h$mids[1:2])*length(vbk_genus)
plot(vbk_x, vbk_y2, col="blue", lwd=2, type="l")
rug(vbk_genus)
vbk_d <- rtnorm(1, mean=mean(vbk_genus), sd=sd(vbk_genus), lower=min(vbk_genus), upper=max(vbk_genus))
abline(v=vbk_d, lty=2, col="red", lwd=2)
abline(v=vbk_toUse, lwd=2)


## using input AgeMax=22, Linf=64, vbk=0.21, Age at maturity(years)=4, asymptotic weight = 2750g.
M_tool <- c(0.29, 0.23, 0.25, 0.25, 0.13, 0.16, 0.32, 0.34, 0.48, 0.41)
M_toUse <- median(M_tool)
M_h <- hist(M_tool, plot=FALSE)
M_x <- seq(min(M_tool), max(M_tool), length=60)
M_y <- dnorm(M_x, mean=mean(M_tool), sd=sd(M_tool))
M_y2 <- M_y*diff(M_h$mids[1:2])*length(M_tool)
plot(M_x, M_y2, col="blue", lwd=2, type="l")
rug(M_tool)
M_d <- rtnorm(1, mean=mean(M_tool), sd=sd(M_tool), lower=min(M_tool), upper=max(M_tool))
abline(v=M_d, lty=2, col="red", lwd=2)
abline(v=M_toUse, lwd=2)
