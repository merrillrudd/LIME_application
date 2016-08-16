rm(list=ls())

###################################
## Packages
###################################

devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dep=TRUE)
library(TMBhelper)

main_dir <- "C:\\Git_Projects\\LIME_application"

###################################
## Costa Rican spotted snapper
###################################

crsnap_dir <- file.path(main_dir, "CRSNAP")
dir.create(crsnap_dir, showWarnings=FALSE)

setwd(crsnap_dir)
source("R_functions\\functions.R")

## life history list
cr_lh <- choose_lh_list(species="CRSNAP", selex="asymptotic")

## ---------------- Model settings and directories ------------------

## data availability scenarios (can also set up for "Catch_LC" when both catch and length comp are available -- just need to make sure the data type is in the name -- with Index, Catch, and LC being the options)
avail_set <- c("Index_LC", "LC")

## estimate variances -- always estimating Recruitment variation (log_sigma_R), but could add on other variance parameters (match variance names exactly ** update manual) -- in this case could estimate the CV for the growth curve 
estsigma_set <- c("RecVar", "RecGrowthVars")

## setup combinations of models to run
modcombos <- as.matrix(expand.grid("Data_avail"=avail_set, "Est_variance"=estsigma_set))

crsnap_res_dir <- file.path(crsnap_dir, "results")

## remove results (for testing purposes)
# unlink(crsnap_res_dir, TRUE)

dir.create(crsnap_res_dir, showWarnings=FALSE)

## transform model combinations into directory names
alldirs <- model_paths(modcombos=modcombos, res_dir=crsnap_res_dir)

## ---------------- Format data ------------------

## below code formats data specifically for my examples, and saves to the folder for easy reading
## functions not included for sharing but contact me if you have questions of how to get data formatted
# obsData <- formatDataList(species="CRSNAP", data_dir=file.path(crsnap_dir, "data"))
# saveRDS(obsData, file.path(crsnap_dir, "data", "CRSNAP_data.rds"))

## get your data in the same form as obsData object
obsData <- readRDS(file.path(crsnap_dir, "data", "CRSNAP_data.rds"))


## ---------------- Assessment model ------------------

start_run <- Sys.time()

dfout <- list()
for(dd in 1:length(alldirs)){

	## get available data types from model path name, used for formatting TMB input
	data_avail <- ifelse(grepl("Index", alldirs[dd]), avail_set[which(grepl("Index", avail_set))], avail_set[which(grepl("Index", avail_set)==FALSE)])

	## get variance parameters to estimate from model name, used for formatting TMB input
	id_sigma <- estsigma_set[which(sapply(1:length(estsigma_set), function(x) grepl(estsigma_set[x], alldirs[dd])))]
	if(id_sigma=="RecVar") est_sigma <- "log_sigma_R"
	if(id_sigma=="RecGrowthVars") est_sigma <- c("log_sigma_R", "log_CV_L")

	## run assessment model and save final gradients, parameter names, and estimates to check convergence

	dfout[[dd]] <- runModel(modpath=alldirs[dd], itervec=NULL, est_sigma=est_sigma, data_avail=data_avail, lh_list=cr_lh, rewrite=TRUE, start_f=0, simulation=FALSE, input_data=obsData)
}

end_run <- Sys.time() - start_run

aic <- calc_AIC(modpath_vec=alldirs)

## move forward with model 3 -- including abundance index and estimating growth variation

## ---------------- explore results ------------------


## ---------------- Figures ------------------

fig_dir <- file.path(crsnap_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

## length frequency
data_raw <- read.csv(file.path(crsnap_dir, "data", "cr_snapper_filtered.csv"))
png(file=file.path(fig_dir, "length_frequency_catch.png"), width=9, height=7, units="in", res=200)
lf <- length_frequency(binwidth=1, linf=cr_lh$linf, lmat=cr_lh$ML50, data=data_raw, plot=TRUE, weight=TRUE)
dev.off()


## model fits
dd <- 4
choose_dir <- alldirs[dd]
	Inputs <- readRDS(file.path(choose_dir, "Inputs2.rds"))
	Report <- readRDS(file.path(choose_dir, "Report.rds"))
	Sdreport <- readRDS(file.path(choose_dir, "Sdreport.rds"))
	Quants <- readRDS(file.path(choose_dir, "Derived_quants.rds"))
	flag <- ifelse(file.exists(file.path(choose_dir, "NAs_final_gradient.txt"))|file.exists(file.path(choose_dir, "high_final_gradient.txt")), TRUE, FALSE)

	LIME_fits(Inputs, Report, Sdreport, save=TRUE)

## fishing mortality with reference points
png(file.path(fig_dir, "base_F.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(1,1))
Mat <- cbind("Year"=obsData$years, "Est"=Report$F_t)
ymax <- 5
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="", ylab="", xaxs="i", yaxs="i")
if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
abline(h=Quants$F30, col="blue", lwd=2)
abline(h=Quants$F40, col="blue", lwd=2, lty=2)
mtext("Year", side=1, cex=1.2, line=2.5)
mtext("Estimated Fishing Mortality", side=2, cex=1.2, line=2.5)
dev.off()

## SPR
png(file.path(fig_dir, "base_SPR.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(1,1))
Mat <- cbind("Year"=obsData$years, "Est"=Report$SPR_t)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="", ylab="", xaxs="i", yaxs="i")
if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)
mtext("Year", side=1, cex=1.2, line=2.5)
mtext("Estimated Spawning Potential Ratio", side=2, cex=1.2, line=2.5)
dev.off()

## kobe plot
png(file.path(fig_dir, "base_kobe.png"), width=8, height=7, res=200, units="in")
plot(Report$SPR_t, Report$F_t/Quants$F30, xlim=c(0,0.6), ylim=c(0, 1.1), pch=19, col="blue", xaxs="i", yaxs="i", xpd=NA, xlab="SPR", ylab=expression("F/F"[30]), type="o")
points(Report$SPR_t[length(obsData$years)], Report$F_t[length(obsData$years)]/Quants$F30, cex=2, col="blue")
# text(Report$SPR_t[length(obsData$years)], (Report$F_t[length(obsData$years)]/Quants$F30)-0.04, "2015")
# text(Report$SPR_t[length(obsData$years)-2], (Report$F_t[length(obsData$years)-2])/Quants$F30, "2013")
# text(Report$SPR_t[length(obsData$years)-4], (Report$F_t[length(obsData$years)-4])/Quants$F30, "2011")
abline(h=1, lty=2)
abline(v=0.3, lty=2)
polygon(x=c(0,1.1,1.1,0), y=c(1,1,1.1,1.1), border=NA, 
	col="#AA000050")
polygon(x=c(0,0.3,0.3,0), y=c(0,0,1.1,1.1), border=NA,
	col="#AA000050")
dev.off()

### plot model fits - length composition
### re-do plots from assessment 
### sensitivities
### other formal model comparison/selection
