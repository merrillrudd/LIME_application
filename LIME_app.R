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
obsData <- formatDataList2(species="CRSNAP", data_dir=file.path(crsnap_dir, "data"), predata_yrs=10)
saveRDS(obsData, file.path(crsnap_dir, "data", "CRSNAP_data.rds"))

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
	if(id_sigma=="RecGrowthIndexVars") est_sigma <- c("log_sigma_R", "log_CV_L", "log_sigma_I")

	## run assessment model and save final gradients, parameter names, and estimates to check convergence

	dfout[[dd]] <- runModel(modpath=alldirs[dd], itervec=NULL, est_sigma=est_sigma, data_avail=data_avail, lh_list=cr_lh, rewrite=TRUE, start_f=0, simulation=FALSE, input_data=obsData)
}

end_run <- Sys.time() - start_run



## ---------------- explore results ------------------

aic <- calc_AIC(modpath_vec=alldirs[-c(5:6)])

## move forward with model 3 -- including abundance index and estimating growth variation

## ---------------- Figures ---------------------------

fig_dir <- file.path(crsnap_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

## length frequency
data_raw <- read.csv(file.path(crsnap_dir, "data", "cr_snapper_filtered.csv"))
png(file=file.path(fig_dir, "length_frequency_catch.png"), width=9, height=7, units="in", res=200)
lf <- length_frequency(binwidth=1, linf=cr_lh$linf, lmat=cr_lh$ML50, data=data_raw, plot=TRUE, weight=TRUE)
dev.off()

png(file=file.path(fig_dir, "cpue.png"), width=9, height=7, units="in", res=200)
par(mfrow=c(1,1))
pIt <- rep(NA, length(obsData$years))
pIt[which(obsData$years_i %in% names(obsData$I_t))] <- obsData$I_t
plot(x=obsData$years, y=pIt, xlim=c(obsData$years[1], max(obsData$years)), ylim=
c(0, max(pIt,na.rm=TRUE)*1.2), pch=19, cex=1.5, xlab="Year", ylab="Catch per unit effort", cex.lab=1.5, cex.axis=1.5, xaxs="i", yaxs="i")
dev.off()

## model fits
dd <- 3
choose_dir <- alldirs[dd]
	Inputs <- readRDS(file.path(choose_dir, "Inputs2.rds"))
	Report <- readRDS(file.path(choose_dir, "Report.rds"))
	Sdreport <- readRDS(file.path(choose_dir, "Sdreport.rds"))
	Quants <- readRDS(file.path(choose_dir, "Derived_quants.rds"))
	flag <- ifelse(file.exists(file.path(choose_dir, "NAs_final_gradient.txt"))|file.exists(file.path(choose_dir, "high_final_gradient.txt")), TRUE, FALSE)

	LIME_fits(Inputs, Report, Sdreport, obsData, save=TRUE)

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

### plot model fits - length composition
png(file.path(fig_dir, "LC_fits.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
LC_yrs <- Inputs$Data$LC_yrs
T_yrs <- Inputs$Data$T_yrs
for(t in 1:length(LC_yrs)){
	plot(obsData$LFprop[t,], pch=19, lwd=2, ylim=c(0, max(obsData$LFprop)))
	lines(Report$plb[LC_yrs[t],], lwd=4, col="red")
}
dev.off()

## kobe plot
png(file.path(fig_dir, "base_kobe.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(1,1))
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

## remove index
png(file.path(fig_dir, "rmIndex.png"), width=8, height=7, units="in", res=200)
par(mfrow=c(1,1))
plot(x=1, y=1, xlim=c(1,max(datyr)), ylim=c(0,10), type="n")
for(dd in 3:4){
	rep <- readRDS(file.path(alldirs[dd], "Report.rds"))
	sdrep <- readRDS(file.path(alldirs[dd], "Sdreport.rds"))
	Mat <- cbind("Year"=1:length(rep$Depl), "Est"=rep$Depl)
	if(dd==3){
		lines(x=Mat[,"Year"], y=Mat[,"Est"], col="red", lwd=4)
	   if(all(is.na(sdrep))==FALSE)if( !("condition" %in% names(attributes(sdrep)))) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.3), border=NA)
	}
	if(dd==4){
	   lines(x=Mat[,"Year"], y=Mat[,"Est"], col="red", lwd=4)
	   if(all(is.na(sdrep))==FALSE)if( !("condition" %in% names(attributes(sdrep)))) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
	} 
}
dev.off()

### retrospective

retro_dir <- file.path(crsnap_dir, "retrospective")
dir.create(retro_dir, showWarnings=FALSE)

datyr <- Inputs$Data$LC_yrs
datyr_i <- 1:Inputs$Data$n_lc
for(i in 1:length(datyr_i)){

	dir_new <- file.path(retro_dir, i)
	dir.create(dir_new, showWarnings=FALSE)

	dat_new <- obsData
	dat_new$years <- obsData$years[1:datyr[i]]
	dat_new$years_i <- obsData$years_i[1:datyr[i]]
	if(i==1){
		dat_new$LF <- t(as.matrix(obsData$LF[1:datyr_i[i],], nrow=i, ncol=ncol(obsData$LF)))
		rownames(dat_new$LF) <- datyr[i]
		dat_new$LFprop <- t(as.matrix(obsData$LFprop[1:datyr_i[i],], nrow=i, ncol=ncol(obsData$LFprop)))
		rownames(dat_new$LFprop) <- datyr[i]
	}
	if(i>1){
		dat_new$LF <- obsData$LF[1:datyr_i[i],]
		dat_new$LFprop <- obsData$LF[1:datyr_i[i],]
	}
	dat_new$I_t <- obsData$I_t[1:datyr_i[i]]
	dat_new$ML_t <- obsData$ML_t[1:datyr_i[i]]
	dat_new$Nyears <- length(dat_new$years)
	dat_new$Nyears_comp <- nrow(dat_new$LFprop)
	dat_new$obs_per_year <- obsData$obs_per_year[1:datyr[i]]

	data_avail <- "Index_LC"
	est_sigma <- c("log_sigma_R", "log_CV_L")

	## run assessment model and save final gradients, parameter names, and estimates to check convergence

	retro <- runModel(modpath=dir_new, itervec=NULL, est_sigma=est_sigma, data_avail=data_avail, lh_list=cr_lh, rewrite=FALSE, start_f=0, simulation=FALSE, input_data=dat_new)
}

png(file.path(fig_dir, "retro_F.png"), width=8, height=7, units="in", res=200)
par(mfrow=c(1,1))
plot(x=1, y=1, xlim=c(1,max(datyr)), ylim=c(0,10), type="n")
for(i in 1:length(datyr_i)){
	rinp <- readRDS(file.path(retro_dir, i, "Inputs2.rds"))
	rrep <- readRDS(file.path(retro_dir, i, "Report.rds"))
	rsdrep <- readRDS(file.path(retro_dir, i, "Sdreport.rds"))
	Mat <- cbind("Year"=1:rinp$Data$n_t, "Est"=rrep$F_t)
	if(i==length(datyr)){
		lines(x=Mat[,"Year"], y=Mat[,"Est"], col="red", lwd=4)
	   if(all(is.na(rsdrep))==FALSE)if( !("condition" %in% names(attributes(rsdrep)))) polygon( y=FUN(summary(rsdrep)[which(rownames(summary(rsdrep))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.3), border=NA)
	}
	if(i != length(datyr)){
	   lines(x=Mat[,"Year"], y=Mat[,"Est"], col="red", lwd=4)
	   if(all(is.na(rsdrep))==FALSE)if( !("condition" %in% names(attributes(rsdrep)))) polygon( y=FUN(summary(rsdrep)[which(rownames(summary(rsdrep))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
	} 
}
dev.off()

png(file.path(fig_dir, "retro_SPR.png"), width=8, height=7, units="in", res=200)
par(mfrow=c(1,1))
plot(x=1, y=1, xlim=c(1,max(datyr)), ylim=c(0,1), type="n")
for(i in 1:length(datyr_i)){
	rinp <- readRDS(file.path(retro_dir, i, "Inputs2.rds"))
	rrep <- readRDS(file.path(retro_dir, i, "Report.rds"))
	rsdrep <- readRDS(file.path(retro_dir, i, "Sdreport.rds"))
	Mat <- cbind("Year"=1:rinp$Data$n_t, "Est"=rrep$SPR_t)
	if(i==length(datyr)){
		lines(x=Mat[,"Year"], y=Mat[,"Est"], col="red", lwd=4)
	   if(all(is.na(rsdrep))==FALSE)if( !("condition" %in% names(attributes(rsdrep)))) polygon( y=FUN(summary(rsdrep)[which(rownames(summary(rsdrep))=="SPR_t"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.3), border=NA)
	}
	if(i != length(datyr)){
	   lines(x=Mat[,"Year"], y=Mat[,"Est"], col="red", lwd=4)
	   if(all(is.na(rsdrep))==FALSE)if( !("condition" %in% names(attributes(rsdrep)))) polygon( y=FUN(summary(rsdrep)[which(rownames(summary(rsdrep))=="SPR_t"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
	} 
}
dev.off()


### re-do plots from assessment 
### sensitivities
### other formal model comparison/selection
### potentially adjust LIME to only estimate 1 parameter selectivity
