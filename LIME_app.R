rm(list=ls())

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


## setup directory for figures
fig_dir <- file.path(crsnap_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

dd <- 4
choose_dir <- alldirs[dd]
	Inputs <- readRDS(file.path(choose_dir, "Inputs2.rds"))
	Report <- readRDS(file.path(choose_dir, "Report.rds"))
	Sdreport <- readRDS(file.path(choose_dir, "Sdreport.rds"))
	flag <- ifelse(file.exists(file.path(choose_dir, "NAs_final_gradient.txt"))|file.exists(file.path(choose_dir, "high_final_gradient.txt")), TRUE, FALSE)

	LIME_fits(Inputs, Report, Sdreport, save=FALSE)


