runAssessment <- function(modpath, lfdata, sp_lfdata, species, obs_mod, data_avail, id_sigma, lh, rewrite=FALSE, sensitivity_inputs=NULL){

	###################################################################
	## ---------------- Format data ------------------
	###################################################################
	years_o <- unique(lfdata$Year)[order(unique(lfdata$Year))]
	years_t <- min(years_o):max(years_o)	

	data_list <- compile_data(data=sp_lfdata, species=species, model_years=years_t, obs=obs_mod)
	
	## get variance parameters to estimate from model name, used for formatting TMB input
	if(id_sigma=="RecVar") est_sigma <- "log_sigma_R"
	if(id_sigma=="RecGrowthVars") est_sigma <- c("log_sigma_R", "log_CV_L")
	if(id_sigma=="RecGrowthIndexVars") est_sigma <- c("log_sigma_R", "log_CV_L", "log_sigma_I")

	## run assessment model and save final gradients, parameter names, and estimates to check convergence to directory
	ignore <- runModel(modpath=modpath, itervec=NULL, est_sigma=est_sigma, data_avail=data_avail, lh_list=lh, rewrite=rewrite, start_f=0, simulation=FALSE, input_data=data_list, sensitivity_inputs=sensitivity_inputs, biascorrect=TRUE)
}