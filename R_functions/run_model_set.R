run_model_set <- function(dir, LC, lh, add_years=0, S_l_input=-1, C_t=NULL, I_t=NULL, LBSPR=TRUE, LIME=TRUE, F_up=5, est_sigma="log_sigma_R", fix_param=FALSE){
	
	years_o <- as.numeric(rownames(LC))
	years_m <- (min(years_o)-add_years):max(years_o)

	input_data <- list("years"=years_m, "LF"=LC, "C_t"=C_t, I_t=I_t)

	if(LBSPR==TRUE){
		run <- run_LBSPR(modpath=dir, lh=lh, species="x", input_data=input_data, rewrite=TRUE, simulation=FALSE, write=TRUE)	
		lbspr_res <- readRDS(file.path(dir, "LBSPR_results.rds"))
	}

	if(LIME==TRUE){
		res <- run_LIME(modpath=dir, write=TRUE, lh=plist, input_data=input_data, est_sigma=est_sigma, data_avail="LC", LFdist=1, simulation=FALSE, rewrite=TRUE, S_l_input=S_l_input, F_up=F_up, fix_param=fix_param)

		if(LBSPR==FALSE & file.exists(file.path(dir, "LBSPR_results.rds"))) lbspr_res <- readRDS(file.path(dir, "LBSPR_results.rds"))
		if(LBSPR==FALSE & file.exists(file.path(dir, "LBSPR_results.rds"))==FALSE) lbspr_res <- NULL
				
		Report <- readRDS(file.path(dir, "Report.rds"))
		Sdreport <- readRDS(file.path(dir, "Sdreport.rds"))
		Inputs <- readRDS(file.path(dir, "Inputs.rds"))

		png(file.path(dir, "output.png"), width=16, height=10, res=200, units="in")
		plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
		dev.off()				

		png(file.path(dir, "LCfits.png"), width=16, height=10, res=200, units="in")
		plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
		dev.off()
	}

}