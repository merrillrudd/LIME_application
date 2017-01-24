runSensitivity <- function(param, val, modpath, lfdata, sp_lfdata, species, obs_mod, data_avail, id_sigma, lh, rewrite=FALSE){

	lh_new <- lh
	lh_new[[param]] <- val

	ignore <- runAssessment(modpath=modpath, lfdata=lfdata, sp_lfdata=sp_lfdata, species=species, obs_mod=obs_mod, data_avail=data_avail, id_sigma=id_sigma, lh=lh_new, rewrite=rewrite)
	
}