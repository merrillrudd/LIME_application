compile_data <- function(data, species, model_years, obs="high", index=NULL, index_years=NULL, catch=NULL, catch_years=NULL, meanlen=NULL, meanlen_years=NULL){
	
	sp_data <- data[[species]]

	years_list <- unique(sp_data$Year)[order(unique(sp_data$Year))]
	years_i <- which(model_years %in% years_list)

	lcomp_list <-  t(sapply(1:length(years_list), function(x) annual_lc(year=years_list[x], data=sp_data)))
	rownames(lcomp_list) <- years_i

	if(obs=="high") obsperyr_list <- rowSums(lcomp_list)
	if(obs=="low") obsperyr_list <- sapply(1:length(years_list), function(x) count_dates(years_list[x], data=sp_data))
	if(obs=="vhigh") obsperyr_list <- rep(1000, length(years_list))
	names(obsperyr_list) <- years_i

	data_list <- formatData(lfreq=lcomp_list, lfreq_years=years_list, obs_per_year=obsperyr_list, lbins=lbins, model_years=model_years, index=index, index_years=index_years, catch=catch, catch_years=catch_years, meanlen=meanlen, meanlen_years=meanlen_years)
	return(data_list)
}