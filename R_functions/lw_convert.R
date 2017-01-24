## convert to weight
lw_convert <- function(species, lengths, lbins, wbins, freq=TRUE){

	lw <- tryCatch(length_weight(species, fields=c("a","b")), error=function(e) NA)
	if(all(is.na(lw))) rep(NA, length(wbins))
	if(all(is.na(lw))==FALSE){
		lwa <- mean(lw$a)
		lwb <- mean(lw$b)	

		length_obs <- unlist(sapply(1:length(lengths), function(x) rep(lbins[x],lengths[x])))
		weight <- round(lwa*as.numeric(length_obs)^lwb)
		wfreq <- make_freq(mvec=weight, bins=wbins)
		if(freq==TRUE) return(wfreq)
		if(freq==FALSE) return(weight)		
	}

}