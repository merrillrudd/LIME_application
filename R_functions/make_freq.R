make_freq <- function(mvec, bins){
	mvec1 <- as.numeric(mvec[which(is.na(mvec)==FALSE)])
	mvec2 <- mvec1[which(is.na(mvec1)==FALSE)]
	outvec <- rep(0, length(bins))
	for(i in 1:length(mvec2)){
		outvec[which(bins>=mvec2[i])[1]] <- outvec[which(bins>=mvec2[i])[1]] + 1
	}
	return(outvec)
}
