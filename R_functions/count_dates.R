count_dates <- function(year, data){
	sub <- data[which(data$Year==year),]
	ndates <- length(unique(sub$Date))
	return(ndates)

}