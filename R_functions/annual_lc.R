annual_lc <- function(year, data){
	sub <- data[which(data$Year==year),]
	sublc <- sub[,grepl("LC_", colnames(sub))]
	lc <- colSums(sublc)
	return(lc)
}
