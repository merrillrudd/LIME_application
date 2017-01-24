freq_summary <- function(year, bins, data, type, wt=FALSE){
	sub <- data[which(data$Year==year),]
	if(nrow(sub)==0) return(NA)
	if(type=="weight") freq <- sub[,grepl("W_", colnames(sub))]
	if(type=="length") freq <- sub[,grepl("LC_", colnames(sub))]
	total <- sum(colSums(freq)*bins)
	if(type=="length") return(total/sum(freq))
	if(type=="weight" & wt==FALSE) return(total)
	if(type=="weight" & wt==TRUE) return(total/length(unique(sub$Date)))
}