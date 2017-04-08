plot_LCdata <- function(lc_list, lc_years, gears, L50=NULL, S50=NULL, ylim=NULL){

	ramp <- colorRamp(c("darkorange", "purple"))
	if(all(is.null(gears))) ncol <- 1
	if(all(is.null(gears))==FALSE) ncol=length(gears)
	col_vec <- rgb( ramp(seq(0, 1, length = ncol)), max = 255)
	names(col_vec) <- gears
	dim <- c(ceiling(sqrt(length(lc_years))), ceiling(sqrt(length(lc_years))))
	par(mfrow=c(dim), mar=c(0,0,0,0), omi=c(1,1,1,1))
	if(all(is.null(ylim))) ylim <- c(0,0.1)
	for(y in 1:length(lc_years)){
		if(length(lc_list)==1){
			lc <- lc_list[[1]]
			bins <- as.numeric(colnames(lc))
			sub <- lc[which(rownames(lc)==lc_years[y]),]/sum(lc[which(rownames(lc)==lc_years[y]),])
			barplot(sub, col=col_vec[1], lwd=4, xaxs="i", yaxs="i", xlim=c(0,length(bins)), ylim=ylim, xaxt="n", yaxt="n", border=NA, space=0)
			text(x=0.2*length(bins), y=0.9*max(ylim), lc_years[y], font=2, cex=2)
			box()
			if(is.null(L50)==FALSE) abline(v=which(bins>=ML50)[1], col="black", lwd=2, lty=2)
			if(is.null(S50)==FALSE) abline(v=which(bins>=S50)[1], col="gray", lwd=2, lty=2)
			labs <- seq(bins[1], max(bins), by=bins[1]*2)
			if(y %in% (length(lc_years)-dim[1]+1):length(lc_years)) axis(1, at=which(bins %in% labs), labels=labs, cex.axis=2)
			if(y %% dim[1]==1) axis(2, pretty(ylim), cex.axis=2, las=2)
		}
		if(length(lc_list)>1){
			barplot(lc_list[[y]], col=col_vec[which(names(col_vec) %in% rownames(lc_list[[y]]))], border=NA, lwd=4, xaxs="i", yaxs="i", xlim=c(0,length(bins)), ylim=c(0, 0.1), xaxt="n", yaxt="n", space=0)
			text(x=0.2*length(bins), y=0.9*max(ylim), lc_years[y], font=2, cex=2)
			box()
			if(is.null(L50)==FALSE) abline(v=which(bins>=ML50)[1], col="black", lwd=2, lty=2)
			if(is.null(S50)==FALSE) abline(v=which(bins>=S50)[1], col="gray", lwd=2, lty=2)
			labs <- seq(bins[1], max(bins), by=bins[1]*2)
			if(y %in% (length(lc_years)-dim[1]+1):length(lc_years)) axis(1, at=which(bins %in% labs), labels=labs, cex.axis=2)
			if(y %% dim[1]==1) axis(2, pretty(ylim), cex.axis=2, las=2)
				gears_name <- gears
				gears_name[which(gears=="")] <- "Unclassified"
			if(y==length(lc_years)) legend("topright", legend=gears_name, col=col_vec, pch=19, cex=1.4)
		}
	}
	mtext("Length bin (cm)", side=1, outer=TRUE, line=4, cex=1.5)
	mtext("Proportion", side=2, outer=TRUE, line=5, cex=1.5)



}