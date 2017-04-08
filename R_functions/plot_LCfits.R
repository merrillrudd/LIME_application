plot_LCfits <- function(Inputs, Inputs2=NULL, Inputs3=NULL, Inputs4=NULL, Report=NULL, Report2=NULL, LBSPR=NULL, true_lc_years=NULL, ylim=NULL){
	# dev.new()

	obs <- Inputs$LF
	lbhighs <- colnames(obs)
	lc_years <- rownames(obs)

	if(all(is.null(Inputs2))==FALSE){
		obs2 <- Inputs2$LF
		lbhighs2 <- colnames(obs2)
		lc_years2 <- rownames(obs2)
	}
	if(all(is.null(Inputs3))==FALSE){
		obs3 <- Inputs3$LF
		lbhighs3 <- colnames(obs3)
		lc_years3 <- rownames(obs3)
	}
	if(all(is.null(Inputs4))==FALSE){
		obs4 <- Inputs4$LF
		lbhighs4 <- colnames(obs4)
		lc_years4 <- rownames(obs4)
	}


	dim <- c(ceiling(sqrt(length(lc_years))), ceiling(sqrt(length(lc_years))))


	if(all(is.null(Report))==FALSE){
		pred <- Report$plb
		Tyrs <- Inputs$T_yrs
	}
	if(all(is.null(Report))){
		pred <- NULL
		Tyrs <- NULL
	}
	if(all(is.null(LBSPR))==FALSE) pred2 <- t(LBSPR$pLF)
	if(all(is.null(LBSPR))) pred2 <- NULL

	par(mfrow=dim, mar=c(0,0,0,0), omi=c(1,1,1,1))

	# find_max <- 0
	# for(i in 1:length(lc_years)){
	# 	yr <- lc_years[i]
	# 		omax <- max(obs[which(lc_years==yr),]/sum(obs[which(lc_years==yr),]))
	# 		if(all(is.null(Inputs2))==FALSE) omax2 <- max(obs2[which(lc_years2==yr),]/sum(obs2[which(lc_years2==yr),]))
	# 		if(all(is.null(Inputs2))) omax2 <- NULL
	# 		if(all(is.null(Report))==FALSE) pmax <- max(pred[which(Tyrs==yr),])
	# 		if(all(is.null(Report))) pmax <- NULL
	# 		max <- max(c(omax,pmax,omax2))
	# 		if(max > find_max) find_max <- max
	# }
	if(all(is.null(ylim))) ylim <- c(0, 0.1)
	for(i in 1:length(lc_years)){
		yr <- lc_years[i]
		barplot(obs[which(lc_years==yr),]/sum(obs[which(lc_years==yr),]), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim, col="#22222250", border=NA, space=0)
		lines(pred[which(Tyrs==yr),], col="blue", lwd=2)
		lines(pred2[which(lc_years==yr),], col="red", lwd=2)
		box()
		if(all(is.null(Inputs2))==FALSE){
			par(new=TRUE)
			barplot(obs2[which(lc_years2==yr),]/sum(obs2[which(lc_years2==yr),]), border=NA, space=0, col="#DD000050", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim)
		}
		if(all(is.null(Inputs3))==FALSE){
			par(new=TRUE)
			barplot(obs3[which(lc_years2==yr),]/sum(obs3[which(lc_years2==yr),]), border=NA, space=0, col="#0000DD50", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim)
		}
		if(all(is.null(Inputs4))==FALSE){
			par(new=TRUE)
			barplot(obs4[which(lc_years2==yr),]/sum(obs4[which(lc_years2==yr),]), border=NA, space=0, col="#00DD0050", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim)
		}

		if(i %in% (length(lc_years)-dim[2]+1):length(lc_years)) axis(1, at=pretty(seq_along(lbhighs)), labels=pretty(as.numeric(lbhighs)), cex.axis=2)
		if(i %in% seq(1,length(lc_years), by=dim[2])) axis(2, at=pretty(ylim), las=2, cex.axis=2)
		if(all(is.null(true_lc_years))==FALSE & length(true_lc_years)!=length(lc_years)){
			text(x=0.2*max(seq_along(lbhighs)), y=0.9*ylim[2], lc_years[i], font=2, cex=2)
			warning("Input years for length composition data do not match number of years in analysis")
		}
		if(all(is.null(true_lc_years))) text(x=0.2*max(seq_along(lbhighs)), y=0.9*ylim[2], lc_years[i], font=2, cex=2)
		if(all(is.null(true_lc_years))==FALSE) text(x=0.2*max(seq_along(lbhighs)), y=0.9*ylim[2], true_lc_years[i], font=2, cex=2)
		if(all(is.null(Report))==FALSE) abline(v=Report$S50, lty=2)
	}
	mtext(side=1, "Length bin (cm)", outer=TRUE, line=4, cex=1.5)
	mtext(side=2, "Proportion", outer=TRUE, line=5, cex=1.5)

}