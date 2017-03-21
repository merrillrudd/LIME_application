plot_output <- function(all_years, lc_years, Inputs, Report, Sdreport, LBSPR=NULL, True=NULL){

# dev.new()
par(mfrow=c(2,3))
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=Inputs$Data$n_a), Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=seq(0,by=1,length.out=Inputs$Data$n_a), Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)

ylim <- c(0, max(1.1, max(Report$F_t)*1.1))
plot(x=seq_along(all_years), y=Report$F_t, lwd=2, col="blue", ylim=ylim, type="l", xaxt="n", ylab="Fishing mortality", xlab="Year", xaxs="i", yaxs="i")
points(x=which(all_years %in% lc_years), y=Report$F_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$F_t, col="black", lwd=2)
  abline(h=F40, lwd=2, lty=2, col=gray(0.4))
  abline(h=F30, lwd=2, lty=2)

  plot(x=seq_along(all_years), y=Report$R_t, lwd=2, col="blue", ylim=c(0, max(Report$R_t)*1.5), type="l", xaxt="n", ylab="Recruitment", xlab="Year", xaxs="i", yaxs="i")
points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$R_t, col="black", lwd=2)

  plot(x=seq_along(all_years), y=Report$L_t_hat, lwd=2, col="blue", ylim=c(0, max(Report$L_t_hat)*1.5), type="l", xaxt="n", ylab="Mean length", xlab="Year", xaxs="i", yaxs="i")
points(x=which(all_years %in% lc_years), y=Report$L_t_hat[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$ML_t, col="black", lwd=2)

  plot(x=seq_along(all_years), y=Report$SB_t, lwd=2, col="blue", ylim=c(0, max(Report$SB_t)*1.5), type="l", xaxt="n", ylab="Spawning biomass", xlab="Year", xaxs="i", yaxs="i")
points(x=which(all_years %in% lc_years), y=Report$SB_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$SB_t, col="black", lwd=2)


  plot(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col="blue", ylim=c(0, 1), type="l", xaxt="n", ylab="SPR", xlab="Year", xaxs="i", yaxs="i")
points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$SPR_t, col="black", lwd=2)
  if(all(is.null(LBSPR))==FALSE){
    lines(x=which(all_years %in% lc_years), y=LBSPR$SPR, col="red", lwd=2)
    points(x=which(all_years %in% lc_years), y=LBSPR$SPR, col="red", pch=19)
  } 
  abline(h=0.4, col=gray(0.4), lwd=2, lty=2)
  abline(h=0.3, lwd=2, lty=2)


  plot(x=1:length(Report$S_l), y=Report$S_l, lwd=2, col="blue", ylim=c(0, 1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i")
  polygon( y=read_sdreport(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),], log=FALSE), x=c(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE), rev(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, at=pretty(seq_along(Inputs$Data$lbhighs)), labels=pretty(Inputs$Data$lbhighs))
  if(all(is.null(True))==FALSE) lines(True$S_l, col="black", lwd=2)
  if(all(is.null(LBSPR))==FALSE){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      LMids <- seq(Inputs$Data$binwidth, by=Inputs$Data$binwidth, length.out=nrow(LBSPR$pLF))
      S_l2 <- 1.0/(1+exp(-log(19)*(LMids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      lines(x=LMids, y=S_l2, col="#AA000050", lwd=2)
    }
  }


}