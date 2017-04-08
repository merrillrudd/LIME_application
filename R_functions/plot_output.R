plot_output <- function(all_years, lc_years, Inputs, Report, Sdreport, LBSPR=NULL, True=NULL, lh){

# dev.new()
par(mfrow=c(2,3), mar=c(4,5,2,2))
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)

if(Inputs$Data$n_s==1){
  xF <- seq_along(all_years)
  xLC <- which(all_years %in% lc_years)
}
if(Inputs$Data$n_s>1){
  xF <- 1:Inputs$Data$n_y
  xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
}
ylim <- c(0, max(1.1, max(Report$F_y)*1.1))
plot(x=xF, y=Report$F_y, lwd=2, col="blue", ylim=ylim, type="l", xaxt="n", ylab="Fishing mortality", xlab="Time", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=xLC, y=Report$F_y[xLC], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
  sd[,2][which(is.na(sd[,2]))] <- 0
    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)  
}
  axis(1, cex.axis=2, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$F_t, col="black", lwd=2)
  abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2, col=gray(0.4))
  abline(h=F30*Inputs$Data$n_s, lwd=2, lty=2)

  plot(x=seq_along(all_years), y=Report$R_t, lwd=2, col="blue", ylim=c(0, max(Report$R_t)*1.5), type="l", xaxt="n", ylab="Recruitment", xlab="Time", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
}
  axis(1, cex.axis=2, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$R_t, col="black", lwd=2)

  plot(x=seq_along(all_years), y=Report$L_t_hat, lwd=2, col="blue", ylim=c(0, max(Report$L_t_hat)*1.5), type="l", xaxt="n", ylab="Mean length", xlab="Time", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$L_t_hat[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
}
  axis(1, cex.axis=2, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$ML_t, col="black", lwd=2)

  plot(x=seq_along(all_years), y=Report$SB_t, lwd=2, col="blue", ylim=c(0, max(Report$SB_t)*1.5), type="l", xaxt="n", ylab="Relative spawning biomass", xlab="Time", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$SB_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
}
  axis(1, cex.axis=2, at=pretty(seq_along(all_years)), labels=pretty(all_years))
  if(all(is.null(True))==FALSE) lines(True$SB_t, col="black", lwd=2)


  plot(x=xF, y=Report$SPR_t, lwd=2, col="blue", ylim=c(0, 1), type="l", xaxt="n", ylab="SPR", xlab="Time", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=xLC, y=Report$SPR_t[xLC], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
} 
  axis(1, cex.axis=2, at=pretty(xF), labels=pretty(xF))
  if(all(is.null(True))==FALSE) lines(True$SPR_t, col="black", lwd=2)
  if(all(is.null(LBSPR))==FALSE){
    lines(x=xLC, y=LBSPR$SPR, col="red", lwd=2)
    points(x=xLC, y=LBSPR$SPR, col="red", pch=19)
  } 
  abline(h=0.4, col=gray(0.4), lwd=2, lty=2)
  abline(h=0.3, lwd=2, lty=2)


  mids <- seq(Inputs$Data$binwidth, by=Inputs$Data$binwidth, length.out=ncol(Inputs$Data$LF))
  plot(x=1:length(mids), y=Report$S_l, lwd=2, col="blue", ylim=c(0, 1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
if(all(is.na(Sdreport))==FALSE)  polygon( y=read_sdreport(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),], log=FALSE), x=c(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE), rev(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, cex.axis=2, at=pretty(seq_along(Inputs$Data$lbhighs)), labels=pretty(Inputs$Data$lbhighs))
  if(all(is.null(True))==FALSE) lines(True$S_l, col="black", lwd=2)
  if(all(is.null(LBSPR))==FALSE){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      lines(x=1:length(mids), y=S_l2, col="#AA000050", lwd=2)
    }
  }
  legend("bottomright", col=c("blue", "red"), lwd=2, legend=c("LIME", "LB-SPR"), cex=1.7)


}