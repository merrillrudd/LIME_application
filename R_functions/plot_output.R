plot_output <- function(all_years, lc_years, Inputs, Report, Sdreport, LBSPR=NULL, lh, true_years){

if(Inputs$Data$n_s==1){
  xY <- seq_along(all_years)
  xLC <- which(all_years %in% lc_years)
}
if(Inputs$Data$n_s>1){
  xY <- 1:Inputs$Data$n_y
  xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
}
par(mfrow=c(2,3), mar=c(4,5,2,2))
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)

ylim <- c(0, max(Report$F_y)*1.5)
plot(x=xY, y=Report$F_y, lwd=2, col="blue", ylim=ylim, type="l", xaxt="n", xaxs="i", yaxs="i", cex.axis=2, ylab="Fishing mortality", xlab="Year", cex.lab=2)
points(x=xLC, y=Report$F_y[xLC], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
  sd[,2][which(is.na(sd[,2]))] <- 0
    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)  
}
by <- round(length(true_years)/5)
lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
ilab <- which(true_years %in% lab)
  axis(1, cex.axis=2, at=ilab, labels=lab)
  abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2)
  abline(h=F30*Inputs$Data$n_s, lwd=2, lty=3)

plot(x=seq_along(all_years), y=Report$R_t, lwd=2, col="blue", ylim=c(0, max(Report$R_t)*1.5), type="l", xaxt="n", ylab="Recruitment", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
}
ilab2 <- sapply(1:length(ilab), function(x){
  sub <- which(Inputs$Data$S_yrs %in% ilab[x])
  return(sub[length(sub)])
})
axis(1, cex.axis=2, at=ilab2, labels=lab)

plot(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col="blue", ylim=c(0, 1), type="l", xaxt="n", ylab="SPR", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
} 
axis(1, cex.axis=2, at=ilab2, labels=lab)
if(all(is.null(LBSPR))==FALSE){
  par(new=TRUE)
  plot(LBSPR$SPR, xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0, 1), xaxt="n", yaxt="n", lwd=2, col="red", type="l")
  points(x=xLC, LBSPR$SPR[xLC], col="red", pch=19)
}
  abline(h=0.4, lwd=2, lty=2)
  abline(h=0.3, lwd=2, lty=3)


plot(x=seq_along(all_years), y=Report$L_t_hat, lwd=2, col="blue", ylim=c(0, max(Report$L_t_hat)*1.5), type="l", xaxt="n", ylab="Mean length", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$L_t_hat[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
}
axis(1, cex.axis=2, at=ilab2, labels=lab)

plot(x=seq_along(all_years), y=Report$D_t, lwd=2, col="blue", ylim=c(0, max(Report$D_t)*1.5), type="l", xaxt="n", ylab="Relative spawning biomass", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
points(x=which(all_years %in% lc_years), y=Report$D_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA)
if(all(is.na(Sdreport))==FALSE){
  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]
sd[,2][which(is.na(sd[,2]))] <- 0
  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
}
axis(1, cex.axis=2, at=ilab2, labels=lab)
  
  mids <- seq(Inputs$Data$binwidth, by=Inputs$Data$binwidth, length.out=ncol(Inputs$Data$LF))
  plot(x=1:length(mids), y=Report$S_l, lwd=2, col="blue", ylim=c(0, 1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
if(all(is.na(Sdreport))==FALSE)  polygon( y=read_sdreport(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),], log=FALSE), x=c(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE), rev(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  axis(1, cex.axis=2, at=pretty(seq_along(Inputs$Data$lbhighs)), labels=pretty(Inputs$Data$lbhighs))
  if(all(is.null(LBSPR))==FALSE){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      lines(x=1:length(mids), y=S_l2, col="#AA000050", lwd=2)
    }
  }
  legend("bottomright", col=c("blue", "red", "black", "black"), lwd=2, legend=c("LIME", "LB-SPR", "SPR 40%", "SPR 30%"), cex=1.7, lty=c(1,1,2,3))


}