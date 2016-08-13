LIME_fits <- function(Inputs, Report, Sdreport, save){

        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1.96*InputMat[,2]), rev(exp(InputMat[,1]+1.96*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1.96*InputMat[,2], rev(InputMat[,1]+1.96*InputMat[,2])))
        } 

      Nyears <- Inputs$Data$n_t
  if(any(is.na(Sdreport))) next

  if(save==TRUE) png(file.path(fig_dir, "Time_series_estimates.png"), width=16, height=10, res=200, units="in")
  
  par(mfrow=c(2,2), mar=c(1,5,0,0), omi=c(0.7,0.5,0.5,0.5))
	Mat <- cbind("Year"=1:Nyears, "Est"=Report$F_t)
	ymax <- 5
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), xaxt="n", yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
    mtext(side=2, "Fishing mortality", line=3, cex=1.5)
  
	Mat <- cbind("Year"=1:Nyears, "Est"=Report$R_t_hat)
	ymax <- 3
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), xaxt="n", yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
    mtext(side=2, "Relative recruitment", line=3, cex=1.5)

	Mat <- cbind("Year"=1:Nyears, "Est"=Report$Depl)
	ymax <- 3
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)	
    mtext(side=2, "Relative biomass", line=3, cex=1.5)

  Mat <- cbind("Year"=1:Nyears, "Est"=Report$SPR_t)
  ymax <- 3
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
    mtext(side=2, "SPR", line=3, cex=1.5)

    mtext("Year", outer=TRUE, side=1, line=2, cex=1.5)


  if(save==TRUE) dev.off()
}