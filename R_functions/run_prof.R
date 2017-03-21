run_prof <- function(vec, param, dir, lh, input_data, rewrite=TRUE){

  for(i in 1:length(vec)){
	out <- file.path(dir, paste0(param, "_", vec[i]))
	dir.create(out, showWarnings=FALSE)

	param_adjust <- param
	val_adjust <- vec[i]

	  lh_new <- lh
      if("ML50" %in% param_adjust){
        lh_new <- create_lh_list(vbk=lh$vbk, linf=lh$linf, lwa=lh$lwa, lwb=lh$lwb, S50=lh$S50, M50=val_adjust, selex_input="age", maturity_input="length", selex_type=lh$selex_type, dome=lh$dome, binwidth=lh$binwidth, t0=lh$t0, CVlen=lh$CVlen, SigmaC=lh$SigmaC, SigmaI=lh$SigmaI, SigmaR=lh$SigmaR, SigmaF=lh$SigmaF, R0=lh$R0,  h=lh$h, qcoef=lh$qcoef, M=lh$M, F1=lh$F1, Fequil=lh$Fequil, Frate=lh$Frate, Fmax=lh$Fmax, start_ages=min(lh$ages), rho=lh$rho, theta=lh$theta)
      }
      if("M50" %in% param_adjust){
          lh_new <- create_lh_list(vbk=lh$vbk, linf=lh$linf, lwa=lh$lwa, lwb=lh$lwb, S50=lh$S50, M50=val_adjust, selex_input="age", maturity_input="age", selex_type=lh$selex_type, dome=lh$dome, binwidth=lh$binwidth, t0=lh$t0, CVlen=lh$CVlen, SigmaC=lh$SigmaC, SigmaI=lh$SigmaI, SigmaR=lh$SigmaR, SigmaF=lh$SigmaF, R0=lh$R0,  h=lh$h, qcoef=lh$qcoef, M=lh$M, F1=lh$F1, Fequil=lh$Fequil, Frate=lh$Frate, Fmax=lh$Fmax, start_ages=min(lh$ages), rho=lh$rho, theta=lh$theta)
      }
      if("SL50" %in% param_adjust){
          lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=val_adjust, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta)
      }
      if("S50" %in% param_adjust){
          lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=val_adjust, M50=lh_new$ML50, selex_input="age", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta)
      }
      lh_new[[param_adjust]] <- val_adjust

	## run models
	res <- run_LIME(modpath=out, write=TRUE, lh=lh_new, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", LFdist=1, simulation=FALSE, rewrite=rewrite, param_adjust=c("SigmaR","SigmaF"), val_adjust=c(0.7,0.2))
	
	res2 <- run_LBSPR(modpath=out, lh=lh_new, species="sig_sut", input_data=input_data, rewrite=rewrite, simulation=FALSE, write=TRUE)

	Report <- readRDS(file.path(out, "Report.rds"))
	Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
	Inputs <- readRDS(file.path(out, "Inputs.rds"))
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))	

	png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
	plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
	dev.off()	

	png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
	plot_LC(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
	dev.off()
  }

  ## plot likelihood profile
	nll_vec <- rep(NA, length(vec))
	col_vec <- rep("black", length(vec))
	for(i in 1:length(vec)){
	out <- file.path(dir, paste0(param, "_", vec[i]))
		Report <- readRDS(file.path(out, "Report.rds"))
		if(round(max(Report$F_t))==10) col_vec[i] <- "gray"
		nll_vec[i] <- Report$jnll
	}
	png(file.path(dir, "likelihood_profile.png"), width=16, height=10, res=200, units="in")
	plot(x=vec, y=nll_vec, xlab=param, ylab="NLL", xaxs="i", yaxs="i", pch=19, cex=2, xpd=NA, ylim=c(0.8*min(nll_vec), 1.2*max(nll_vec)), col=col_vec)
	dev.off()


}