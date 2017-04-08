run_prof <- function(vec, param, dir, LC, lh, add_years=0, S_l_input=-1, C_t=NULL, I_t=NULL, F_up=5, est_sigma="log_sigma_R", fix_param=FALSE, param_adjust=FALSE, val_adjust=FALSE, data_avail="LC", rerun=TRUE, theta_type=0){

  for(i in 1:length(vec)){
	out <- file.path(dir, paste0(param, "_", vec[i]))
	dir.create(out, showWarnings=FALSE)

	param_adjust <- param
	val_adjust <- vec[i]

	### need to do the sensitivity out here so it also applies to LBSPR
	  lh_new <- lh
	  lh_new[[param_adjust]] <- val_adjust
      if("ML50" %in% param_adjust){
        lh_new <- create_lh_list(vbk=lh$vbk, linf=lh$linf, lwa=lh$lwa, lwb=lh$lwb, S50=lh$SL50, M50=val_adjust, selex_input="length", maturity_input="length", selex_type=lh$selex_type, dome=lh$dome, binwidth=lh$binwidth, t0=lh$t0, CVlen=lh$CVlen, SigmaC=lh$SigmaC, SigmaI=lh$SigmaI, SigmaR=lh$SigmaR, SigmaF=lh$SigmaF, R0=lh$R0,  h=lh$h, qcoef=lh$qcoef, M=lh$M, F1=lh$F1, Fequil=lh$Fequil, Frate=lh$Frate, Fmax=lh$Fmax, start_ages=min(lh$ages), rho=lh$rho, theta=lh$theta, nseasons=1)
      }
      if("M50" %in% param_adjust){
          lh_new <- create_lh_list(vbk=lh$vbk, linf=lh$linf, lwa=lh$lwa, lwb=lh$lwb, S50=lh$SL50, M50=val_adjust, selex_input="length", maturity_input="age", selex_type=lh$selex_type, dome=lh$dome, binwidth=lh$binwidth, t0=lh$t0, CVlen=lh$CVlen, SigmaC=lh$SigmaC, SigmaI=lh$SigmaI, SigmaR=lh$SigmaR, SigmaF=lh$SigmaF, R0=lh$R0,  h=lh$h, qcoef=lh$qcoef, M=lh$M, F1=lh$F1, Fequil=lh$Fequil, Frate=lh$Frate, Fmax=lh$Fmax, start_ages=min(lh$ages), rho=lh$rho, theta=lh$theta, nseasons=1)
      }
      if("SL50" %in% param_adjust){
          lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=val_adjust, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      if("S50" %in% param_adjust){
          lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=val_adjust, M50=lh_new$ML50, selex_input="age", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      if("linf" %in% param_adjust){
      	    lh_new <- create_lh_list(vbk=lh_new$vbk, linf=val_adjust, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=lh_new$SL50, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      if("vbk" %in% param_adjust){
      	    lh_new <- create_lh_list(vbk=val_adjust, linf=lh_new$linf, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=lh_new$SL50, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      if("M" %in% param_adjust){
      	    lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=lh_new$lwa, lwb=lh_new$lwb, S50=lh_new$SL50, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=val_adjust, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      if("lwa" %in% param_adjust){
      	    lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=val_adjust, lwb=lh_new$lwb, S50=lh_new$SL50, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      if("lwb" %in% param_adjust){
      	    lh_new <- create_lh_list(vbk=lh_new$vbk, linf=lh_new$linf, lwa=lh_new$lwa, lwb=val_adjust, S50=lh_new$SL50, M50=lh_new$ML50, selex_input="length", maturity_input="length", selex_type=lh_new$selex_type, dome=lh_new$dome, binwidth=lh_new$binwidth, t0=lh_new$t0, CVlen=lh_new$CVlen, SigmaC=lh_new$SigmaC, SigmaI=lh_new$SigmaI, SigmaR=lh_new$SigmaR, SigmaF=lh_new$SigmaF, R0=lh_new$R0,  h=lh_new$h, qcoef=lh_new$qcoef, M=lh_new$M, F1=lh_new$F1, Fequil=lh_new$Fequil, Frate=lh_new$Frate, Fmax=lh_new$Fmax, start_ages=min(lh_new$ages), rho=lh_new$rho, theta=lh_new$theta, nseasons=1)
      }
      

	years_o <- as.numeric(rownames(LC))
	years_m <- (min(years_o)-add_years):max(years_o)

	input_data <- list("years"=years_m, "LF"=LC, "C_t"=C_t, I_t=I_t)


		if(file.exists(file.path(dir, "LBSPR_results.rds"))==FALSE | rerun==TRUE){
				run <- run_LBSPR(modpath=out, lh=lh_new, species="x", input_data=input_data, rewrite=TRUE, simulation=FALSE)	
				lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))
		}

		if(file.exists(file.path(out, "Report.rds"))==FALSE | rerun==TRUE){
			res <- run_LIME(modpath=out, lh=lh_new, input_data=input_data, est_sigma=est_sigma, data_avail=data_avail, LFdist=1, simulation=FALSE, rewrite=TRUE, S_l_input=S_l_input, F_up=F_up, fix_param=fix_param, param_adjust=param_adjust, val_adjust=val_adjust, theta_type=theta_type)
		}

		if(file.exists(file.path(out, "Report.rds"))){
			if(file.exists(file.path(out, "LBSPR_results.rds"))) lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))
			if(file.exists(file.path(out, "LBSPR_results.rds"))==FALSE) lbspr_res <- NULL
				
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data$years, lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
			dev.off()				

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res)
			dev.off()
		}

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
	plot(x=vec, y=nll_vec, xlab=param, ylab="NLL", xaxs="i", yaxs="i", pch=19, cex=2, xpd=NA, ylim=c(0.95*min(nll_vec), 1.05*max(nll_vec)), col=col_vec)
	dev.off()


}