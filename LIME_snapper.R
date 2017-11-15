rm(list=ls())

## load LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

## load LBSPR package
devtools::install_github("adrianhordyk/LBSPR", build.vignettes=TRUE, dependencies=TRUE)
library(LBSPR)

## load FishLife
devtools::install_github("james-thorson/FishLife")
library(FishLife)

###################################
## Directories
###################################
main_dir <- "C:\\merrill\\LIME_application"

R_dir <- file.path(main_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

###################################
## Species directories
###################################

sp_dir <- file.path(main_dir, "costa_rican_snapper")
dir.create(sp_dir, showWarnings=FALSE)

data_dir <- file.path(sp_dir, "data")
dir.create(data_dir, showWarnings=FALSE)

res_dir <- file.path(sp_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

figs_dir <- file.path(sp_dir, "figures")
dir.create(figs_dir, showWarnings=FALSE)

###################################
## Life history info
###################################

## add life history information for species chosen for assessment
## siganus sutor
plist <- readRDS(file.path(data_dir, "CRSNAP_life_history_annual.rds"))
plist_q <- readRDS(file.path(data_dir, "CRSNAP_life_history_quarterly.rds"))
plist_m <- readRDS(file.path(data_dir, "CRSNAP_life_history_monthly.rds"))

###################################
## Load data
###################################

## LC separate by gears
LC_bline <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_bottomlongline.rds"))
LCq_bline <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_bottomlongline_quarterly.rds"))

LC_gnet <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_gillnet.rds"))
LCq_gnet <- readRDS(file.path(data_dir, "Lutjanus_guttatus_LF_gillnet_quarterly.rds"))

#########################################
## run assessments -- life history fixed
#########################################

years_o <- as.numeric(unique(c(rownames(LC_bline), rownames(LC_gnet))))
years_t <- min(years_o):max(years_o)

quarters_o <- as.numeric(unique(c(rownames(LCq_bline), rownames(LCq_gnet))))
quarters_t <- min(quarters_o):max(quarters_o)

q_all <- as.vector(sapply(1:length(years_t), function(x) rep(years_t[x],4)))
true_q <- q_all[quarters_o]


LC <- array(0, dim=c(length(years_t), ncol(LC_bline), 2))
rownames(LC) <- years_t
colnames(LC) <- colnames(LC_bline)
LC[which(rownames(LC) %in% rownames(LC_bline)),,1] <- LC_bline
LC[which(rownames(LC) %in% rownames(LC_gnet)),,2] <- LC_gnet

LCq <- array(0, dim=c(length(quarters_t), ncol(LC_bline), 2))
rownames(LCq) <- quarters_t
colnames(LCq) <- colnames(LC_bline)
LCq[which(rownames(LCq) %in% rownames(LCq_bline)),,1] <- LCq_bline
LCq[which(rownames(LCq) %in% rownames(LCq_gnet)),,2] <- LCq_gnet

input_data_y <- list("years"=years_t, "LF"=LC)
input_data_q <- list("years"=quarters_t, "LF"=LCq)

# LCtest <- array(0, dim=c(1,ncol(LC[,,1]),1))
# LCtest[,,1] <- LC[1,,1]
# rownames(LCtest) <- years_t[1]
# colnames(LCtest) <- colnames(LC)

# input_data_test <- list("years"=years_t[1], "LF"=LCtest)

src <- "C:\\merrill\\LIME\\src"
setwd(src)
dyn.unload("LIME.dll")
file.remove(paste0("LIME",c(".dll",".o"), sep=""))
compile("LIME.cpp")
dyn.load(file.path(src, "LIME.dll"))

## annual LIME
out <- file.path(res_dir, "LCy")
dir.create(out, showWarnings=FALSE)

	res <- run_LIME(modpath=out, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, L_a_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

	modpath=out
	lh=plist
	input_data=input_data_y
	est_sigma="log_sigma_R"
	data_avail="LC"
	itervec=NULL
	rewrite=FALSE
	simulation=FALSE
	f_true=FALSE
	C_opt=0
	LFdist=1
	param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI")
	val_adjust=c(0.2,0.737,0.2,0.2)
	F_up=10
	S_l_input=-1
	L_a_input=-1
	fix_param=FALSE
	theta_type=1
	randomR=TRUE
	SigRpen=1
	fix_param_t=FALSE
	iter=1

	df <- readRDS(file.path(out, "check_convergence.rds"))
	Report <- readRDS(file.path(out, "Report.rds"))
	Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
	Inputs <- readRDS(file.path(out, "Inputs.rds"))
	
	png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
	plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_y$years, lc_years=rownames(input_data_y$LF), LBSPR=NULL, true_years=years_t, lh=plist)
	dev.off()	
	
	png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
	plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_y$LF), LBSPR=NULL, ylim=c(0,0.3))
	dev.off()


############### weighted length comps ##########################
## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

	png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
	plot_LCfits(Inputs=input_data_y, Report=NULL, true_lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, ylim=c(0,0.2), ML50=plist$ML50, SL50=lbspr_res$SL50)
	dev.off()

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_highN")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_up, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_highN <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_bottomlongline")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_bline, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_bl <- readRDS(file.path(out, "LBSPR_results.rds"))

## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy_altM")
dir.create(out, showWarnings=FALSE)

plist_new <- plist
plist_new$M <- 0.43
	run <- run_LBSPR(modpath=out, lh=plist_new, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res_altM <- readRDS(file.path(out, "LBSPR_results.rds"))




## quarterly LIME
out <- file.path(res_dir, "LCq")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			base_rep <- Report

			## previous F = 0.34
			calc_ref(ages=plist_q$ages, Mat_a=plist_q$Mat_a, W_a=plist_q$W_a, M=plist_q$M, S_a=Report$S_a, F=0.34/4)

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_q$years, lc_years=rownames(input_data_q$LF), LBSPR=lbspr_res, true_years=years_t, lh=plist_q)
			dev.off()	
			 F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_q$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
			 Report$F_t/F40
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=TRUE)[c(length(Report$F_y), length(Report$F_y)+1)]
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=FALSE)[c(length(Report$SPR_t), length(Report$SPR_t)+1)]
    		pthresh <- pnorm(0.3, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		ptarget <- pnorm(0.45, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]) - pnorm(0.35, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=true_q, LBSPR=NULL, ylim=c(0,0.3), ML50=plist_q$ML50, SL50=Report$S50)
			dev.off()




### weighted length comp - counts only - all years - gear1
out <- file.path(res_dir, "LCq_bottomlongline")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_bline_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))
			 
			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=FALSE)[c(length(Report$SPR_t), length(Report$SPR_t)+1)]
    		pthresh <- pnorm(0.3, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		ptarget <- pnorm(0.45, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]) - pnorm(0.35, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_bline_q$years, lc_years=rownames(input_data_bline_q$LF), LBSPR=lbspr_res_bl, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_bline_q$LF), LBSPR=NULL, ylim=c(0,0.3), ML50=plist$ML50, SL50=Report$S50)
			dev.off()

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCq_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_up_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=FALSE)[c(length(Report$SPR_t), length(Report$SPR_t)+1)]
    		pthresh <- pnorm(0.3, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		ptarget <- pnorm(0.45, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]) - pnorm(0.35, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		pup <- (1-pnorm(0.4, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up_q$years, lc_years=rownames(input_data_up_q$LF), LBSPR=lbspr_res_highN, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up_q$LF), LBSPR=NULL, ylim=c(0,0.3), ML50=plist$ML50, SL50=Report$S50)
			dev.off()


### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCq_altM")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_up_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI","M"), val_adjust=c(0.2,0.737,0.2,0.2,0.43), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			 sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
  			sd[,2][which(is.na(sd[,2]))] <- 0
    		read_sdreport(sd, log=FALSE)[c(length(Report$SPR_t), length(Report$SPR_t)+1)]
    		pthresh <- pnorm(0.3, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		ptarget <- pnorm(0.45, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]) - pnorm(0.35, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)])
    		pup <- (1-pnorm(0.4, Report$SPR_t[length(Report$SPR_t)], sd[,2][length(Report$SPR_t)]))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up_q$years, lc_years=rownames(input_data_up_q$LF), LBSPR=lbspr_res_highN, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up_q$LF), LBSPR=NULL, ylim=c(0,0.3), ML50=plist$ML50, SL50=Report$S50)
			dev.off()



## dome-shaped selex
plist_q_lowdome <- create_lh_list(vbk=plist$vbk, linf=plist$linf, t0=plist$t0, M=plist$M, lwa=plist$lwa, lwb=plist$lwb, S50=base_rep$S50, S95=NULL, M50=plist$ML50, M95=plist$ML95, selex_input="length", maturity_input="length", nseasons=4, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=51)
plist_q_highdome <- create_lh_list(vbk=plist$vbk, linf=plist$linf, t0=plist$t0, M=plist$M, lwa=plist$lwa, lwb=plist$lwb, S50=base_rep$S50, S95=NULL, M50=plist$ML50, M95=plist$ML95, selex_input="length", maturity_input="length", nseasons=4, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=28)


par(mfrow=c(1,1), mar=c(5,5,4,4))
plot(base_rep$S_l, type="o", lwd=2, pch=19, cex.axis=2, xlab="Length (cm)", ylab="Proportion", cex.lab=2)
lines(plist_q_lowdome$S_l, type="o", lwd=2, lty=2, col="steelblue", pch=19)
lines(plist_q_highdome$S_l, type="o", lwd=2, lty=3, col="tomato", pch=19)

## low-dome
out <- file.path(res_dir, "LCq_lowdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_q, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=plist_q_lowdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_q$years, lc_years=rownames(input_data_q$LF), LBSPR=lbspr_res, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_q$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


## high-dome
out <- file.path(res_dir, "LCq_highdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_q, input_data=input_data_q, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2), F_up=10, S_l_input=plist_q_highdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_q$years, lc_years=rownames(input_data_q$LF), LBSPR=lbspr_res, lh=plist_q, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_q$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()



###########################
## likelihood profiles
###########################

## sigma R
out <- file.path(res_dir, "sens_SigmaR")
dir.create(out, showWarnings=FALSE)

sigR_vec <- seq(0.05, 0.95, length=10)

for(i in 1:length(sigR_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma=FALSE, data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR"), val_adjust=c(0.2,sigR_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
sigR_like <- rep(NA, length(sigR_vec))
col_conv <- rep("black", length(sigR_vec))
for(i in 1:length(sigR_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	sigR_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "sigR_likeprof.png"), height=10, width=14, res=200, units="in")
plot(sigR_vec, sigR_like, pch=19, cex=2, xlab="Fixed SigmaR", ylab="NLL", col=col_conv)
points(rep1$sigma_R, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


## check fishbase and natural mortality toolfor life history distributions
proper_name <- "Lutjanus guttatus"
genus <- strsplit(proper_name, " ")[[1]][1]

# growth2 <- popgrowth(species_list(Genus=genus))
# saveRDS(growth2, file.path(data_dir, "genus_info.rds"))

## linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

linf_genus <- growth2$Loo
linf_vec <- seq(min(linf_genus), max(linf_genus), length=10)

for(i in 1:length(linf_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","linf"), val_adjust=c(0.2,0.737,linf_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
linf_like <- rep(NA, length(linf_vec))
col_conv <- rep("black", length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	linf_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "linf_likeprof.png"), height=10, width=14, res=200, units="in")
plot(linf_vec, linf_like, pch=19, cex=2, xlab="Assumed Linf", ylab="NLL", col=col_conv)
points(plist_q$linf, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

vbk_genus <- growth2$K
vbk_vec <- seq(min(vbk_genus), max(vbk_genus), length=10)

for(i in 1:length(vbk_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","vbk"), val_adjust=c(0.2,0.737,vbk_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
vbk_like <- rep(NA, length(vbk_vec))
col_conv <- rep("black", length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	vbk_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "vbk_likeprof.png"), height=10, width=14, res=200, units="in")
plot(vbk_vec, vbk_like, pch=19, cex=2, xlab="Assumed vbk", ylab="NLL", col=col_conv)
points(plist_q$vbk, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## M
out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

M_tool <- c(0.29, 0.23, 0.25, 0.25, 0.13, 0.16, 0.32, 0.34, 0.48, 0.41)
M_vec <- seq(min(M_tool), max(M_tool), length=10)

for(i in 1:length(M_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

			res <- run_LIME(modpath=out2, lh=plist_q, input_data=input_data_q, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI", "M"), val_adjust=c(0.2,0.737,0.2,0.2, M_vec[i]), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
M_like <- rep(NA, length(M_vec))
col_conv <- rep("black", length(M_vec))
for(i in 1:length(M_vec)){
	inp <- readRDS(file.path(out, i, "Inputs.rds"))
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	M_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
png(file.path(figs_dir, "M_likeprof.png"), height=10, width=14, res=200, units="in")
plot(M_vec, M_like, pch=19, cex=2, xlab="Assumed M", ylab="NLL", col=col_conv)
points(plist_q$M*4, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


png(file.path(figs_dir, "like_prof_compare.png"), height=10, width=15, units="in", res=200)

mat <- matrix(c(1,2,3,4,
				1,2,3,4,
				5,6,7,8,
				5,6,7,8,
				9,9,9,9,
				10,11,12,13,
				10,11,12,13), nrow=7, ncol=4, byrow=TRUE)
nf <- layout(mat)
layout.show(nf)

par(mar=c(0,4,0,0), omi=c(0.5,0.5,0.5,0.5))
rep1 <- readRDS(file.path(res_dir, "LCq", "Report.rds"))
growth2 <- readRDS(file.path(data_dir, "genus_info.rds"))

## using input AgeMax=4, Linf=36.2, vbk=0.87, Age at maturity(years)=1
M_tool <- c(0.29, 0.23, 0.25, 0.25, 0.13, 0.16, 0.32, 0.34, 0.48, 0.41)
M_toUse <- median(M_tool)
M_h <- hist(M_tool, plot=FALSE)
M_x <- seq(min(M_tool), max(M_tool), length=60)
M_y <- dnorm(M_x, mean=mean(M_tool), sd=sd(M_tool))
M_y2 <- (M_y*diff(M_h$mids[1:2])*length(M_tool))/sum((M_y*diff(M_h$mids[1:2])*length(M_tool)))
plot(M_x, M_y2, col="black", lwd=2, type="l", xaxt="n", xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
rug(M_tool, ticksize=0.15, lwd=2)
# M_d <- rtnorm(1, mean=mean(M_tool), sd=sd(M_tool), lower=min(M_tool), upper=max(M_tool))
# abline(v=M_d, lty=2, col="red", lwd=2)
abline(v=plist$M, lwd=2, col="blue")
mtext(side=3, "M", cex=2.5, line=1)
mtext("Probability", side=2, line=3, cex=2)
print.letter("(a)", xy=c(0.9,0.95), cex=2)

linf_genus <- growth2$Loo
linf_h <- hist(linf_genus, plot=FALSE)
linf_x <- seq(min(linf_genus), max(linf_genus), length=60)
linf_y <- dnorm(linf_x, mean=mean(linf_genus), sd=sd(linf_genus))
linf_y2 <- linf_y*diff(linf_h$mids[1:2])*length(linf_genus)/sum(linf_y*diff(linf_h$mids[1:2])*length(linf_genus))
plot(linf_x, linf_y2, col="black", lwd=2, type="l", xaxt="n",xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
rug(linf_genus, ticksize=0.15, lwd=2)
# linf_d <- rtnorm(1, mean=mean(linf_genus), sd=sd(linf_genus), lower=min(linf_genus), upper=max(linf_genus))
# abline(v=linf_d, lty=2, col="red", lwd=2)
abline(v=plist$linf, lwd=2, col="blue")
mtext(side=3, "Linf", cex=2.5, line=1)
print.letter("(b)", xy=c(0.9,0.95), cex=2)

vbk_genus <- growth2$K
vbk_h <- hist(vbk_genus, plot=FALSE)
vbk_x <- seq(min(vbk_genus), max(vbk_genus), length=60)
vbk_y <- dnorm(vbk_x, mean=mean(vbk_genus), sd=sd(vbk_genus))
vbk_y2 <- vbk_y*diff(vbk_h$mids[1:2])*length(vbk_genus)/sum(vbk_y*diff(vbk_h$mids[1:2])*length(vbk_genus))
plot(vbk_x, vbk_y2, col="black", lwd=2, type="l", xaxt="n",xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
rug(vbk_genus, ticksize=0.15, lwd=2)
# vbk_d <- rtnorm(1, mean=mean(vbk_genus), sd=sd(vbk_genus), lower=min(vbk_genus), upper=max(vbk_genus))
# abline(v=vbk_d, lty=2, col="red", lwd=2)
abline(v=plist$vbk, lwd=2, col="blue")
mtext(side=3, "k", cex=2.5, line=1)
print.letter("(c)", xy=c(0.9,0.95), cex=2)

sigR_x <- seq(0, 1.5, length=60)
sigR_y <- dnorm(sigR_x, mean=mean(0.737), sd=0.353)/sum(dnorm(sigR_x, mean=mean(0.737), sd=0.353))
plot(sigR_x, sigR_y, col="black", lwd=2, type="l", xaxt="n",xlab="", ylab="", xaxs="i", yaxs="i", cex.axis=2)
# sigR_d <- rtnorm(1, mean=mean(sigR_genus), sd=sd(sigR_genus), lower=min(sigR_genus), upper=max(sigR_genus))
# abline(v=sigR_d, lty=2, col="red", lwd=2)
abline(v=0.737, lwd=2, col="blue")
mtext(side=3, "SigmaR", cex=2.5, line=1)
print.letter("(d)", xy=c(0.9,0.95), cex=2)



out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

M_vec <- seq(min(M_tool), max(M_tool), length=10)

M_like <- rep(NA, length(M_vec))
col_conv <- rep("black", length(M_vec))
for(i in 1:length(M_vec)){
	inp <- readRDS(file.path(out, i, "Inputs.rds"))
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	M_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv[i] <- "gray"
}
i1 <- 1
i2 <- 8
col_conv[i1] <- "orange"
col_conv[i2] <- "orangered"
plot(x=M_vec, y=M_like, pch=19, cex=3, col=col_conv, xpd=NA, cex.axis=2, xlab="", ylab="", yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(M_like),max(M_like)))[-length(pretty(c(min(M_like),max(M_like))))])
points(plist_q$M*4, rep1$jnll, pch=19, cex=3, col="blue")
mtext(side=2, "NLL", cex=2, line=3)
print.letter("(e)", xy=c(0.9,0.95), cex=2)

### linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

linf_genus <- growth2$Loo
linf_vec <- seq(min(linf_genus), max(linf_genus), length=10)

linf_like <- rep(NA, length(linf_vec))
col_conv_linf <- rep("black", length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	linf_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv_linf[i] <- "gray"
}
i1 <- which(col_conv_linf=="black")[1]
i2 <- 8
col_conv_linf[i1] <- "orange"
col_conv_linf[i2] <- "orangered"
plot(x=linf_vec, y=linf_like, pch=19, cex=3, col=col_conv_linf, xpd=NA, cex.axis=2, xlab="", ylab="", yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(linf_like),max(linf_like)))[-length(pretty(c(min(linf_like),max(linf_like))))])
points(plist_q$linf, rep1$jnll, pch=19, cex=3, col="blue")
print.letter("(f)", xy=c(0.9,0.95), cex=2)

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

vbk_genus <- growth2$K
vbk_vec <- seq(min(vbk_genus), max(vbk_genus), length=10)

vbk_like <- rep(NA, length(vbk_vec))
col_conv_vbk <- rep("black", length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	vbk_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv_vbk[i] <- "gray"
}
i1 <- 2
i2 <- 5
col_conv_vbk[i1] <- "orange"
col_conv_vbk[i2] <- "orangered"
plot(x=vbk_vec, y=vbk_like, pch=19, cex=3, col=col_conv_vbk, xpd=NA, cex.axis=2, xlab="", ylab="", yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(vbk_like),max(vbk_like)))[-length(pretty(c(min(vbk_like),max(vbk_like))))])
points(plist_q$vbk, rep1$jnll, pch=19, cex=3, col="blue")
print.letter("(g)", xy=c(0.9,0.95), cex=2)


## sigma R
out <- file.path(res_dir, "sens_SigmaR")
dir.create(out, showWarnings=FALSE)

sigR_vec <- seq(0.05, 0.95, length=10)

sigR_like <- rep(NA, length(sigR_vec))
col_conv_sigR <- rep("black", length(sigR_vec))
for(i in 1:length(sigR_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	sdrep <- readRDS(file.path(out, i, "Sdreport.rds"))
	sigR_like[i] <- rep$jnll
	if(file.exists(file.path(out, i, "high_final_gradient.txt")) | all(is.na(sdrep))) col_conv_sigR[i] <- "gray"
}
i1 <- 4
i2 <- 7
col_conv_sigR[i1] <- "orange"
col_conv_sigR[i2] <- "orangered"
plot(x=sigR_vec, y=sigR_like, pch=19, cex=3, col=col_conv_sigR, xpd=NA, cex.axis=2, xlab="", ylab="", xlim=c(min(sigR_x), max(sigR_x)), yaxt="n")
axis(2, cex.axis=2, at=pretty(c(min(sigR_like),max(sigR_like)))[-length(pretty(c(min(sigR_like),max(sigR_like))))])
points(rep1$sigma_R, rep1$jnll, pch=19, cex=3, col="blue")
print.letter("(h)", xy=c(0.9,0.95), cex=2)


plot(x=1,y=1,type="n", axes=F, ann=F)
### SPR
out <- file.path(res_dir, "sens_M")
i1 <- 1
i2 <- 8
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=2, "SPR", line=3, cex=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(i)", xy=c(0.9,0.95), cex=2)

out <- file.path(res_dir, "sens_linf")
i1 <- which(col_conv_linf=="black")[1]
i2 <- 8
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(j)", xy=c(0.9,0.95), cex=2)

out <- file.path(res_dir, "sens_vbk")
i1 <- 2
i2 <- 5
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(k)", xy=c(0.9,0.95), cex=2)

out <- file.path(res_dir, "sens_SigmaR")
i1 <- 4
i2 <- 7
replow <- readRDS(file.path(out, i1, "Report.rds"))
rephigh <- readRDS(file.path(out, i2, "Report.rds"))
plot(x=years_t, y=unique(rep1$SPR_t), col="blue", lwd=2, type="l", ylim=c(0,1.1), xlab="", ylab="", cex.axis=2)
lines(x=years_t, y=unique(replow$SPR_t), col="orange", lwd=2)
lines(x=years_t, y=unique(rephigh$SPR_t), col="orangered", lwd=2)
abline(h=0.4, lty=2)
mtext(side=1, "Year", line=3, cex=2)
print.letter("(l)", xy=c(0.9,0.95), cex=2)

dev.off()


## model fits
out <- file.path(res_dir, "LCq")
dir.create(out, showWarnings=FALSE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			base_rep <- Report

out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

	obs1 <- Inputs$Data$LF
	all_mos <- 1:(length(years_t)*4)
	mo_yrs <- as.vector(sapply(1:length(years_t), function(x) rep(years_t[x], 4)))
	choose_mo <- rep(NA, length(years_o))
	for(i in 1:length(years_o)){
		sub_mo <- quarters_o[which(quarters_o %in% all_mos[which(mo_yrs==years_o[i])])]
		sub_mo2 <- which(rownames(obs1) %in% sub_mo)
		choose_mo[i] <- sub_mo2[1]
	}

	obs2 <- input_data_y$LF

png(file.path(figs_dir, "CompareModelFits.png"), height=15, width=8, res=200, units="in")
par(mfcol=c(length(years_o),2), mar=c(0,0,0,0), omi=c(1,1,1,1))
for(i in 1:length(years_o)){
	barplot(obs1[which(rownames(obs1)==choose_mo[i]),]/sum(obs1[which(rownames(obs1)==choose_mo[i]),]), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,0.2), col=gray(0.6), border=NA, space=0)
	lines(base_rep$plb[choose_mo[i],], lwd=4, col="blue")
	if(i==length(years_o)) axis(1, cex.axis=2)
	axis(2, cex.axis=2, at=c(0,0.15), las=2)
	abline(v=plist$ML50, lty=2)
	box()
}
for(i in 1:length(years_o)){
	barplot(obs2[which(rownames(obs2)==years_o[i]),]/sum(obs2[which(rownames(obs2)==years_o[i]),]), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,0.2), col=gray(0.6), border=NA, space=0)
	lines(lbspr_res$pLF[,i], lwd=4, col="red")
	if(i==length(years_o)) axis(1, cex.axis=2)
	print.letter(years_o[i], xy=c(0.8,0.8), cex=3)
	abline(v=plist$ML50, lty=2)
	box()
}
mtext(side=1, "Length bin (cm)", cex=2, line=4, outer=TRUE)
mtext(side=2, "Proportions", cex=2, line=5, outer=TRUE)
dev.off()























####################################
### MONTE CARLO
####################################
n <- 1000
dist <- data.frame("linf"=rlnorm(n,log(plist$linf),0.25), "vbk"=rlnorm(n,log(plist$vbk),0.25), "M"=rlnorm(n,log(plist$M),0.25))


par(mfrow=c(1,3), omi=c(0.5,1,0.5,0.5))
hist(dist$linf, col="goldenrod", xlab="", ylab="", main="", cex.axis=2, cex.lab=2, border=NA, yaxs="i")
mtext(side=3, "Asymptotic length (Linf)", cex=2)
abline(v=plist$linf, lwd=2, col="blue")
box()

hist(dist$vbk, col="goldenrod", xlab="", ylab="", main="", cex.axis=2, cex.lab=2, border=NA, yaxs="i")
mtext(side=3, "von Bertalanffy (k)", cex=2)
abline(v=plist$vbk, lwd=2, col="blue")
box()

hist(dist$M, col="goldenrod", xlab="", ylab="", main="", cex.axis=2, cex.lab=2, border=NA, yaxs="i")
mtext(side=3, "Natural mortality (M)", cex=2)
abline(v=plist$M, lwd=2, col="blue")
box()
mtext("Frequency", side=2, outer=TRUE, cex=2)



############################
### life history figures

linf_genus <- growth2$Loo
linf_h <- hist(linf_genus, plot=FALSE)
linf_x <- seq(min(linf_genus), max(linf_genus), length=60)
linf_y <- dnorm(linf_x, mean=mean(linf_genus), sd=sd(linf_genus))
linf_y2 <- linf_y*diff(linf_h$mids[1:2])*length(linf_genus)
plot(linf_x, linf_y2, col="blue", lwd=2, type="l")
rug(linf_genus)
linf_d <- rtnorm(1, mean=mean(linf_genus), sd=sd(linf_genus), lower=min(linf_genus), upper=max(linf_genus))
abline(v=linf_d, lty=2, col="red", lwd=2)
abline(v=linf_toUse, lwd=2)

vbk_genus <- growth2$K
vbk_h <- hist(vbk_genus, plot=FALSE)
vbk_x <- seq(min(vbk_genus), max(vbk_genus), length=60)
vbk_y <- dnorm(vbk_x, mean=mean(vbk_genus), sd=sd(vbk_genus))
vbk_y2 <- vbk_y*diff(vbk_h$mids[1:2])*length(vbk_genus)
plot(vbk_x, vbk_y2, col="blue", lwd=2, type="l")
rug(vbk_genus)
vbk_d <- rtnorm(1, mean=mean(vbk_genus), sd=sd(vbk_genus), lower=min(vbk_genus), upper=max(vbk_genus))
abline(v=vbk_d, lty=2, col="red", lwd=2)
abline(v=vbk_toUse, lwd=2)


## using input AgeMax=22, Linf=64, vbk=0.21, Age at maturity(years)=4, asymptotic weight = 2750g.
M_tool <- c(0.29, 0.23, 0.25, 0.25, 0.13, 0.16, 0.32, 0.34, 0.48, 0.41)
M_toUse <- median(M_tool)
M_h <- hist(M_tool, plot=FALSE)
M_x <- seq(min(M_tool), max(M_tool), length=60)
M_y <- dnorm(M_x, mean=mean(M_tool), sd=sd(M_tool))
M_y2 <- M_y*diff(M_h$mids[1:2])*length(M_tool)
plot(M_x, M_y2, col="blue", lwd=2, type="l")
rug(M_tool)
M_d <- rtnorm(1, mean=mean(M_tool), sd=sd(M_tool), lower=min(M_tool), upper=max(M_tool))
abline(v=M_d, lty=2, col="red", lwd=2)
abline(v=M_toUse, lwd=2)
