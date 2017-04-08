rm(list=ls())

## load LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

## load LBSPR package
library(LBSPR)

## load Rfishbase
library(rfishbase)

###################################
## Directories
###################################
main_dir <- "C:\\Git_Projects\\LIME_application"

R_dir <- file.path(main_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

region_dir <- file.path(main_dir, "kenyan_reef_fish")

######################################
## Load full dataset
######################################

data_dir <- file.path(region_dir, "data")
dir.create(data_dir, showWarnings=FALSE)

alldata <- read.csv(file.path(data_dir, "LengthDataChecked.csv"), stringsAsFactors=FALSE)

######################################
## Species to assess
######################################
### Goal:
## Lethrinus lentjan
## Siganus sutor
## Leptoscarus vaigiensis

## names of species to search for in full dataset for assessment
name <- c("siganus sutor") 
proper_name <- c("Siganus sutor")

## name to save files
save_name <- c("Siganus_sutor")

## create directory for this assessment
stock_dir <- file.path(region_dir, save_name)
dir.create(stock_dir, showWarnings=FALSE)

## data to be used for the assessment
data <- alldata[which(alldata$Name %in% name),]

## years in the dataset
years <- unique(data$Year)[order(unique(data$Year))]
all_months <- as.vector(sapply(1:length(years), function(x) paste0(1:12, "/", years[x])))

## months in the dataset
data$Month <- sapply(1:nrow(data), function(x){
	mo <- strsplit(data$Date[x],"/")[[1]][2]
	if(mo=="?") return(NA)
	if(mo %in% paste0(0,1:9)) mo <- strsplit(mo,"0")[[1]][2]
	if(as.numeric(mo) > 13) mo <- strsplit(data$Date[x],"/")[[1]][1]
	return(as.numeric(mo))
})

data <- data[which(data$Month<13),]

## order by month
data_new <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
colnames(data_new) <- colnames(data)
data_new <- as.data.frame(data_new)
ilast <- 1
for(y in 1:length(years)){
	sub <- data[which(data$Year==years[y]),]
	mos <- unique(sub$Month)[order(unique(sub$Month))]
	for(m in 1:length(mos)){
		sub2 <- sub[which(sub$Month==mos[m]),]
		index <- ilast:(ilast+nrow(sub2)-1)
		data_new[index,] <- sub2
		ilast <- max(index)+1
		rm(sub2)
	}
}

data_new$MoYr <- sapply(1:nrow(data_new), function(x) paste0(data_new$Month[x], "/", data_new$Year[x]))
data_new$MoTotal <- sapply(1:nrow(data_new), function(x) which(all_months==data_new$MoYr[x]))
moyrs <- unique(data_new$MoYr)
mototal <- unique(data_new$MoTotal)

## number of sampling days per year
ndays <- sapply(1:length(years), function(x){
	sub <- data_new[which(data_new$Year==years[x]),]
	return(length(unique(sub$Date)))
})
d_percentiles <- quantile(ndays, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))

ndays_mo <- sapply(1:length(mototal), function(x){
	sub <- data_new[which(data_new$MoTotal==mototal[x]),]
	return(length(unique(sub$Date)))
})
dm_percentiles <- quantile(ndays_mo, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))


## upper 50th percentile of days sampled per year
years_up <- years[which(ndays > median(ndays))]
years_mid <- years[which(ndays > quantile(ndays, prob=0.25))]

months_up <- mototal[which(ndays_mo > median(ndays_mo))]
months_mid <- mototal[which(ndays_mo > quantile(ndays_mo, prob=0.25))]

## number of samples per year
nsamp <- sapply(1:length(years), function(x){
	sub <- data_new[which(data_new$Year==years[x]),]
	return(nrow(sub))
})
nsamp_up <- nsamp[which(years %in% years_up)]
nsamp_mid <- nsamp[which(years %in% years_mid)]

######################################
## Species life history
######################################

## check fishbase
lw <- length_weight(proper_name)
lwa_med <- median(lw$a, na.rm=TRUE)
lwb_med <- median(lw$b, na.rm=TRUE)

growth <- popgrowth(proper_name)
linf_med <- median(growth$TLinfinity, na.rm=TRUE)
vbk_med <- median(growth$K, na.rm=TRUE)
winf_med <- median(growth$Winfinity, na.rm=TRUE)
temp_med <- median(growth$Temperature)

## use values from Mombasa area
linf_toUse <- 42.7
vbk_toUse <- 0.87
t0_toUse <- -0.24

## use natural mortality tool
M_toUse <- 1.42

## add life history information for species chosen for assessment
## siganus sutor
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=1, binwidth=1)
plist_m <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1)

######################################
## Length composition data_new
######################################

## length composition
bins <- seq(plist$binwidth, 1.5*plist$linf, by=plist$binwidth)
LC <- matrix(0, nrow=length(years), ncol=length(bins))
	rownames(LC) <- years
	colnames(LC) <- bins
for(y in 1:length(years)){
	sub <- data_new[which(data_new$Year==years[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$value<=bins[x] & sub$value>(bins[x]-bins[1]))))
	LC[y,] <- freq
}
LC_up <- LC[which(years %in% years_up),]
LC_mid <- LC[which(years %in% years_mid),]

## put in chronological order
LCm <- matrix(0, nrow=length(mototal), ncol=length(bins))
	rownames(LCm) <- mototal
	colnames(LCm) <- bins
for(y in 1:length(mototal)){
	sub <- data_new[which(data_new$MoTotal==mototal[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$value<=bins[x] & sub$value>(bins[x]-bins[1]))))
	LCm[y,] <- freq
}

LC_lbspr <- t(LC)
colnames(LC_lbspr) <- NULL
rownames(LC_lbspr) <- bins - plist$binwidth/2

write.csv(LC, file.path(data_dir, paste0(save_name, "_LC.csv")))
write.csv(LCm, file.path(data_dir, paste0(save_name, "_LCmonthly.csv")))
write.csv(LC_up, file.path(data_dir, paste0(save_name, "_LC_up.csv")))
write.csv(LC_mid, file.path(data_dir, paste0(save_name, "_LC_mid.csv")))
write.table(LC_lbspr, file.path(data_dir, paste0(save_name, "_LC_forLBSPR.csv")), col.names=FALSE, sep=",")

## check all observations are accounted for in length composition
sum(LC)==nrow(data)
sum(LCm)==nrow(data)

## length comp proportion
LCprop <- LC/rowSums(LC)
write.csv(LCprop, file.path(data_dir, paste0(save_name, "_LCprop.csv")))

LCmprop <- LCm/rowSums(LCm)
write.csv(LCmprop, file.path(data_dir, paste0(save_name, "_LCmonthlyprop.csv")))

## plot length comp
plot_LCdata(lc_list=list("LF"=LCm), lc_years=mototal, gears=NULL, ylim=c(0,0.3))

LCm_up <- LCm[which(rownames(LCm) %in% months_up),]
LCm_mid <- LCm[which(rownames(LCm) %in% months_mid),]

plot_LCdata(lc_list=list("LF"=LCm_up), lc_years=months_up, gears=NULL, ylim=c(0,0.3))
plot_LCdata(lc_list=list("LF"=LCm_mid), lc_years=months_mid, gears=NULL, ylim=c(0,0.3))


## samples by gear annually
gears_raw <- unique(data$Fishgear)[order(unique(data$Fishgear))]
gearmat_raw <- matrix(NA, nrow=length(years), ncol=length(gears_raw))
rownames(gearmat_raw) <- years
colnames(gearmat_raw) <- gears_raw
for(i in 1:length(years)){
	sub1 <- data[which(data$Year==years[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Fishgear==gears_raw[j]),]
		gearmat_raw[i,j] <- nrow(sub2)
	}
}
## in order of highest to lowest samples
gearmat <- gearmat_raw[,rev(order(colSums(gearmat_raw)))]
gear_totals <- colSums(gearmat)[rev(order(colSums(gearmat)))]
gear_prop <- gear_totals/nrow(data)
gears <- colnames(gearmat)

## samples in each length bin by gear
lc_gear <- list()
for(g in 1:length(gears)){
	sub <- data[which(data$Fishgear==gears[g]),]
	subyears <- unique(sub$Year)[order(unique(sub$Year))]
	mat <- matrix(0, nrow=length(subyears), ncol=length(bins))
			rownames(mat) <- subyears
			colnames(mat) <- bins
	for(y in 1:length(subyears)){
		sub2 <- sub[which(sub$Year==subyears[y]),]
		freq <- sapply(1:length(bins), function(x) length(which(sub2$value<=bins[x] & sub2$value>(bins[x]-bins[1]))))
		mat[y,] <- freq
	}
	if(sum(mat)!=nrow(sub)) stop("Not all lengths included in LC")
	lc_gear[[g]] <- mat
}
names(lc_gear) <- gears

## proportion of samples in each length bin by gear
lcprop_gear <- lapply(1:length(lc_gear), function(x) lc_gear[[x]]/rowSums(lc_gear[[x]]))
names(lcprop_gear) <- gears


####################################
## run assessments
####################################

res_dir <- file.path(stock_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LC))
years_t <- min(years_o):max(years_o)

months_o <- as.numeric(rownames(LCm))
months_t <- min(months_o):max(months_o)

input_data_y <- list("years"=years_t, "LF"=LC)
input_data_m <- list("years"=months_t, "LF"=LCm)

############### weighted length comps ##########################

out <- file.path(res_dir, "LC_counts")
dir.create(out, showWarnings=FALSE)

				run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
				lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

			res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_y$years, lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data$LF), LBSPR=lbspr_res, ylim=c(0,0.3))
			dev.off()


### weighted length comp - counts only - all months
out <- file.path(res_dir, "LC_counts_monthly")
dir.create(out, showWarnings=FALSE)


				run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
				lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=NULL, lh=plist_m)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only
out <- file.path(res_dir, "LC_counts_monthly_highsampling")
dir.create(out, showWarnings=FALSE)

months_o <- as.numeric(rownames(LCm_up))
months_t <- min(months_o):max(months_o)

input_data_m <- list("years"=months_t, "LF"=LCm_up)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=NULL, lh=plist_m)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


base_rep <- readRDS(file.path(res_dir, "LC_counts_monthly", "Report.rds"))
df <- readRDS(file.path(res_dir, "LC_counts_monthly", "check_convergence.rds"))
base_sdrep <- readRDS(file.path(res_dir, "LC_counts_monthly", "Sdreport.rds"))
base_inputs <- readRDS(file.path(res_dir, "LC_counts_monthly", "Inputs.rds"))
	### base model doesn't converge with abundance index
read_sdreport(summary(base_sdrep)[which(rownames(summary(base_sdrep))=="SPR_t"),], log=FALSE)
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, ref=0.3)$root, error=function(e) NA)

F_est <- rep(0, length(years)*12)
F_est[1:length(base_rep$F_t)] <- base_rep$F_t
F_t <- sapply(1:(length(years)), function(x) sum(F_est[(1:12)+12*(x-1)]))
SPR <- sapply(1:length(F_t), function(x) calc_ref(ages=plist_m$ages, Mat_a=base_rep$Mat_a, W_a=base_rep$W_a, M=base_rep$M, S_a=base_rep$S_a, F=F_t[x]/12))

par(mfrow=c(1,3))

plot(base_rep$R_t, type="l", lwd=4)

plot(F_t, type="l", lwd=4, ylim=c(0, max(F_t)))
abline(h=F40*12, lty=2)
abline(h=F30*12, lty=2)

plot(SPR, type="l", lwd=4, ylim=c(0, 1))
lines(lbspr_res$SPR, col="red", lwd=3)
abline(h=0.4, lty=2)
abline(h=0.3, lty=2)




