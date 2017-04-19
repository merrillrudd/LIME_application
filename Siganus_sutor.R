rm(list=ls())

## load LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

## load LBSPR package
library(LBSPR)

## load Rfishbase
library(rfishbase)
library(gplots)


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

fig_dir <- file.path(stock_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)


## data to be used for the assessment
data <- alldata[which(alldata$Name %in% name),]

## years in the dataset
years <- unique(data$Year)[order(unique(data$Year))]
years_t <- min(years):max(years)
all_months <- as.vector(sapply(1:length(years_t), function(x) paste0(1:12, "/", years_t[x])))

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
for(y in 1:length(years_t)){
	sub <- data[which(data$Year==years_t[y]),]
	if(nrow(sub)==0) next
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

plist <- readRDS(file.path(data_dir, "Siganus_sutor_life_history_annual.rds"))
plist_m <- readRDS(file.path(data_dir, "Siganus_sutor_life_history_monthly.rds"))

## M/K
plist$M/plist$vbk


######################################
## Length composition data_new
######################################

## unweighted length composition
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

## samples by gear annually
gears_raw <- unique(data_new$Fishgear)[order(unique(data_new$Fishgear))]
gearmat_raw <- matrix(NA, nrow=length(years), ncol=length(gears_raw))
rownames(gearmat_raw) <- years
colnames(gearmat_raw) <- gears_raw
for(i in 1:length(years)){
	sub1 <- data_new[which(data_new$Year==years[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Fishgear==gears_raw[j]),]
		gearmat_raw[i,j] <- nrow(sub2)
	}
}

gearmat_raw2 <- matrix(NA, nrow=length(mototal), ncol=length(gears_raw))
rownames(gearmat_raw2) <- mototal
colnames(gearmat_raw2) <- gears_raw
for(i in 1:length(mototal)){
	sub1 <- data_new[which(data_new$MoTotal==mototal[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Fishgear==gears_raw[j]),]
		gearmat_raw2[i,j] <- nrow(sub2)
	}
}

## in order of highest to lowest samples
gearmat <- gearmat_raw[,rev(order(colSums(gearmat_raw)))]
gearmat_mo <- gearmat_raw2[,rev(order(colSums(gearmat_raw2)))]
gear_totals <- colSums(gearmat)[rev(order(colSums(gearmat)))]
gear_prop <- gear_totals/nrow(data_new)
gears <- colnames(gearmat)

## samples in each length bin by gear
lc_gear <- list()
for(g in 1:length(gears)){
	sub <- data_new[which(data_new$Fishgear==gears[g]),]
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

## gears by month
lc_gear_month <- list()
for(g in 1:length(gears)){
	sub <- data_new[which(data_new$Fishgear==gears[g]),]
	subyears <- unique(sub$MoTotal)[order(unique(sub$MoTotal))]
	mat <- matrix(0, nrow=length(subyears), ncol=length(bins))
			rownames(mat) <- subyears
			colnames(mat) <- bins
	for(y in 1:length(subyears)){
		sub2 <- sub[which(sub$MoTotal==subyears[y]),]
		freq <- sapply(1:length(bins), function(x) length(which(sub2$value<=bins[x] & sub2$value>(bins[x]-bins[1]))))
		mat[y,] <- freq
	}
	if(sum(mat)!=nrow(sub)) stop("Not all lengths included in LC")
	lc_gear_month[[g]] <- mat
}
names(lc_gear_month) <- gears


## samples in each length bin by year
lc_year <- list()
for(y in 1:length(years)){
	sub <- data_new[which(data_new$Year==years[y]),]
	subgears <- as.character(unique(sub$Fishgear))
	mat <- matrix(0, nrow=length(subgears), ncol=length(bins))
		rownames(mat) <- subgears
		colnames(mat) <- bins
		for(g in 1:length(subgears)){
			sub2 <- lc_gear[[which(names(lc_gear)==subgears[g])]]
			mat[g,] <- sub2[which(rownames(sub2)==years[y]),]
		}
	lc_year[[y]] <- mat
}
names(lc_year) <- years

lcprop_year <- lapply(1:length(lc_year), function(x) lc_year[[x]]/rowSums(lc_year[[x]]))
names(lcprop_year) <- years


## samples in each length bin by month
lc_month <- list()
for(y in 1:length(mototal)){
	sub <- data_new[which(data_new$MoTotal==mototal[y]),]
	subgears <- as.character(unique(sub$Fishgear))
	mat <- matrix(0, nrow=length(subgears), ncol=length(bins))
		rownames(mat) <- subgears
		colnames(mat) <- bins
		for(g in 1:length(subgears)){
			sub2 <- lc_gear_month[[which(names(lc_gear_month)==subgears[g])]]
			mat[g,] <- sub2[which(rownames(sub2)==mototal[y]),]
		}
	lc_month[[y]] <- mat
}
names(lc_month) <- mototal

lcprop_month <- lapply(1:length(lc_month), function(x) lc_month[[x]]/rowSums(lc_month[[x]]))
names(lcprop_month) <- mototal

## weighted length comp
lc_weight_year <- list()
lcprop_weight_year <- list()
for(y in 1:length(years)){
	sub <- lcprop_year[[y]]
	lc_weight_year[[y]] <- sub*gearmat[which(rownames(gearmat)==years[y]),which(colnames(gearmat) %in% rownames(sub))]
	lcprop_weight_year[[y]] <- lc_weight_year[[y]]/sum(lc_weight_year[[y]])
}
names(lcprop_weight_year) <- names(lc_weight_year) <- years

lc_weight_month <- list()
lcprop_weight_month <- list()
for(y in 1:length(mototal)){
	sub <- lcprop_month[[y]]
	lc_weight_month[[y]] <- sub*gearmat_mo[which(rownames(gearmat_mo)==mototal[y]),which(colnames(gearmat_mo) %in% rownames(sub))]
	lcprop_weight_month[[y]] <- lc_weight_month[[y]]/sum(lc_weight_month[[y]])
}
names(lcprop_weight_month) <- names(lc_weight_month) <- mototal


lc_weight <- t(sapply(1:length(years), function(x) colSums(lc_weight_year[[x]])))
rownames(lc_weight) <- years

lc_weight_mo <- t(sapply(1:length(mototal), function(x) colSums(lc_weight_month[[x]])))
rownames(lc_weight_mo) <- mototal


## length composition
LC_bt <- matrix(0, nrow=length(years), ncol=length(bins))
	rownames(LC_bt) <- years
	colnames(LC_bt) <- bins
data_bt <- data_new[which(data_new$Fishgear==gears[1]),]
for(y in 1:length(years)){
	sub <- data_bt[which(data_bt$Year==years[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$value<=bins[x] & sub$value>(bins[x]-bins[1]))))
	LC_bt[y,] <- freq
}

mototal_bt <- unique(data_bt$MoTotal)
LCm_bt <- matrix(0, nrow=length(mototal_bt), ncol=length(bins))
	rownames(LCm_bt) <- mototal_bt
	colnames(LCm_bt) <- bins
for(y in 1:length(mototal_bt)){
	sub <- data_bt[which(data_bt$MoTotal==mototal_bt[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$value<=bins[x] & sub$value>(bins[x]-bins[1]))))
	LCm_bt[y,] <- freq
}


## plot length comp
png(file.path(fig_dir, "LC_annual_raw.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LC), lc_years=years, gears=NULL, ylim=c(0,0.25))
dev.off()

png(file.path(fig_dir, "LC_monthly_raw.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LCm), lc_years=mototal, gears=NULL, ylim=c(0,0.4))
dev.off()


png(file.path(fig_dir, "LC_annual_weighted.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=list("LF"=lc_weight), lc_years=years, gears=NULL, ylim=c(0,0.25))
dev.off()

png(file.path(fig_dir, "LC_monthly_weighted.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=list("LF"=lc_weight_mo), lc_years=mototal, gears=NULL, ylim=c(0,0.4))
dev.off()

png(file.path(fig_dir, "LC_annual_weight_byGear.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=lcprop_weight_year, lc_years=years, gears=gears, ylim=c(0,0.3), colors=c("red", "orange", "yellow", "green","blue","purple"))
dev.off()

png(file.path(fig_dir, "LC_monthly_weight_byGear.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=lcprop_weight_month, lc_years=mototal, gears=gears, ylim=c(0,0.3), colors=c("red", "orange", "yellow", "green","blue","purple"), legend=FALSE)
dev.off()

ramp <- colorRamp(c("red", "orange", "yellow", "green","blue","purple"))
col_vec <- rgb( ramp(seq(0, 1, length = length(gears))), max = 255)
plot(x=1,y=1,ylim=c(0,50),type="n",axes=F,ann=F)
gears[which(gears=="")] <- "unclassified"
legend("topright", legend=gears, col=col_vec, lwd=4)

png(file.path(fig_dir, "LC_annual_baskettrap.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LC_bt), lc_years=years, gears=NULL, ylim=c(0,0.25), colors="red")
dev.off()

png(file.path(fig_dir, "LC_monthly_baskettrap.png"), height=16, width=16, res=200, units="in")
plot_LCdata(lc_list=list("LF"=LCm_bt), lc_years=mototal_bt, gears=NULL, ylim=c(0,0.4), colors="red")
dev.off()

LCm_up <- LCm[which(rownames(LCm) %in% months_up),]
LCm_mid <- LCm[which(rownames(LCm) %in% months_mid),]

####################################
## run assessments
####################################

res_dir <- file.path(stock_dir, "results")
dir.create(res_dir, showWarnings=FALSE)

years_o <- as.numeric(rownames(LC))
years_t <- min(years_o):max(years_o)

months_o <- as.numeric(rownames(LCm))
months_t <- min(months_o):max(months_o)

input_data_y <- list("years"=years_t, "LF"=lc_weight)
input_data_m <- list("years"=months_t, "LF"=lc_weight_mo)
input_data_bt <- list("years"=months_t, "LF"=LCm_bt)
input_data_up <- list("years"=months_t, "LF"=lc_weight_mo[which(rownames(lc_weight_mo) %in% months_up),])

############### weighted length comps ##########################
## annual LBPSR
out <- file.path(res_dir, "LBSPR_LCy")
dir.create(out, showWarnings=FALSE)

	run <- run_LBSPR(modpath=out, lh=plist, species="x", input_data=input_data_y, rewrite=TRUE, simulation=FALSE)	
	lbspr_res <- readRDS(file.path(out, "LBSPR_results.rds"))


## annual LIME
out <- file.path(res_dir, "LCy")
dir.create(out, showWarnings=FALSE)

			res <- run_LIME(modpath=out, lh=plist, input_data=input_data_y, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_y$years, lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_y$LF), LBSPR=lbspr_res, ylim=c(0,0.3))
			dev.off()


### weighted length comp - counts only - all months
out <- file.path(res_dir, "LCm")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))


			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_m$years, lc_years=rownames(input_data_m$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_m$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - all months - basket trap only
out <- file.path(res_dir, "LCm_baskettrap")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_bt, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_bt$years, lc_years=rownames(input_data_bt$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_bt$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()

### weighted length comp - counts only - months with high sample size
out <- file.path(res_dir, "LCm_highN")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_up, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


## dome-shaped selex
plist_m_lowdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=41)
plist_m_highdome <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2, selex_type="dome", dome_sd=22)


par(mfrow=c(1,1), mar=c(5,5,4,4))
plot(plist_m$S_l, type="o", lwd=2, pch=19, cex.axis=2, xlab="Length (cm)", ylab="Proportion", cex.lab=2)
lines(plist_m_lowdome$S_l, type="o", lwd=2, lty=2, col="steelblue", pch=19)
lines(plist_m_highdome$S_l, type="o", lwd=2, lty=3, col="tomato", pch=19)

## low-dome
out <- file.path(res_dir, "LCm_lowdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=plist_m_lowdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()


## low-dome
out <- file.path(res_dir, "LCm_highdome")
dir.create(out, showWarnings=FALSE)

				res <- run_LIME(modpath=out, lh=plist_m, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=TRUE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=plist_m_highdome$S_l, fix_param=FALSE, theta_type=1, randomR=TRUE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))

			png(file.path(out, "output.png"), width=16, height=10, res=200, units="in")
			plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, all_years=input_data_up$years, lc_years=rownames(input_data_up$LF), LBSPR=NULL, lh=plist_m, true_years=years_t)
			dev.off()	

			png(file.path(out, "LCfits.png"), width=16, height=10, res=200, units="in")
			plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=rownames(input_data_up$LF), LBSPR=NULL, ylim=c(0,0.3))
			dev.off()



### likelihood profiles
## check fishbase
growth <- popgrowth(proper_name)
# linf_vec <- seq(min(min(growth$Loo),0.75*linf_toUse),max(max(growth$Loo),1.25*linf_toUse),length.out=10)
linf_vec <- seq(0.75*linf_toUse, 1.25*linf_toUse, length.out=10)
# vbk_vec <- seq(min(min(growth$K),0.75*vbk_toUse), max(max(growth$K),1.25*vbk_toUse), length.out=10)
vbk_vec <- seq(0.75*vbk_toUse, 1.25*vbk_toUse, length.out=10)
# M_vec <- seq(min(min(M_dist),0.75*M_toUse), max(max(M_dist),1.25*M_toUse), length.out=10)
M_vec <- seq(0.75*M_toUse, 1.25*M_toUse, length.out=10)
ML50_vec <- seq(0.75*plist$ML50, 1.25*plist$ML50, length.out=10)


## linf
out <- file.path(res_dir, "sens_linf")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(linf_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_toUse, linf=linf_vec[i], t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
linf_nll <- rep(NA, length(linf_vec))
for(i in 1:length(linf_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	linf_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "linf_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(linf_vec, linf_nll, pch=19, cex=2, xlab="Assumed Linf", ylab="NLL", ylim=c(0.8*min(linf_nll), max(linf_nll)*1.2))
points(linf_toUse, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## vbk
out <- file.path(res_dir, "sens_vbk")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(vbk_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_vec[i], linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
vbk_nll <- rep(NA, length(vbk_vec))
for(i in 1:length(vbk_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	vbk_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "vbk_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(vbk_vec, vbk_nll, pch=19, cex=2, xlab="Assumed vbk", ylab="NLL", ylim=c(0.8*min(vbk_nll), max(vbk_nll)*1.2))
points(vbk_toUse, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## M
out <- file.path(res_dir, "sens_M")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(M_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_vec[i], lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=20.2, M95=25, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
M_nll <- rep(NA, length(M_vec))
for(i in 1:length(M_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	M_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "M_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(M_vec, M_nll, pch=19, cex=2, xlab="Assumed M", ylab="NLL", ylim=c(0.8*min(M_nll), max(M_nll)*1.2))
points(M_toUse, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()

## ML50
out <- file.path(res_dir, "sens_ML50")
dir.create(out, showWarnings=FALSE)

for(i in 1:length(ML50_vec)){
	out2 <- file.path(out, i)
	dir.create(out2, showWarnings=FALSE)

	plist_m_new <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, t0=t0_toUse, M=M_toUse, lwa=lwa_med, lwb=lwb_med, S50=11, S95=15, M50=ML50_vec[i], M95=NULL, selex_input="length", maturity_input="length", nseasons=12, binwidth=1, SigmaR=0.777, SigmaF=0.2)

	res <- run_LIME(modpath=out2, lh=plist_m_new, input_data=input_data_m, est_sigma="log_sigma_R", data_avail="LC", itervec=NULL, rewrite=FALSE, simulation=FALSE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2), F_up=10, S_l_input=-1, fix_param=FALSE, theta_type=1, randomR=TRUE)
}
rep1 <- readRDS(file.path(res_dir, "LCm", "Report.rds"))
ML50_nll <- rep(NA, length(ML50_vec))
for(i in 1:length(ML50_vec)){
	rep <- readRDS(file.path(out, i, "Report.rds"))
	ML50_nll[i] <- rep$jnll
}
png(file.path(fig_dir, "ML50_NLLprof.png"), height=10, width=14, res=200, units="in")
plot(ML50_vec, ML50_nll, pch=19, cex=2, xlab="Assumed ML50", ylab="NLL", ylim=c(0.8*min(ML50_nll), max(ML50_nll)*1.2))
points(plist_m$ML50, rep1$jnll, pch=19, cex=2, col="blue")
dev.off()


#### output
## base model
out <- file.path(res_dir, "LCm")
dir.create(out, showWarnings=FALSE)

			df <- readRDS(file.path(out, "check_convergence.rds"))
			Report <- readRDS(file.path(out, "Report.rds"))
			Sdreport <- readRDS(file.path(out, "Sdreport.rds"))
			Inputs <- readRDS(file.path(out, "Inputs.rds"))
			LBSPR <- readRDS(file.path(res_dir, "LBSPR_LCy", "LBSPR_results.rds"))

	F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
  	F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)
  	Report$F_y[length(Report$F_y)]
  	Report$F_y[length(Report$F_y)]/(F40*Inputs$Data$n_s)
  	Report$F_y[length(Report$F_y)]/(F30*Inputs$Data$n_s)

  	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	SPR_ci <- read_sdreport(sd, log=FALSE)
	SPR_ci[length(Report$SPR)]
	SPR_ci[length(Report$SPR)+1]
	Report$SPR[length(Report$SPR)]
	Report$sigma_R



# dev.new()
	all_years <- months_t
	lc_years <- months_o
	true_years <- years_t

png(file.path(fig_dir, "BaseModelOutput.png"), height=10, width=18, units="in", res=200)
plot_output(all_years=months_t, lc_years=months_o, true_years=years_t, Report=Report, Inputs=Inputs, Sdreport=Sdreport, LBSPR=LBSPR, lh=plist_m)
dev.off()

png(file.path(fig_dir, "BaseModelLCfits.png"), width=16, height=14, res=200, units="in")
plot_LCfits(Inputs=Inputs$Data, Report=Report, true_lc_years=NULL, LBSPR=NULL, ylim=c(0,0.3))
dev.off()

## compare models
png(file.path(fig_dir, "CompareModels.png"), height=8, width=15, res=200, units="in")
dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "LCm_highN")
dir_list[[2]] <- file.path(res_dir, "LCm_baskettrap")
dir_list[[3]] <- file.path(res_dir, "LCm")

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}


col_list <- list()
col_list[[1]] <- col2hex("orange")
col_list[[2]] <- col2hex("orangered")
col_list[[3]] <- "#0000AA"

by <- round(length(true_years)/5)
lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
ilab <- which(true_years %in% lab)
ilab2 <- sapply(1:length(ilab), function(x){
  sub <- which(Inputs$Data$S_yrs %in% ilab[x])
  return(sub[length(sub)])
})

par(mfrow=c(1,3), mar=c(4,5,1,1), omi=c(0.2,0.2,0.2,0.2))
plot(x=1,y=1,type="n", xlim=c(1,length(true_years)), ylim=c(0,3.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Fishing mortality", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	if(Inputs$Data$n_s==1){
	  xY <- seq_along(all_years)
	  xLC <- which(all_years %in% lc_years)
	}
	if(Inputs$Data$n_s>1){
	  xY <- 1:Inputs$Data$n_y
	  xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
	}
	  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
	  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)	

	lines(x=xY, y=Report$F_y, lwd=2, col=col_list[[i]])
	points(x=xLC, y=Report$F_y[xLC], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
	  sd[,2][which(is.na(sd[,2]))] <- 0
	    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)  
	}
	if(i==1) axis(1, cex.axis=2, at=ilab, labels=lab)
	  if(i==length(dir_list)){
	  	abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2)
	  	abline(h=F30*Inputs$Data$n_s, lwd=2, lty=3)
	  }
if(i==1) legend("bottomleft", legend=c("All data", "Basket trap only", "High sampling months"), col=rev(unlist(col_list)), lty=1, lwd=4, cex=2, bty="n")
}

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,2.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Recruitment", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$R_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)	
}	

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,0.6), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="SPR", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
dev.off()

## compare selectivity models
png(file.path(fig_dir, "CompareSelex.png"), height=8, width=15, res=200, units="in")
dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "LCm_highdome")
dir_list[[2]] <- file.path(res_dir, "LCm_lowdome")
dir_list[[3]] <- file.path(res_dir, "LCm")

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}


col_list <- list()
col_list[[1]] <- col2hex("orange")
col_list[[2]] <- col2hex("orangered")
col_list[[3]] <- "#0000AA"

par(mfrow=c(2,2), mar=c(4,5,1,1), omi=c(0.2,0.2,0.2,0.2))
plot(x=1,y=1,type="n", xlim=c(1, length(plist_m$S_l)), ylim=c(0,1.1), xaxs="i", yaxs="i", cex.axis=2, xlab="Length bin (cm)", ylab="Selectivity", cex.lab=2, xaxt="n", yaxt="n")
axis(2, cex.axis=2, at=seq(0,1,by=0.25))
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=1:length(Report$S_l), y=Report$S_l, lwd=4, col=col_list[[i]])
	if(all(is.na(Sdreport))==FALSE){
	  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),]
	  sd[,2][which(is.na(sd[,2]))] <- 0
	    polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)  
	}
	if(i==1) axis(1, cex.axis=2, at=pretty(c(min(plist_m$highs),max(plist_m$highs))))

}

plot(x=1,y=1,type="n", xlim=c(1,length(true_years)), ylim=c(0,3.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Fishing mortality", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	if(Inputs$Data$n_s==1){
	  xY <- seq_along(all_years)
	  xLC <- which(all_years %in% lc_years)
	}
	if(Inputs$Data$n_s>1){
	  xY <- 1:Inputs$Data$n_y
	  xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
	}
	  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
	  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=plist_m$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)	

	lines(x=xY, y=Report$F_y, lwd=2, col=col_list[[i]])
	points(x=xLC, y=Report$F_y[xLC], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
	  sd[,2][which(is.na(sd[,2]))] <- 0
	    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)  
	}
	if(i==1) axis(1, cex.axis=2, at=ilab, labels=lab)
	  if(i==length(dir_list)){
	  	abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2)
	  	abline(h=F30*Inputs$Data$n_s, lwd=2, lty=3)
	  }
if(i==1) legend("topright", legend=c("Logistic", "Low dome", "High dome"), col=rev(unlist(col_list)), lty=1, lwd=4, cex=2, bty="n")

}

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,2.5), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="Recruitment", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$R_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)	
}	

plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,0.6), xaxs="i", yaxs="i", cex.axis=2, xlab="Year", ylab="SPR", cex.lab=2, xaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}

dev.off()


## sensitivity to biological inputs
png(file.path(fig_dir, "SensitivityBioInputs.png"), height=8, width=15, res=200, units="in")
dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_linf",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_linf",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}


col_list <- list()
col_list[[1]] <- col2hex("orange")
col_list[[2]] <- "#0000AA"
col_list[[3]] <- col2hex("orangered")

by <- round(length(true_years)/5)
lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
ilab <- which(true_years %in% lab)
ilab2 <- sapply(1:length(ilab), function(x){
  sub <- which(Inputs$Data$S_yrs %in% ilab[x])
  return(sub[length(sub)])
})

par(mfrow=c(2,2), mar=c(0,0,0,0), omi=c(1,1,1,1))
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	# if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
print.letter("Asymptotic length", xy=c(0.25,0.9), cex=2, font=2)
axis(2, seq(0,1,by=0.5), cex.axis=2)

dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_vbk",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_vbk",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	# if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
print.letter("von Bertalanffy k", xy=c(0.25,0.9), cex=2, font=2)

dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_M",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_M",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
print.letter("Natural mortality", xy=c(0.25,0.9), cex=2, font=2)
axis(2, seq(0,1,by=0.5), cex.axis=2)

dir_list <- list()
dir_list[[1]] <- file.path(res_dir, "sens_ML50",1)
dir_list[[2]] <- file.path(res_dir, "LCm")
dir_list[[3]] <- file.path(res_dir, "sens_ML50",10)

rep_list <- sdrep_list <- inp_list <- list()
for(i in 1:length(dir_list)){
	rep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Report.rds"))
	sdrep_list[[i]] <- readRDS(file.path(dir_list[[i]], "Sdreport.rds"))
	inp_list[[i]] <- readRDS(file.path(dir_list[[i]], "Inputs.rds"))
}
plot(x=1,y=1,type="n", xlim=c(1,length(all_years)), ylim=c(0,1.2), xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", cex.lab=2, xaxt="n", yaxt="n")
for(i in 1:length(dir_list)){
	Inputs <- inp_list[[i]]
	Report <- rep_list[[i]]
	Sdreport <- sdrep_list[[i]]

	lines(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col=col_list[[i]])
	points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col=col_list[[i]], pch=19, xpd=NA)
	if(all(is.na(Sdreport))==FALSE){
	sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
	sd[,2][which(is.na(sd[,2]))] <- 0
	  polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_list[[i]],"35"), border=NA)
	}
	if(i==length(dir_list)){
	  abline(h=0.4, lwd=2, lty=2)
	  abline(h=0.3, lwd=2, lty=3)		
	}
	if(i==1) axis(1, cex.axis=2, at=ilab2, labels=lab)		
}
legend("topright", legend=c("-25%", "base", "+25%"), col=unlist(col_list), lwd=4, bty="n", cex=2, title="Input value")
print.letter("Length at 50% maturity", xy=c(0.25,0.9), cex=2, font=2)
mtext(side=1, "Year", cex=2, line=4, outer=TRUE)
mtext(side=2, "SPR", cex=2, line=4, outer=TRUE)
dev.off()