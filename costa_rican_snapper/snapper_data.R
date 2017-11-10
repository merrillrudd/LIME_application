rm(list=ls())
library(LIME)
library(tidyr)
library(dplyr)
library(purrrlyr)

###################################
## Directories
###################################

main_dir <- "C:\\merrill\\LIME_application"
Rmain_dir <- file.path(main_dir, "R_functions")
funs <- list.files(Rmain_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(Rmain_dir, funs[x])))

init_dir <- file.path(main_dir, "costa_rican_snapper")
save_name <- "Lutjanus_guttatus"

R_dir <- file.path(init_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

data_dir <- file.path(init_dir, "data")

figs_dir <- file.path(init_dir, "figures")

###############################
## life history
###############################

plist <- readRDS(file.path(data_dir, "CRSNAP_life_history_annual.rds"))
plist_q <- readRDS(file.path(data_dir, "CRSNAP_life_history_quarterly.rds"))
plist_m <- readRDS(file.path(data_dir, "CRSNAP_life_history_monthly.rds"))

###############################
## length composition data
###############################
# data_raw <- read.csv(file.path(data_dir, "cr_snapper_database.csv"), stringsAsFactors=FALSE)

# data_raw$Year <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x], out="year"))
# data_raw$Month <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x], out="month"))
# data_raw$TL_cm <- as.numeric(data_raw$TL_cm)
# data_raw$W_g <- as.numeric(data_raw$W_g)
# 	data_raw$W_g[which(data_raw$W_g==0)] <- NA
# data_raw$W_kg <- data_raw$W_g/1000

# #############################
# ## subset Lutjanus guttatus
# #############################

# data <- data_raw[which(data_raw$Species=="Lutjanus guttatus" & data_raw$Year>1950),]
# write.csv(data, file.path(data_dir, "cr_snapper_filtered.csv"))

data <- read.csv(file.path(data_dir, "cr_snapper_filtered.csv"), row=1, stringsAsFactors=FALSE)
max(data$TL_cm, na.rm=TRUE)
data$TL_cm[which(data$TL_cm==372)] <- 37.2
max(data$TL_cm, na.rm=TRUE)


#############################
## subset by gears
#############################
data_bl <- data[which(data$Gear=="Bottom Longline"),]
data_g <- data[which(data$Gear=="Gillnet"),]

data2 <- data[which(data$Gear=="Bottom Longline" | data$Gear=="Gillnet"),]

#################
## time steps
#################

binwidth <- plist$binwidth
bins <- 1:ceiling(1.5*plist$linf)
years <- unique(data2$Year)[order(unique(data2$Year))]
years_t <- min(as.numeric(data2$Year)):max(as.numeric(data2$Year))
all_months <- as.vector(sapply(1:length(years_t), function(x) paste0(1:12, "/", years_t[x])))
all_quarters <- as.vector(sapply(1:length(years_t), function(x) paste0(1:4, "/", years_t[x])))

## months in the dataset
data2$Month <- sapply(1:nrow(data2), function(x){
	mo <- strsplit(data2$Date[x],"/")[[1]][1]
	if(mo=="?") return(NA)
	if(mo %in% paste0(0,1:9)) mo <- strsplit(mo,"0")[[1]][2]
	if(as.numeric(mo) > 13) mo <- strsplit(data2$Date[x],"/")[[1]][1]
	return(as.numeric(mo))
})

## months in the dataset
data2$Quarter <- sapply(1:nrow(data2), function(x){
	mo <- strsplit(data2$Date[x],"/")[[1]][1]
	if(mo=="?") return(NA)
	if(mo %in% paste0(0,1:9)) mo <- strsplit(mo,"0")[[1]][2]
	if(as.numeric(mo) > 13) mo <- strsplit(data2$Date[x],"/")[[1]][1]
	if(mo %in% 1:3) qu <- 1
	if(mo %in% 4:6) qu <- 2
	if(mo %in% 7:9) qu <- 3
	if(mo %in% 10:12) qu <- 4
	return(as.numeric(qu))
})

## order by month
data_new <- matrix(NA, nrow=nrow(data2), ncol=ncol(data2))
colnames(data_new) <- colnames(data2)
data_new <- as.data.frame(data_new)
ilast <- 1
for(y in 1:length(years_t)){
	sub <- data2[which(data2$Year==years_t[y]),]
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
data_new$mototal <- sapply(1:nrow(data_new), function(x) which(all_months==data_new$MoYr[x]))
moyrs <- unique(data_new$MoYr)
mototal <- unique(data_new$mototal)
moobs <- NULL
for(i in 1:length(mototal)){
	sub <- data_new[which(data_new$mototal==mototal[i]),]
	if(all(is.na(sub$TL_cm))==FALSE) moobs <- c(moobs, mototal[i])
}

data_new$QuYr <- sapply(1:nrow(data_new), function(x) paste0(data_new$Quarter[x], "/", data_new$Year[x]))
data_new$qutotal <- sapply(1:nrow(data_new), function(x) which(all_quarters==data_new$QuYr[x]))
quyrs <- unique(data_new$QuYr)
qutotal <- unique(data_new$qutotal)
quobs <- NULL
for(i in 1:length(qutotal)){
	sub <- data_new[which(data_new$qutotal==qutotal[i]),]
	if(all(is.na(sub$TL_cm))==FALSE) quobs <- c(quobs, qutotal[i])
}

## number of sampling days per year
ndays <- sapply(1:length(years), function(x){
	sub <- data_new[which(data_new$Year==years[x]),]
	return(length(unique(sub$Date)))
})
d_percentiles <- quantile(ndays, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))

ndays_mo <- sapply(1:length(moobs), function(x){
	sub <- data_new[which(data_new$mototal==moobs[x]),]
	return(length(unique(sub$Date)))
})
dm_percentiles <- quantile(ndays_mo, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))

ndays_qu <- sapply(1:length(quobs), function(x){
	sub <- data_new[which(data_new$qutotal==quobs[x]),]
	return(length(unique(sub$Date)))
})
dq_percentiles <- quantile(ndays_qu, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))


## upper 50th percentile of days sampled per year
years_up <- years[which(ndays > median(ndays))]
years_mid <- years[which(ndays > quantile(ndays, prob=0.25))]

months_up <- moobs[which(ndays_mo > median(ndays_mo))]
months_mid <- moobs[which(ndays_mo > quantile(ndays_mo, prob=0.25))]

quarters_up <- quobs[which(ndays_qu > median(ndays_qu))]
quarters_mid <- quobs[which(ndays_qu > quantile(ndays_qu, prob=0.25))]

## number of samples per year
nsamp <- sapply(1:length(years), function(x){
	sub <- data_new[which(data_new$Year==years[x]),]
	return(nrow(sub))
})
nsamp_up <- nsamp[which(years %in% years_up)]
nsamp_mid <- nsamp[which(years %in% years_mid)]

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
	freq <- sapply(1:length(bins), function(x) length(which(sub$TL_cm<=bins[x] & sub$TL_cm>(bins[x]-bins[1]))))
	LC[y,] <- freq
}
LC_up <- LC[which(years %in% years_up),]
LC_mid <- LC[which(years %in% years_mid),]

## put in chronological order
LCm <- matrix(0, nrow=length(moobs), ncol=length(bins))
	rownames(LCm) <- moobs
	colnames(LCm) <- bins
for(y in 1:length(moobs)){
	sub <- data_new[which(data_new$mototal==moobs[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$TL_cm<=bins[x] & sub$TL_cm>(bins[x]-bins[1]))))
	LCm[y,] <- freq
}
LCm_up <- LCm[which(moobs %in% months_up),]
LCm_mid <- LCm[which(moobs %in% months_mid),]

## by quarter
LCq <- matrix(0, nrow=length(quobs), ncol=length(bins))
	rownames(LCq) <- quobs
	colnames(LCq) <- bins
for(y in 1:length(quobs)){
	sub <- data_new[which(data_new$qutotal==quobs[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$TL_cm<=bins[x] & sub$TL_cm>(bins[x]-bins[1]))))
	LCq[y,] <- freq
}
LCq_up <- LCq[which(quobs %in% months_up),]
LCq_mid <- LCq[which(quobs %in% months_mid),]

LC_lbspr <- t(LC)
colnames(LC_lbspr) <- NULL
rownames(LC_lbspr) <- bins - plist$binwidth/2

write.csv(LC, file.path(data_dir, paste0(save_name, "_LC.csv")))
write.csv(LCm, file.path(data_dir, paste0(save_name, "_LCmonthly.csv")))
write.csv(LC_up, file.path(data_dir, paste0(save_name, "_LC_up.csv")))
write.csv(LC_mid, file.path(data_dir, paste0(save_name, "_LC_mid.csv")))
write.csv(LCm_up, file.path(data_dir, paste0(save_name, "_LCmonthly_up.csv")))
write.csv(LCm_mid, file.path(data_dir, paste0(save_name, "_LCmonthly_mid.csv")))
write.table(LC_lbspr, file.path(data_dir, paste0(save_name, "_LC_forLBSPR.csv")), col.names=FALSE, sep=",")
write.csv(LCq, file.path(data_dir, paste0(save_name, "_LCquarterly.csv")))
write.csv(LCq_up, file.path(data_dir, paste0(save_name, "_LCquarterly_up.csv")))
write.csv(LCq_mid, file.path(data_dir, paste0(save_name, "_LCquarterly_mid.csv")))



## check all observations are accounted for in length composition
sum(LC)==nrow(data2[which(is.na(data2$TL_cm)==FALSE),])
sum(LCm)==nrow(data2[which(is.na(data2$TL_cm)==FALSE),])
sum(LCq)==nrow(data2[which(is.na(data2$TL_cm)==FALSE),])

## length comp proportion
LCprop <- LC/rowSums(LC)
write.csv(LCprop, file.path(data_dir, paste0(save_name, "_LCprop.csv")))

LCmprop <- LCm/rowSums(LCm)
write.csv(LCmprop, file.path(data_dir, paste0(save_name, "_LCmonthlyprop.csv")))

## samples by gear annually
gears_raw <- unique(data_new$Gear)[order(unique(data_new$Gear))]
gearmat_raw <- matrix(NA, nrow=length(years), ncol=length(gears_raw))
rownames(gearmat_raw) <- years
colnames(gearmat_raw) <- gears_raw
for(i in 1:length(years)){
	sub1 <- data_new[which(data_new$Year==years[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Gear==gears_raw[j]),]
		gearmat_raw[i,j] <- nrow(sub2)
	}
}

gearmat_raw2 <- matrix(NA, nrow=length(moobs), ncol=length(gears_raw))
rownames(gearmat_raw2) <- moobs
colnames(gearmat_raw2) <- gears_raw
for(i in 1:length(moobs)){
	sub1 <- data_new[which(data_new$mototal==moobs[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Gear==gears_raw[j]),]
		gearmat_raw2[i,j] <- nrow(sub2)
	}
}

gearmat_raw3 <- matrix(NA, nrow=length(quobs), ncol=length(gears_raw))
rownames(gearmat_raw3) <- quobs
colnames(gearmat_raw3) <- gears_raw
for(i in 1:length(quobs)){
	sub1 <- data_new[which(data_new$qutotal==quobs[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Gear==gears_raw[j]),]
		gearmat_raw3[i,j] <- nrow(sub2)
	}
}

## in order of highest to lowest samples
gearmat <- gearmat_raw[,rev(order(colSums(gearmat_raw)))]
gearmat_mo <- gearmat_raw2[,rev(order(colSums(gearmat_raw2)))]
gearmat_qu <- gearmat_raw3[,rev(order(colSums(gearmat_raw3)))]
gear_totals <- colSums(gearmat)[rev(order(colSums(gearmat)))]
gear_prop <- gear_totals/nrow(data_new)
gears <- colnames(gearmat)

## samples in each length bin by gear
lc_gear <- list()
for(g in 1:length(gears)){
	sub <- data_new[which(data_new$Gear==gears[g]),]
	subyears <- unique(sub$Year)[order(unique(sub$Year))]
	mat <- matrix(0, nrow=length(subyears), ncol=length(bins))
			rownames(mat) <- subyears
			colnames(mat) <- bins
	for(y in 1:length(subyears)){
		sub2 <- sub[which(sub$Year==subyears[y]),]
		freq <- sapply(1:length(bins), function(x) length(which(sub2$TL_cm<=bins[x] & sub2$TL_cm>(bins[x]-bins[1]))))
		mat[y,] <- freq
	}
	if(sum(mat)!=nrow(sub[which(is.na(sub$TL_cm)==FALSE),])) stop("Not all lengths included in LC")
	lc_gear[[g]] <- mat
}
names(lc_gear) <- gears

## proportion of samples in each length bin by gear
lcprop_gear <- lapply(1:length(lc_gear), function(x) lc_gear[[x]]/rowSums(lc_gear[[x]]))
names(lcprop_gear) <- gears

## gears by month
lc_gear_month <- list()
for(g in 1:length(gears)){
	sub <- data_new[which(data_new$Gear==gears[g]),]
	subyears <- unique(sub$mototal)[order(unique(sub$mototal))]
	mat <- matrix(0, nrow=length(subyears), ncol=length(bins))
			rownames(mat) <- subyears
			colnames(mat) <- bins
	for(y in 1:length(subyears)){
		sub2 <- sub[which(sub$mototal==subyears[y]),]
		freq <- sapply(1:length(bins), function(x) length(which(sub2$TL_cm<=bins[x] & sub2$TL_cm>(bins[x]-bins[1]))))
		mat[y,] <- freq
	}
	if(sum(mat)!=nrow(sub[which(is.na(sub$TL_cm)==FALSE),])) stop("Not all lengths included in LC")
	lc_gear_month[[g]] <- mat
}
names(lc_gear_month) <- gears

## gears by quarter
lc_gear_quarter <- list()
for(g in 1:length(gears)){
	sub <- data_new[which(data_new$Gear==gears[g]),]
	subyears <- unique(sub$qutotal)[order(unique(sub$qutotal))]
	mat <- matrix(0, nrow=length(subyears), ncol=length(bins))
			rownames(mat) <- subyears
			colnames(mat) <- bins
	for(y in 1:length(subyears)){
		sub2 <- sub[which(sub$qutotal==subyears[y]),]
		freq <- sapply(1:length(bins), function(x) length(which(sub2$TL_cm<=bins[x] & sub2$TL_cm>(bins[x]-bins[1]))))
		mat[y,] <- freq
	}
	if(sum(mat)!=nrow(sub[which(is.na(sub$TL_cm)==FALSE),])) stop("Not all lengths included in LC")
	lc_gear_quarter[[g]] <- mat
}
names(lc_gear_quarter) <- gears



## samples in each length bin by year
lc_year <- list()
for(y in 1:length(years)){
	sub <- data_new[which(data_new$Year==years[y]),]
	subgears <- as.character(unique(sub$Gear))
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
for(y in 1:length(moobs)){
	sub <- data_new[which(data_new$mototal==moobs[y]),]
	subgears <- as.character(unique(sub$Gear))
	mat <- matrix(0, nrow=length(subgears), ncol=length(bins))
		rownames(mat) <- subgears
		colnames(mat) <- bins
		for(g in 1:length(subgears)){
			sub2 <- lc_gear_month[[which(names(lc_gear_month)==subgears[g])]]
			mat[g,] <- sub2[which(rownames(sub2)==moobs[y]),]
		}
	keep_i <- which(rowSums(mat)!=0)
	if(length(keep_i)==1){
		mat <- t(as.matrix(mat[keep_i,]))
		rownames(mat) <- subgears[keep_i]
		colnames(mat) <- bins
	}
	if(length(keep_i)>1){
		mat <- mat[keep_i,]
		rownames(mat) <- subgears[keep_i]
		colnames(mat) <- bins
	}
	lc_month[[y]] <- mat
}
names(lc_month) <- moobs

lcprop_month <- lapply(1:length(lc_month), function(x) lc_month[[x]]/rowSums(lc_month[[x]]))
names(lcprop_month) <- moobs

## samples in each length bin by quarter
lc_quarter <- list()
for(y in 1:length(quobs)){
	sub <- data_new[which(data_new$qutotal==quobs[y]),]
	subgears <- as.character(unique(sub$Gear))
	mat <- matrix(0, nrow=length(subgears), ncol=length(bins))
		rownames(mat) <- subgears
		colnames(mat) <- bins
		for(g in 1:length(subgears)){
			sub2 <- lc_gear_quarter[[which(names(lc_gear_quarter)==subgears[g])]]
			mat[g,] <- sub2[which(rownames(sub2)==quobs[y]),]
		}
	keep_i <- which(rowSums(mat)!=0)
	if(length(keep_i)==1){
		mat <- t(as.matrix(mat[keep_i,]))
		rownames(mat) <- subgears[keep_i]
		colnames(mat) <- bins
	}
	if(length(keep_i)>1){
		mat <- mat[keep_i,]
		rownames(mat) <- subgears[keep_i]
		colnames(mat) <- bins
	}
	lc_quarter[[y]] <- mat
}
names(lc_quarter) <- quobs

lcprop_quarter <- lapply(1:length(lc_quarter), function(x) lc_quarter[[x]]/rowSums(lc_quarter[[x]]))
names(lcprop_quarter) <- quobs

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
for(y in 1:length(moobs)){
	sub <- lcprop_month[[y]]
	lc_weight_month[[y]] <- sub*gearmat_mo[which(rownames(gearmat_mo)==moobs[y]),which(colnames(gearmat_mo) %in% rownames(sub))]
	lcprop_weight_month[[y]] <- lc_weight_month[[y]]/sum(lc_weight_month[[y]])
}
names(lcprop_weight_month) <- names(lc_weight_month) <- moobs

lc_weight_quarter <- list()
lcprop_weight_quarter <- list()
for(y in 1:length(quobs)){
	sub <- lcprop_quarter[[y]]
	lc_weight_quarter[[y]] <- sub*gearmat_qu[which(rownames(gearmat_qu)==quobs[y]),which(colnames(gearmat_qu) %in% rownames(sub))]
	lcprop_weight_quarter[[y]] <- lc_weight_quarter[[y]]/sum(lc_weight_quarter[[y]])
}
names(lcprop_weight_quarter) <- names(lc_weight_quarter) <- quobs


lc_weight <- t(sapply(1:length(years), function(x) colSums(lc_weight_year[[x]])))
rownames(lc_weight) <- years

lcprop_weight <- lc_weight/rowSums(lc_weight)

# for(i in 1:length(moobs)){
# 	xx <- colSums(lc_weight_month[[i]])
# 	if(any(is.na(xx))) break
# }

lc_weight_mo <- t(sapply(1:length(moobs), function(x) colSums(lc_weight_month[[x]])))
rownames(lc_weight_mo) <- moobs

lcprop_weight_mo <- lc_weight_mo/rowSums(lc_weight_mo)

lc_weight_qu <- t(sapply(1:length(quobs), function(x) colSums(lc_weight_quarter[[x]])))
rownames(lc_weight_qu) <- quobs

lcprop_weight_qu <- lc_weight_qu/rowSums(lc_weight_qu)


## length composition
LC_g1 <- matrix(0, nrow=length(years), ncol=length(bins))
	rownames(LC_g1) <- years
	colnames(LC_g1) <- bins
data_g1 <- data_new[which(data_new$Gear==gears[1]),]
for(y in 1:length(years)){
	sub <- data_g1[which(data_g1$Year==years[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$TL_cm<=bins[x] & sub$TL_cm>(bins[x]-bins[1]))))
	LC_g1[y,] <- freq
}

moobs_g1 <- unique(data_g1$mototal)
LCm_g1 <- matrix(0, nrow=length(moobs_g1), ncol=length(bins))
	rownames(LCm_g1) <- moobs_g1
	colnames(LCm_g1) <- bins
for(y in 1:length(moobs_g1)){
	sub <- data_g1[which(data_g1$mototal==moobs_g1[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$TL_cm<=bins[x] & sub$TL_cm>(bins[x]-bins[1]))))
	LCm_g1[y,] <- freq
}

quobs_g1 <- unique(data_g1$qutotal)
LCq_g1 <- matrix(0, nrow=length(quobs_g1), ncol=length(bins))
	rownames(LCq_g1) <- quobs_g1
	colnames(LCq_g1) <- bins
for(y in 1:length(quobs_g1)){
	sub <- data_g1[which(data_g1$qutotal==quobs_g1[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(sub$TL_cm<=bins[x] & sub$TL_cm>(bins[x]-bins[1]))))
	LCq_g1[y,] <- freq
}


## unweighted
saveRDS(LC, file.path(data_dir, "Lutjanus_guttatus_LF_raw.rds"))
saveRDS(LCprop, file.path(data_dir, "Lutjanus_guttatus_LFprop_raw.rds"))

## proportions by year and gear
saveRDS(lc_weight_year, file.path(data_dir, "Lutjanus_guttatus_LF_gear_year.rds"))
saveRDS(lcprop_weight_year, file.path(data_dir, "Lutjanus_guttatus_LFprop_gear_year.rds"))

saveRDS(lc_weight_month, file.path(data_dir, "Lutjanus_guttatus_LF_gear_monthly.rds"))
saveRDS(lcprop_weight_month, file.path(data_dir, "Lutjanus_guttatus_LFprop_gear_monthly.rds"))

saveRDS(lc_weight_quarter, file.path(data_dir, "Lutjanus_guttatus_LF_gear_quarterly.rds"))
saveRDS(lcprop_weight_quarter, file.path(data_dir, "Lutjanus_guttatus_LFprop_gear_quarterly.rds"))


## weighted
saveRDS(lc_weight, file.path(data_dir, "Lutjanus_guttatus_LF_weighted.rds"))
saveRDS(lcprop_weight, file.path(data_dir, "Lutjanus_guttatus_LFprop_weighted.rds"))

saveRDS(lc_weight_mo, file.path(data_dir, "Lutjanus_guttatus_LF_weighted_monthly.rds"))
saveRDS(lcprop_weight_mo, file.path(data_dir, "Lutjanus_guttatus_LFprop_weighted_monthly.rds"))

saveRDS(lc_weight_qu, file.path(data_dir, "Lutjanus_guttatus_LF_weighted_quarterly.rds"))
saveRDS(lcprop_weight_qu, file.path(data_dir, "Lutjanus_guttatus_LFprop_weighted_quarterly.rds"))

saveRDS(lc_weight[which(years_t %in% years_up),], file.path(data_dir, "Lutjanus_guttatus_LF_weighted_Nup.rds"))
saveRDS(lc_weight_mo[which(moobs %in% months_up),], file.path(data_dir, "Lutjanus_guttatus_LFmonthly_weighted_Nup.rds"))
saveRDS(lc_weight_qu[which(quobs %in% quarters_up),], file.path(data_dir, "Lutjanus_guttatus_LFquarterly_weighted_Nup.rds"))



## gears separate
saveRDS(lc_gear[["Bottom Longline"]], file.path(data_dir, "Lutjanus_guttatus_LF_bottomlongline.rds"))
saveRDS(lc_gear_quarter[["Bottom Longline"]], file.path(data_dir, "Lutjanus_guttatus_LF_bottomlongline_quarterly.rds"))

saveRDS(lc_gear[["Gillnet"]], file.path(data_dir, "Lutjanus_guttatus_LF_gillnet.rds"))
saveRDS(lc_gear_quarter[["Gillnet"]], file.path(data_dir, "Lutjanus_guttatus_LF_gillnet_quarterly.rds"))

## matrix with samples per  gear each year
saveRDS(gearmat, file.path(data_dir, "Samples_per_gear.rds"))

lc_list <- list("LF"=LCprop)
plot_LCfits(Inputs=lc_list)

lc_list <- list("LF"=lc_weight)
plot_LCdata(lc_list=lc_list, lc_years=years, gears=gears)

plot_LCfits(Inputs=list("LF"=LCprop), Inputs2=list("LF"=lc_weight))
plot_LCfits(Inputs=list("LF"=lc_gear[["Bottom Longline"]]/rowSums(lc_gear[["Bottom Longline"]])), Inputs2=list("LF"=lc_gear[["Gillnet"]]/rowSums(lc_gear[["Gillnet"]])))


#############################
## CPUE
#############################

idata <- data
bldata <- idata[which(idata$Gear=="Bottom Longline"),]
bldata$Observation_type[which(bldata$Observation_type=="Dockside ")] <- "Dockside"
data_avail <- bldata[which(is.na(bldata$Line_n_hooks)==FALSE),]

mysum <- aggregate(data_avail[,"Year"],by=list(Year=data_avail$Year),length)
names(mysum) <- c("Year", "Nobs")

years <- mysum$Year

effort <- rep(0, length(years))
nodata <- NULL
for(y in 1:length(years)){
	sub <- data_avail[which(data_avail$Year==years[y]),]
	dates <- unique(sub$Date)
	for(d in 1:length(dates)){
		sub2 <- sub[which(sub$Date==dates[d]),]
		passes <- unique(sub2$Pass)
		for(p in 1:length(passes)){
			if(is.na(passes[p])) next
			sub3 <- sub2[which(sub2$Pass==passes[p]),]
			if(length(unique(sub3$Bait))>1) stop("multiple baits per pass")
			if(length(unique(sub3$Observation_type))>1) stop("Multiple observation types in one pass")
			if(length(unique(sub3$Line_n_hooks))>1) stop("More than one number of hooks specified per pass")
			effort[y] <- effort[y] + sub3$Line_n_hooks[1]
		if(is.na(effort[y])) stop("Effort NA")
		}
	}
}

mysum$Effort <- effort
mysum$CPUE <- mysum$Nobs/mysum$Effort

png(file.path(figs_dir, "nominal_CPUE.png"), height=8, width=12, units="in", res=200)
par(mfrow=c(3,1), mar=c(0,4,0,0), omi=c(1,1,1,1))
plot(mysum$Year, mysum$Nobs, ylab="Number of fish observed", xlab="",type="o", lwd=2, xaxt="n", xaxs="i", yaxs="i", xpd=NA, ylim=c(0, max(mysum$Nobs)*1.1))
plot(mysum$Year, mysum$Effort,xlab="",ylab="Number of hooks", type="o", lwd=2, xaxt="n", xaxs="i", yaxs="i", xpd=NA, ylim=c(0, max(mysum$Effort)*1.1))
plot(mysum$Year, mysum$CPUE, xlab="", ylab="CPUE (fish per hook)", type="o", lwd=2, xaxs="i", yaxs="i", xpd=NA, ylim=c(0, max(mysum$CPUE)*1.1))
mtext("Year", side=1, line=3)
dev.off()

cpuemat <- mysum[which(mysum$Effort!=0),]
cpuemat$logCPUE <- log(cpuemat$CPUE)

## save nominal CPUE
I_t <- cpuemat$CPUE
names(I_t) <- cpuemat$Year
saveRDS(I_t, file.path(data_dir, "Nominal_CPUE.rds"))


## CPUE standardization
par(mfrow=c(1,2))
hist(cpuemat$CPUE, main="CPUE", xlab="CPUE")
hist(cpuemat$logCPUE, main="log(CPUE)", xlab="log(CPUE)")

par(mfrow=c(1,2))
qqnorm(cpuemat$logCPUE)
qqline(cpuemat$logCPUE, col="red")
boxplot(logCPUE ~ as.factor(Year), data=cpuemat, na.rm=TRUE, main="Boxplot nominal logCPUE by year", ylab="log(CPUE)", col="light grey", boxwex=0.65)

par(mfrow=c(1,2))
qqnorm(cpuemat$CPUE)
qqline(cpuemat$CPUE, col="red")
boxplot(CPUE ~ as.factor(Year), data=cpuemat, na.rm=TRUE, main="Boxplot nominal CPUE by year", ylab="CPUE", col="light grey", boxwex=0.65)

## use logCPUE
par(mfrow=c(1,1))
myfactors <- c("Year", "Observation_type", "Boat", "Captain", "Observer", "Location", "Fishing_zone", "Bait", "Line_hook_type")
tmp <- bldata[,which(names(bldata) %in% myfactors)]
