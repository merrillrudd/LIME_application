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

LCm_up <- LCm[which(rownames(LCm) %in% months_up),]
LCm_mid <- LCm[which(rownames(LCm) %in% months_mid),]


saveRDS(LC, file.path(data_dir, "Siganus_sutor_LCraw.rds"))
saveRDS(LCm, file.path(data_dir, "Siganus_sutor_LCraw_monthly.rds"))
saveRDS(lc_weight, file.path(data_dir, "Siganus_sutor_LC_weighted.rds"))
saveRDS(lc_weight_mo, file.path(data_dir, "Siganus_sutor_LCmonthly_weighted.rds"))
saveRDS(LCm_bt, file.path(data_dir, "Siganus_sutor_LCmonthly_baskettrap.rds"))
saveRDS(lc_weight_mo[which(rownames(lc_weight_mo) %in% months_up),], file.path(data_dir, "Siganus_sutor_LCmonthly_weighted_Nup.rds"))








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

