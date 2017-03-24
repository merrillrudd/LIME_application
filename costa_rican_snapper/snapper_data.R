rm(list=ls())
library(LIME)
library(tidyr)

###################################
## Directories
###################################

main_dir <- "C:\\Git_Projects\\LIME_application"
Rmain_dir <- file.path(main_dir, "R_functions")
funs <- list.files(Rmain_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(Rmain_dir, funs[x])))

init_dir <- file.path(main_dir, "costa_rican_snapper")

R_dir <- file.path(init_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

data_dir <- file.path(init_dir, "data")

###############################
## life history 
###############################
## choose values (from Bystrom thesis)
linf_toUse <- 64.58
vbk_toUse <- 0.21
t0_toUse <- -0.01
lwa_toUse <- 0.0245
lwb_toUse <- 2.790

## use natural mortality tool
M_toUse <- 0.43

## add life history information for species chosen for assessment
## siganus sutor
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwb_toUse, S50=32, S95=39, selex_input="length", M50=34, maturity_input="length", SigmaR=0.3, M=M_toUse, F1=0.34, CVlen=0.1, nseasons=1, binwidth=1)


###############################
## length composition data
###############################
# data_raw <- read.csv(file.path(data_dir, "cr_snapper_database.csv"), 
# 	stringsAsFactors=FALSE)

# data_raw$Year <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x], out="year"))
# data_raw$Month <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x], out="month"))
# data_raw$TL_cm <- as.numeric(data_raw$TL_cm)
# data_raw$W_g <- as.numeric(data_raw$W_g)
# 	data_raw$W_g[which(data_raw$W_g==0)] <- NA
# data_raw$W_kg <- data_raw$W_g/1000


#############################
## subset Lutjanus guttatus
#############################

# data <- data_raw[which(data_raw$Species=="Lutjanus guttatus" & data_raw$Year>1950),]
# write.csv(data, file.path(data_dir, "cr_snapper_filtered.csv"))

data <- read.csv(file.path(data_dir, "cr_snapper_filtered.csv"), row=1)
max(data$TL_cm, na.rm=TRUE)
data$TL_cm[which(data$TL_cm==372)] <- 37.2
max(data$TL_cm, na.rm=TRUE)




#############################
## subset by gears
#############################
data_bl <- data[which(data$Gear=="Bottom Longline"),]
data_g <- data[which(data$Gear=="Gillnet"),]

data2 <- data[which(data$Gear=="Bottom Longline" | data$Gear=="Gillnet"),]

binwidth <- plist$binwidth
bins <- 1:ceiling(1.5*plist$linf)
years <- min(as.numeric(data2$Year)):max(as.numeric(data2$Year))

### raw length composition
LC <- matrix(0, nrow=length(years), ncol=length(bins))
rownames(LC) <- years
colnames(LC) <- bins
for(y in 1:length(years)){
	sub <- data2[which(data2$Year==years[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(ceiling(sub$TL_cm)==bins[x])))
	LC[y,] <- freq
}

## proportions raw length compositions
LCprop <- t(sapply(1:nrow(LC), function(x) LC[x,]/sum(LC[x,])))
rownames(LCprop) <- rownames(LC)

### weight by catch
### don't remove observations without length data2 - use them to weight by gear
gears_raw <- unique(data2$Gear)[order(unique(data2$Gear))]
gearmat_raw <- matrix(NA, nrow=length(years), ncol=length(gears_raw))
rownames(gearmat_raw) <- years
colnames(gearmat_raw) <- gears_raw
for(i in 1:length(years)){
	sub1 <- data2[which(data2$Year==years[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Gear==gears_raw[j]),]
		gearmat_raw[i,j] <- nrow(sub2)
	}
}
## in order of highest to lowest samples
gearmat <- gearmat_raw[,rev(order(colSums(gearmat_raw)))]
gear_totals <- colSums(gearmat)[rev(order(colSums(gearmat)))]
gear_prop <- gear_totals/nrow(data2)
gears <- colnames(gearmat)

## samples in each length bin by gear
lc_gear <- list()
for(g in 1:length(gears)){
	sub <- data2[which(data2$Gear==gears[g]),]
	subyears <- unique(sub$Year)[order(unique(sub$Year))]
	mat <- matrix(0, nrow=length(subyears), ncol=length(bins))
			rownames(mat) <- subyears
			colnames(mat) <- bins
	for(y in 1:length(subyears)){
		sub2 <- sub[which(sub$Year==subyears[y]),]
		freq <- sapply(1:length(bins), function(x) length(which(ceiling(sub2$TL_cm)==bins[x])))
		mat[y,] <- freq
	}
	lc_gear[[g]] <- mat
}
names(lc_gear) <- gears

## proportion of samples in each length bin by gear
lcprop_gear <- lapply(1:length(lc_gear), function(x) lc_gear[[x]]/rowSums(lc_gear[[x]]))
names(lcprop_gear) <- gears

## samples in each length bin by year
lc_year <- list()
for(y in 1:length(years)){
	sub <- data2[which(data2$Year==years[y]),]
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

## proportion of samples in each length bin by year
## all length proportions for each year/gear combo add up to 1
lcprop_year <- lapply(1:length(lc_year), function(x) lc_year[[x]]/rowSums(lc_year[[x]]))
names(lcprop_year) <- years

## weight proportions by number of samples
## proportions for each year should add up to 1
lc_weight_year <- list()
lcprop_weight_year <- list()
for(y in 1:length(years)){
	sub <- lcprop_year[[y]]
	lc_weight_year[[y]] <- sub*gearmat[which(rownames(gearmat)==years[y]),which(colnames(gearmat) %in% rownames(sub))]
	lcprop_weight_year[[y]] <- lc_weight_year[[y]]/sum(lc_weight_year[[y]])
}
names(lcprop_weight_year) <- names(lc_weight_year) <- years

lc_weight <- t(sapply(1:length(years), function(x) colSums(lc_weight_year[[x]])))
lcprop_weight <- t(sapply(1:length(years), function(x) colSums(lcprop_weight_year[[x]])))
rownames(lc_weight) <- rownames(lcprop_weight) <- years
colnames(lc_weight) <- colnames(lcprop_weight) <- bins

## unweighted
saveRDS(LC, file.path(data_dir, "CR_Snapper_LF_raw.rds"))
saveRDS(LCprop, file.path(data_dir, "CR_Snapper_LFprop_raw.rds"))

## proportions by year and gear
saveRDS(lc_weight_year, file.path(data_dir, "CR_Snapper_LF_gear_year.rds"))
saveRDS(lcprop_weight_year, file.path(data_dir, "CR_Snapper_LFprop_gear_year.rds"))

## weighted
saveRDS(lc_weight, file.path(data_dir, "CR_Snapper_LF_weighted.rds"))
saveRDS(lcprop_weight, file.path(data_dir, "CR_Snapper_LFprop_weighted.rds"))

## gears separate
saveRDS(lc_gear[["Bottom Longline"]], file.path(data_dir, "CR_Snapper_LF_bottom_longline.rds"))
saveRDS(lc_gear[["Gillnet"]], file.path(data_dir, "CR_Snapper_LF_gillnet.rds"))

## matrix with samples per  gear each year
saveRDS(gearmat, file.path(data_dir, "Samples_per_gear.rds"))

lc_list <- list("LF"=LCprop)
plot_LCdata(lc_list=lc_list, lc_years=years, gears=gears, L50=plist$ML50, S50=plist$S50)

lc_list <- lc_weight_year
plot_LCdata(lc_list=lc_list, lc_years=years, gears=gears, L50=plist$ML50, S50=plist$S50)

lc_list <- list("LF"=lc_weight)
plot_LCdata(lc_list=lc_list, lc_years=years, gears=gears, L50=plist$ML50, S50=plist$S50)

plot_LCfits(Inputs=list("LF"=LCprop), Inputs2=list("LF"=lc_weight))


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

par(mfrow=c(3,1), mar=c(0,4,0,0))
plot(mysum$Year, mysum$Nobs, ylab="Number of fish observed", xlab="",type="o", lwd=2, xaxt="n", xaxs="i", yaxs="i", xpd=NA, ylim=c(0, max(mysum$Nobs)*1.1))
plot(mysum$Year, mysum$Effort,xlab="",ylab="Number of hooks", type="o", lwd=2, xaxt="n", xaxs="i", yaxs="i", xpd=NA, ylim=c(0, max(mysum$Effort)*1.1))
plot(mysum$Year, mysum$CPUE, xlab="", ylab="CPUE (fish per hook)", type="o", lwd=2, xaxs="i", yaxs="i", xpd=NA, ylim=c(0, max(mysum$CPUE)*1.1))

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
