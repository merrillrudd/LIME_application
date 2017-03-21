rm(list=ls())
library(LIME)

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


# #############################
# ## subset Lutjanus guttatus
# #############################

# data <- data_raw[which(data_raw$Species=="Lutjanus guttatus" & data_raw$Year>1950),]
# write.csv(data, file.path(data_dir, "cr_snapper_filtered.csv"))

data <- read.csv(file.path(data_dir, "cr_snapper_filtered.csv"))
max(data$TL_cm, na.rm=TRUE)
data$TL_cm[which(data$TL_cm==372)] <- 37.2
max(data$TL_cm, na.rm=TRUE)

#############################
## subset by gears
#############################
data_bl <- data[which(data$Gear=="Bottom Longline"),]
data_g <- data[which(data$Gear=="Gillnet"),]

binwidth <- plist$binwidth
bins <- 1:ceiling(1.5*plist$linf)
years <- min(as.numeric(data$Year)):max(as.numeric(data$Year))

### raw length composition
LC <- matrix(0, nrow=length(years), ncol=length(bins))
rownames(LC) <- years
colnames(LC) <- bins
for(y in 1:length(years)){
	sub <- data[which(data$Year==years[y]),]
	freq <- sapply(1:length(bins), function(x) length(which(ceiling(sub$TL_cm)==bins[x])))
	LC[y,] <- freq
}

## proportions raw length compositions
LCprop <- t(sapply(1:nrow(LC), function(x) LC[x,]/sum(LC[x,])))
rownames(LCprop) <- rownames(LC)

### weight by catch
### don't remove observations without length data - use them to weight by gear
gears_raw <- unique(data$Gear)[order(unique(data$Gear))]
gearmat_raw <- matrix(NA, nrow=length(years), ncol=length(gears_raw))
rownames(gearmat_raw) <- years
colnames(gearmat_raw) <- gears_raw
for(i in 1:length(years)){
	sub1 <- data[which(data$Year==years[i]),]
	for(j in 1:length(gears_raw)){
		sub2 <- sub1[which(sub1$Gear==gears_raw[j]),]
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
	sub <- data[which(data$Gear==gears[g]),]
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
	sub <- data[which(data$Year==years[y]),]
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
