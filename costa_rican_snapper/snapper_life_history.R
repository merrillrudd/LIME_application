rm(list=ls())
library(LIME)
library(FishLife)
library(msm)

main_dir <- "C:\\merrill\\LIME_application"
init_dir <- file.path(main_dir, "costa_rican_snapper")
data_dir <- file.path(init_dir, "data")

#######################################
## life history parameter distributions
#######################################

## initial values (from Bystrom thesis)
linf_toUse <- 64.58
vbk_toUse <- 0.21
t0_toUse <- -0.01
lwa_toUse <- 0.0245
lwb_toUse <- 2.790
M_init <- 0.43
ML50_toUse <- 34

## FishLife
Predict <- Plot_taxa(Search_species(Genus = "Lutjanus", Species = "guttatus")$match_taxonomy, mfrow=c(2,2))
pred <- exp(Predict[[1]]$Mean_pred)
var <- Predict[[1]]$Cov_pred


## check fishbase and natural mortality toolfor life history distributions
proper_name <- "Lutjanus guttatus"
genus <- strsplit(proper_name, " ")[[1]][1]

growth <- popgrowth(proper_name)
growth2 <- popgrowth(species_list(Genus=genus))

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


## add life history information for species chosen for assessment
## siganus sutor
plist <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwa_toUse, S50=32, S95=39, selex_input="length", M50=ML50_toUse, maturity_input="length", SigmaR=0.7, SigmaF=0.2, M=M_toUse, F1=0.34, CVlen=0.1, nseasons=1, binwidth=1, AgeMax=22)
plist_q <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwa_toUse, S50=32, S95=39, selex_input="length", M50=ML50_toUse, maturity_input="length", SigmaR=0.7, SigmaF=0.2, M=M_toUse, F1=0.34, CVlen=0.1, nseasons=4, binwidth=1, AgeMax=22)
plist_m <- create_lh_list(vbk=vbk_toUse, linf=linf_toUse, lwa=lwa_toUse, lwb=lwa_toUse, S50=32, S95=39, selex_input="length", M50=ML50_toUse, maturity_input="length", SigmaR=0.7, SigmaF=0.2, M=M_toUse, F1=0.34, CVlen=0.1, nseasons=12, binwidth=1, AgeMax=22)


saveRDS(plist, file.path(data_dir, "CRSNAP_life_history_annual.rds"))
saveRDS(plist_m, file.path(data_dir, "CRSNAP_life_history_monthly.rds"))
saveRDS(plist_q, file.path(data_dir, "CRSNAP_life_history_quarterly.rds"))
