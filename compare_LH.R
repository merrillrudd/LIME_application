rm(list=ls())
###################################
main_dir <- "C:\\Git_Projects\\LIME_application"

R_dir <- file.path(main_dir, "R_functions")
funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

gen_fig_dir <- file.path(main_dir, "figures")
dir.create(gen_fig_dir, showWarnings = FALSE)

data_rab <- file.path(main_dir, "kenyan_reef_fish\\data")
data_snap <- file.path(main_dir, "costa_rican_snapper\\data")

plist_rab <- readRDS(file.path(data_rab, "Siganus_sutor_life_history_annual.rds"))
plist_q_rab <- readRDS(file.path(data_rab, "Siganus_sutor_life_history_quarterly.rds"))
plist_m_rab <- readRDS(file.path(data_rab, "Siganus_sutor_life_history_monthly.rds"))

plist_snap <- readRDS(file.path(data_snap, "CRSNAP_life_history_annual.rds"))
plist_q_snap <- readRDS(file.path(data_snap, "CRSNAP_life_history_quarterly.rds"))
plist_m_snap <- readRDS(file.path(data_snap, "CRSNAP_life_history_monthly.rds"))

### ages at length

png(file.path(gen_fig_dir, "Age_length_YM.png"), height=8, width=15, units="in", res=200)

ramp <- colorRamp(c("purple4", "darkorange"))
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,2,1,1))

plba <- age_length(highs=plist_snap$highs, lows=plist_snap$lows, L_a=plist_snap$L_a, CVlen=0.1)
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
matplot(t(plba), type="l", lty=1, lwd=3, col=col_vec, cex.axis=2, xaxs="i", yaxs="i", ylim=c(0, 0.5), xaxt="n")
print.letter(xy=c(0.9,0.9), "(a)", cex=2)
mtext(side=2, expression(italic("Lutjanus guttatus")), line=6, cex=2.5)
mtext(side=3, "Annual time-step", line=1, cex=2.5)

plba <- age_length(highs=plist_q_snap$highs, lows=plist_q_snap$lows, L_a=plist_q_snap$L_a, CVlen=0.1)
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
matplot(t(plba), type="l", lty=1, lwd=3, col=col_vec, cex.axis=2, xaxs="i", yaxs="i", yaxt="n", ylim=c(0, 0.5), xaxt="n")
print.letter(xy=c(0.9,0.9), "(b)", cex=2)
mtext(side=3, "Quarterly time-step", line=1, cex=2.5)

plba <- age_length(highs=plist_m_snap$highs, lows=plist_m_snap$lows, L_a=plist_m_snap$L_a, CVlen=0.1)
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
matplot(t(plba), type="l", lty=1, lwd=3, col=col_vec, yaxt="n", cex.axis=2, xaxs="i", yaxs="i", ylim=c(0, 0.5), xaxt="n")
print.letter(xy=c(0.9,0.9), "(c)", cex=2)
mtext(side=3, "Monthly time-step", line=1, cex=2.5)

plba <- age_length(highs=plist_rab$highs, lows=plist_rab$lows, L_a=plist_rab$L_a, CVlen=0.1)
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
matplot(t(plba), type="l", lty=1, lwd=3, col=col_vec, cex.axis=2, xaxs="i", yaxs="i", ylim=c(0, 0.5))
print.letter(xy=c(0.9,0.9), "(d)", cex=2)
mtext(side=2, expression(italic("Siganus sutor")), line=6, cex=2.5)

plba <- age_length(highs=plist_q_rab$highs, lows=plist_q_rab$lows, L_a=plist_q_rab$L_a, CVlen=0.1)
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
matplot(t(plba), type="l", lty=1, lwd=3, col=col_vec, yaxt="n", cex.axis=2, xaxs="i", yaxs="i", ylim=c(0, 0.5))
print.letter(xy=c(0.9,0.9), "(e)", cex=2)

plba <- age_length(highs=plist_m_rab$highs, lows=plist_m_rab$lows, L_a=plist_m_rab$L_a, CVlen=0.1)
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
matplot(t(plba), type="l", lty=1, lwd=3, col=col_vec, yaxt="n", cex.axis=2, xaxs="i", yaxs="i", ylim=c(0, 0.5))
print.letter(xy=c(0.9,0.9), "(f)", cex=2)


mtext(side=2, "Probability of being an length given age", outer=TRUE, line=3, cex=2)
mtext(side=1, "Length (cm)", outer=TRUE, line=3, cex=2)

dev.off()

