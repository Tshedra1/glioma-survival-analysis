############################################################
# Glioma Survival Analysis
# Author: Tamara
# Description:
#   Kaplan–Meier curves, log-rank tests, Cox models
#   (unadjusted, stratified, adjusted, interaction)
############################################################

# Load packages
library(survival)
library(coin)

#----------------------------------------------------------
# 1. Load Data
#----------------------------------------------------------

# The glioma dataset is included in the 'coin' package
data("glioma")

# Ensure factors are correctly formatted
glioma$group     <- factor(glioma$group, levels = c("Control", "RIT"))
glioma$histology <- factor(glioma$histology, levels = c("GBM", "Grade3"))
glioma$sex       <- factor(glioma$sex, levels = c("Female", "Male"))

# Subsets
g3 <- subset(glioma, histology == "Grade3")
g4 <- subset(glioma, histology == "GBM")

#----------------------------------------------------------
# 2. Kaplan–Meier Curves (Base R)
#----------------------------------------------------------
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Grade III
png("figures/km_grade3.png", width = 800, height = 600)

plot(
  survfit(Surv(time, event) ~ group, data = g3),
  main = "Grade III Glioma",
  lty = 2,
  col = c("red", "green"),
  xlab = "Survival Time (months)",
  ylab = "Survival Probability"
)
legend("topright", legend = levels(g3$group), col = c("red", "green"), lty = 2, bty = "n")

dev.off()
# GBM
png("figures/km_gbm.png", width = 800, height = 600)

plot(
  survfit(Surv(time, event) ~ group, data = g4),
  main = "GBM",
  lty = 2,
  col = c("red", "green"),
  xlab = "Survival Time (months)",
  ylab = "Survival Probability"
)
legend("topright", legend = levels(g4$group), col = c("red", "green"), lty = 2, bty = "n")

dev.off()

#overall KM plot

png("figures/km_overall.png", width = 800, height = 600)

plot(
  survfit(Surv(time, event) ~ group, data = glioma),
  main = "Overall Survival",
  lty = 2,
  col = c("red", "green"),
  xlab = "Survival Time (months)",
  ylab = "Survival Probability"
)
legend("topright", legend = levels(glioma$group), col = c("red", "green"), lty = 2, bty = "n")

dev.off()

#----------------------------------------------------------
# 3. Log-Rank Tests
#----------------------------------------------------------

# Grade III
logrank_g3 <- survdiff(Surv(time, event) ~ group, data = g3)
logrank_g3

# GBM
logrank_g4 <- survdiff(Surv(time, event) ~ group, data = g4)
logrank_g4

# Stratified log-rank (survival package)
survdiff(Surv(time, event) ~ group + strata(histology), data = glioma)

# Stratified log-rank (coin package)
logrank_test(
  Surv(time, event) ~ group | histology,
  data = glioma,
  distribution = approximate(nresample = 1e6)
)

#----------------------------------------------------------
# 4. Cox Models
#----------------------------------------------------------

# Grade III only
cox_g3 <- coxph(Surv(time, event) ~ group, data = g3)
summary(cox_g3)

# GBM only
cox_gbm <- coxph(Surv(time, event) ~ group, data = g4)
summary(cox_gbm)

# Stratified Cox model
cox_strat <- coxph(Surv(time, event) ~ group + strata(histology), data = glioma)
summary(cox_strat)

# Adjusted Cox model
cox_adj <- coxph(Surv(time, event) ~ group + age + sex + histology, data = glioma)
summary(cox_adj)

# Interaction model
cox_int <- coxph(Surv(time, event) ~ group * histology + age + sex, data = glioma)
summary(cox_int)

#----------------------------------------------------------
# 5. Save Outputs
#----------------------------------------------------------

sink("results/cox_results.txt")
cat("Cox Model Results\n\n")
print(summary(cox_adj))
print(summary(cox_int))
sink()
