# glioma-survival-analysis

# Glioma Survival Analysis

## Overview
This project analyzes survival outcomes among patients with Grade III glioma and glioblastoma (GBM) treated with either standard care or RIT. The analysis includes Kaplan–Meier curves, log-rank tests, Cox proportional hazards models, and interaction modeling to evaluate treatment effects and histologic differences.

The dataset comes from the `coin` R package and includes 37 patients with survival time, event status, treatment group, age, sex, and histology.

## Objectives
- Compare survival between RIT and Control groups
- Evaluate differences within Grade III and GBM subgroups
- Fit unadjusted and adjusted Cox models
- Test for treatment–histology interaction
- Produce reproducible survival analysis code

## Methods

### 1. Kaplan–Meier Estimation
- Base R `survfit()` used to generate survival curves
- Separate curves for Grade III and GBM
- Visual comparison of treatment groups

### 2. Log-Rank Tests
- Grade III only
- GBM only
- Stratified by histology (survival + coin packages)

### 3. Cox Proportional Hazards Models
- Unadjusted models for each histology
- Stratified Cox model
- Multivariable model adjusting for age, sex, histology
- Interaction model (treatment × histology)

### 4. Outputs
- KM plots (base R)
- Cox model summaries (`results/cox_results.txt`)
- Reproducible R script (`scripts/glioma_survival_analysis.R`)

## Repository Structure
