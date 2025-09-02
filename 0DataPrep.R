################################################################################
################################################################################
##### Refining methods for attributing health 
##### impacts to climate change: a heat-mortality case study in Zürich
#####     by Stuart-Smith et al. 2025 Climatic Change 
################################################################################
################################################################################

################################################################################

# 0. DATA LOADING AND PREPARATION

################################################################################
# 0. LOAD PACKAGES
library(lubridate); library(dplyr) ; library(dlnm); library(splines) ;
library(MASS) ; library(scales) ; library(ggplot2)
#library(RcppRoll)

################################################################################
# 00. LOAD DATA - DOI:
dta.attr <- read.csv("dtafinal_forpub.csv")

# TRANSFORM INTO DATE VARIABLE
dta.attr$date <- as.Date(dta.attr$date)
