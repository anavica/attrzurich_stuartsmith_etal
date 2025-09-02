################################################################################
################################################################################
##### Refining methods for attributing health 
##### impacts to climate change: a heat-mortality case study in Zürich
#####     by Stuart-Smith et al. 2025 Climatic Change 
################################################################################
################################################################################

################################################################################

# 1. ESTIMATION ER - OVERALL PERIOD

################################################################################

################################################################################
# ESTIMATION OF THE EXPOSURE-RESPONSE FUNCTION
# We estimate the overall ER (the average risk across the whole study period)
# using the standard time-series analysis and
# distributed-lag non-linear models.

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
varper <- c(50,90)

# SPECIFICATION OF THE LAG FUNCTION 
lag <- 7
lagnk <- 2

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 4

# DEGREE OF FREEDOM FOR TREND
dftrend <- 1

# MODEL FORMULA 
formula <- deaths ~ cb + dow + ns(yday,df=dfseas):factor(year) + ns(date,df=round(length(unique(year))/dftrend/10))

#### WHOLE PERIOD ANALYSIS
# We run a full-period analysis (1969 - 2017)

# CREATE CB
argvar <- list(fun=varfun,knots=quantile(dta.attr$tmean,varper/100,na.rm=T), 
                 Bound=range(dta.attr$tmean,na.rm=T))
arglag <- list(knots=logknots(lag,lagnk))

cb <- crossbasis(dta.attr$tmean,lag=lag,
                 argvar=argvar,
                 arglag=arglag, group=dta.attr$year) 

# RUN MODEL & PREDICT
model <- glm(formula,dta.attr,family=quasipoisson,na.action="na.exclude")
cp <- crossreduce(cb,model,cen=mean(dta.attr$tmean,na.rm=T), by=0.1)

# RE-CENTER IN MMT & STORE PREDICTIONS
mmt_overall <- cp$predvar[which.min(cp$fit)]
cp_overall <- crossreduce(cb,model,cen=mmt_overall, by=0.1)
