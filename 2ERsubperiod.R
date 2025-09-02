################################################################################
################################################################################
##### Refining methods for attributing health 
##### impacts to climate change: a heat-mortality case study in Zürich
#####     by Stuart-Smith et al. 2025 Climatic Change 
################################################################################
################################################################################

################################################################################

# 2. ESTIMATION ER - BY SUBPERIOD

################################################################################

################################################################################


################################################################################
# ESTIMATION OF THE EXPOSURE-RESPONSE FUNCTION BY SUBPERIOD
# We estimate the overall ER by aprox. 15-year periods, 
# using the standard time-series analysis and
# distributed-lag non-linear models.

# We use the cb specifications as the ER overall.

#### SUB-PERIOD ANALYSIS

# CREATE ER BY DECADE
d1 <- c(1969, 1986, 2004)
d2 <- c(1985, 2003, 2018)

# CREATE EMPTY OBJECTS TO STORE THE PREDICTIONS AND THE MMT/MMP
cp_list <- cp_comp_list <- list()
mmt <- mmp <- nrow_sub <- vector()

# LOOP BY SUB-PERIOD
for (i in seq(length(d1))){
  
  # SELECT THE PERIOD
  dta.sub <- subset(dta.attr, year %in% d1[i]:d2[i])
  
  # DFINE VAR WITH N ROWS PER SUBPERIOD
  nrow_sub[i] <- nrow(dta.sub)
  
  # CREATE CB
  argvar <- list(fun=varfun,knots=quantile(dta.sub$tmean,varper/100,na.rm=T), 
                 Bound=range(dta.sub$tmean,na.rm=T))
  arglag <- list(knots=logknots(lag,lagnk))

  # CREATE CROSSBASIS
  cb <- crossbasis(dta.sub$tmean,lag=lag,
                   argvar=argvar,
                   arglag=arglag, group=dta.sub$year) 
  
  # MODEL FORMULA 
  formula <- deaths ~ cb + dow + ns(yday,df=dfseas):factor(year) + ns(date,df=round(length(unique(year))/dftrend/10))

  # RUN MODEL & PREDICT
  model <- glm(formula,dta.sub,family=quasipoisson,na.action="na.exclude")
  cp <- crossreduce(cb,model, by=0.1, at=c(quantile(dta.sub$tmean, 0.20, na.rm=T):max(dta.sub$tmean, na.rm=T)))
  # RE-CENTER IN MMT
  mmt[i] <- cp$predvar[which.min(cp$fit)]
  mmp[i] <- ecdf(dta.sub$tmean) (mmt[i])
  
  # STORE PREDICTIONS
  cp_list[[i]] <- crossreduce(cb,model,cen=mmt[i], by=0.1)
  cp_comp_list[[i]] <- crosspred(cb,model,cen=mmt[i], by=0.1)

}