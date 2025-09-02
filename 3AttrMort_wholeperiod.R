
################################################################################
################################################################################
##### Refining methods for attributing health 
##### impacts to climate change: a heat-mortality case study in Zürich
#####     by Stuart-Smith et al. 2025 Climatic Change 
################################################################################
################################################################################

################################################################################

# 3. ATTRIBUTABLE HEAT-MORTALITY - WHOLE PERIOD

################################################################################


################################################################################
# QUANTIFICATION OF ATTRIBUTABLE MORTALITY UNDER EACH SCENARIO
# We calculate the heat-attributable mortality using the different temperatures 
# series (daily): 1) factual scenario (Tmean) 2) counterfactual scenario (tmean 
# output from different GCMs, and observational temperatures).

# For the estimation of the attributable mortality, we define 3 scenarios:
# 1) Accounting for adaptation: use period specific ER functions for both scenarios
# 2) No adaptation: use the ER of the 2nd period for both scenarios in 3rd period
# 3) Average risk: Use the overall ER estimate using the whole period


# CREATE OBJECT WITH THE NAMES OF THE SERIES TO SELECT IN THE LOOP
temp_sim <- names(dta.attr)[c(5, 9:22)]

# SET NUMBER OF SIMULATIONS
nsim <- 1000

## SCENARIO 1: ADAPTATION
# CREATE EMPTY ARRAY TO STORE RESULTS
ancitysim_tot_1 <- array(NA,dim=list(length(temp_sim),length(d1),2,nsim+1),
                         dimnames=list(temp_sim,c(paste0("dec",seq(length(d1)))),
                                       c("abs","dif"),
                                       c("est",paste0("sim",seq(nsim)))))

anannual_tot_1 <- array(NA,dim=list(length(unique(year(dta.attr$date))),
                                     length(temp_sim),2),
                         dimnames=list(unique(year(dta.attr$date)),temp_sim,
                                       c("abs","dif")))

tmeanannual_1 <- array(NA,dim=list(length(unique(year(dta.attr$date))),
                                     length(temp_sim),2),
                         dimnames=list(unique(year(dta.attr$date)),temp_sim,
                                       c("abs","dif")))

# CREATE VECTOR TO STORE RESULTS
deaths_tot <- vector()

# LOOP BY TEMPERARTURE SERIES ACROSS SUB-PERIODS
for (g in seq(length(temp_sim))){
  
  for (i in seq(length(d1))){
    
    # SELECT THE PERIOD
    dta.sub <- subset(dta.attr, year %in% d1[i]:d2[i])
    
    # EXTRACT PARAMETERS
    coef <- cp_list[[i]]$coefficients
    vcov <- cp_list[[i]]$vcov
    
    # DEFINE ARGVAR (AS IN ESTIMATION), CENTERING AND COEF-VCOV
    argvar <- list(fun=varfun,knots=quantile(dta.sub$tmean_obs,varper/100,
            na.rm=T),Bound=range(dta.sub$tmean_obs,na.rm=T))

    # DERIVE THE CENTERED BASIS
    bvar <- do.call(onebasis,c(list(x=dta.sub[,temp_sim[g]]),argvar))
    cenvec <- do.call(onebasis,c(list(x=mmt[i]),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
      
    # INDICATOR FOR HEAT DAYS
    indheat <- dta.sub[,temp_sim[g]]>mmt[i]
       
    # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS (forward perspective, 
    # we use the moving average of the deaths observed the next 7 days)
    an <- (1-exp(-bvarcen%*%coef))*dta.sub$mvdeaths
    
    # STORE NUMBER OF ANNUAL DEATHS
    seqselyear <- which(unique(year(dta.attr$date)) %in% c(d1[i]:d2[i]))
    anannual_tot_1[seqselyear,g,1] <- tapply(an[indheat], year(dta.sub$date[indheat]), sum, na.rm=T)
    
    # COMPUTE AVERAGE ANNUAL TEMPERATURE
    tmeanannual_1[seqselyear,g,1] <- tapply(dta.sub[,temp_sim[g]],year(dta.sub$date),mean, na.rm=T)
    
    # SUM 
    ancitysim_tot_1[g,i,1,1] <- sum(an[indheat], na.rm=T)

    # COMPUTE SUM 
    deaths_tot[i] <- sum(dta.sub$mvdeaths[!is.na(dta.sub[,temp_sim[g]]) & !is.na(an)])

    # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975)
    coefsim <- mvrnorm(nsim,coef,vcov)
        
    # LOOP ACROSS ITERATIONS
    for(s in seq(nsim)) {
          
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*dta.sub$mvdeaths
                  
      # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
      # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
      ancitysim_tot_1[g,i,1,s+1] <- sum(an[indheat], na.rm=T)
    }
  }  
}  

## SCENARIO 2: NO ADAPTATION
# CREATE EMPTY ARRAY TO STORE RESULTS
ancitysim_tot_2 <- array(NA,dim=list(length(temp_sim),length(d1),2,nsim+1),
                         dimnames=list(temp_sim,c(paste0("dec",seq(length(d1)))),
                                       c("abs","dif"),
                                       c("est",paste0("sim",seq(nsim)))))

anannual_tot_2 <- array(NA,dim=list(length(unique(year(dta.attr$date))),
                                     length(temp_sim),2),
                         dimnames=list(unique(year(dta.attr$date)),temp_sim,
                                       c("abs","dif")))

tmeanannual_2 <- array(NA,dim=list(length(unique(year(dta.attr$date))),
                                     length(temp_sim),2),
                         dimnames=list(unique(year(dta.attr$date)),temp_sim,
                                       c("abs","dif")))

# CREATE VECTOR TO STORE RESULTS
deaths_tot <- vector()

# LOOP BY TEMPERATURE SERIES AND ACROSS SUBPERIODS
for (g in seq(length(temp_sim))){
  
  for (i in seq(length(d1))){
    
    # SELECT THE PERIOD
    dta.sub <- subset(dta.attr, year %in% d1[i]:d2[i])
    
    # EXTRACT PARAMETERS - use ER of the 2nd decade
    if(i==3){ coef <- cp_list[[2]]$coefficients
    vcov <- cp_list[[2]]$vcov} else {
      
    coef <- cp_list[[i]]$coefficients
    vcov <- cp_list[[i]]$vcov
    }

    # DEFINE ARGVAR (AS IN ESTIMATION), CENTERING AND COEF-VCOV
   if(i==3){ dta.argvar <- subset(dta.attr, year %in% d1[2]:d2[2])
   } else {dta.argvar <- subset(dta.attr, year %in% d1[i]:d2[i])}
    
    argvar <- list(fun=varfun,knots=quantile(dta.argvar$tmean_obs,varper/100,
            na.rm=T),Bound=range(dta.argvar$tmean_obs,na.rm=T))

    # DERIVE THE CENTERED BASIS
    if(i==3){ bvar <- do.call(onebasis,c(list(x=dta.sub[,temp_sim[g]]),argvar))
    cenvec <- do.call(onebasis,c(list(x=mmt[2]),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
    
    # INDICATOR FOR HEAT DAYS
    indheat <- dta.sub[,temp_sim[g]]>mmt[2]
    }else{
    bvar <- do.call(onebasis,c(list(x=dta.sub[,temp_sim[g]]),argvar))
    cenvec <- do.call(onebasis,c(list(x=mmt[i]),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
    
    # INDICATOR FOR HEAT DAYS
    indheat <- dta.sub[,temp_sim[g]]>mmt[i]
    }
      
    # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
    an <- (1-exp(-bvarcen%*%coef))*dta.sub$mvdeaths
  
    # STORE NUMBER OF ANNUAL DEATHS
    seqselyear <- which(unique(year(dta.attr$date)) %in% c(d1[i]:d2[i]))
    anannual_tot_2[seqselyear,g,1] <- tapply(an[indheat], year(dta.sub$date[indheat]), sum, na.rm=T)
    
    # COMPUTE AVERAGE ANNUAL TEMPERATURE
    tmeanannual_2[seqselyear,g,1] <- tapply(dta.sub[,temp_sim[g]],year(dta.sub$date),mean, na.rm=T)
    
    # SUM 
    ancitysim_tot_2[g,i,1,1] <- sum(an[indheat], na.rm=T)

    # COMPUTE SUM 
    deaths_tot[i] <- sum(dta.sub$mvdeaths[!is.na(dta.sub[,temp_sim[g]]) & !is.na(an)])

    # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975)
    coefsim <- mvrnorm(nsim,coef,vcov)
        
    # LOOP ACROSS ITERATIONS
    for(s in seq(nsim)) {
          
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*dta.sub$mvdeaths
                  
      # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
      # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
      ancitysim_tot_2[g,i,1,s+1] <- sum(an[indheat], na.rm=T)
    }
  }  
} 

## SCENARIO 3: AVERAGE RISK
# CREATE EMPTY ARRAY TO STORE RESULTS
ancitysim_tot_3 <- array(NA,dim=list(length(temp_sim),length(d1),2,nsim+1),
                         dimnames=list(temp_sim,c(paste0("dec",seq(length(d1)))),
                                       c("abs","dif"),
                                       c("est",paste0("sim",seq(nsim)))))

anannual_tot_3 <- array(NA,dim=list(length(unique(year(dta.attr$date))),
                                     length(temp_sim),2),
                         dimnames=list(unique(year(dta.attr$date)),temp_sim,
                                       c("abs","dif")))

tmeanannual_3 <- array(NA,dim=list(length(unique(year(dta.attr$date))),
                                     length(temp_sim),2),
                         dimnames=list(unique(year(dta.attr$date)),temp_sim,
                                       c("abs","dif")))

# CREATE VECTOR TO STORE RESULTS
deaths_tot <- vector()

# RUN CODE
for (g in seq(length(temp_sim))){
  
  for (i in seq(length(d1))){
    
    # SELECT THE PERIOD
    dta.sub <- subset(dta.attr, year %in% d1[i]:d2[i])
    
    # EXTRACT PARAMETERS
    coef <- cp_overall$coefficients
    vcov <- cp_overall$vcov
    
    # DEFINE ARGVAR (AS IN ESTIMATION), CENTERING AND COEF-VCOV
    argvar <- list(fun=varfun,knots=quantile(dta.attr$tmean_obs,varper/100,
            na.rm=T),Bound=range(dta.attr$tmean_obs,na.rm=T))

    # DERIVE THE CENTERED BASIS
    bvar <- do.call(onebasis,c(list(x=dta.sub[,temp_sim[g]]),argvar))
    cenvec <- do.call(onebasis,c(list(x=mmt_overall),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
      
    # INDICATOR FOR HEAT DAYS
    indheat <- dta.sub[,temp_sim[g]]>mmt_overall
      
    # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
    an <- (1-exp(-bvarcen%*%coef))*dta.sub$mvdeaths
    
    # STORE NUMBER OF ANNUAL DEATHS
    seqselyear <- which(unique(year(dta.attr$date)) %in% c(d1[i]:d2[i]))
    anannual_tot_3[seqselyear,g,1] <- tapply(an[indheat], year(dta.sub$date[indheat]), sum, na.rm=T)
    
    # COMPUTE AVERAGE ANNUAL TEMPERATURE
    tmeanannual_3[seqselyear,g,1] <- tapply(dta.sub[,temp_sim[g]],year(dta.sub$date),mean, na.rm=T)
  
    # SUM 
    ancitysim_tot_3[g,i,1,1] <- sum(an[indheat], na.rm=T)

    # COMPUTE SUM 
    deaths_tot[i] <- sum(dta.sub$mvdeaths[!is.na(dta.sub[,temp_sim[g]]) & !is.na(an)])

    # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975)
    coefsim <- mvrnorm(nsim,coef,vcov)
        
    # LOOP ACROSS ITERATIONS
    for(s in seq(nsim)) {
          
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*dta.sub$mvdeaths
                  
      # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
      # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
      ancitysim_tot_3[g,i,1,s+1] <- sum(an[indheat], na.rm=T)
    }
  }  
} 


# EXTRACT THE SUMMARIES

# ESTIMATE DIFF BETWEEN MODELLED (COUNTERFACT) AND THE OBSERVED SERIES (FACT)
ancitysim_tot_1[,,2,] <- ancitysim_tot_1[rep(1,length(temp_sim)),,1,] - ancitysim_tot_1[,,1,] 
ancitysim_tot_2[,,2,] <- ancitysim_tot_2[rep(1,length(temp_sim)),,1,] - ancitysim_tot_2[,,1,] 
ancitysim_tot_3[,,2,] <- ancitysim_tot_3[rep(1,length(temp_sim)),,1,] - ancitysim_tot_3[,,1,] 

anannual_tot_1[,,2] <- anannual_tot_1[,rep(1,length(temp_sim)),1] - anannual_tot_1[,,1]
anannual_tot_2[,,2] <- anannual_tot_2[,rep(1,length(temp_sim)),1] - anannual_tot_2[,,1]
anannual_tot_3[,,2] <- anannual_tot_3[,rep(1,length(temp_sim)),1] - anannual_tot_3[,,1]

tmeanannual_1[,,2] <- tmeanannual_1[,rep(1,length(temp_sim)),1] - tmeanannual_1[,,1]
tmeanannual_2[,,2] <- tmeanannual_2[,rep(1,length(temp_sim)),1] - tmeanannual_2[,,1]
tmeanannual_3[,,2] <- tmeanannual_3[,rep(1,length(temp_sim)),1] - tmeanannual_3[,,1]

annualdeaths <- tapply(dta.attr$deaths, year(dta.attr$date), sum,na.rm=T)


# SERIES-SPECIFIC RESULTS
# CREATE EMPTY ARRAY TO STORE SUMMARY RESULTS
an_tot1 <- af_tot1 <- an_tot2 <- af_tot2 <- an_tot3 <- 
  af_tot3 <- array(0,dim=c(length(temp_sim),length(d1),3,2)
    ,dimnames=list(c(temp_sim),paste0("dec", seq(length(d1))),
    c("est","ci.l","ci.u"),
    c("abs","dif")))

an_tot_over1 <- af_tot_over1 <- an_tot_over2 <- af_tot_over2 <- an_tot_over3 <- 
  af_tot_over3 <- array(0,dim=c(length(temp_sim),3,2)
    ,dimnames=list(c(temp_sim),c("est","ci.l","ci.u"),
    c("abs","dif")))


# TOTAL
an_tot1[,,1,] <- ancitysim_tot_1[,,,1]
an_tot1[,,2,] <- apply(ancitysim_tot_1[,,,-1],1:3,quantile,0.025, na.rm=T)
an_tot1[,,3,] <- apply(ancitysim_tot_1[,,,-1],1:3,quantile,0.975, na.rm=T)

an_tot2[,,1,] <- ancitysim_tot_2[,,,1]
an_tot2[,,2,] <- apply(ancitysim_tot_2[,,,-1],1:3,quantile,0.025, na.rm=T)
an_tot2[,,3,] <- apply(ancitysim_tot_2[,,,-1],1:3,quantile,0.975, na.rm=T)

an_tot3[,,1,] <- ancitysim_tot_3[,,,1]
an_tot3[,,2,] <- apply(ancitysim_tot_3[,,,-1],1:3,quantile,0.025, na.rm=T)
an_tot3[,,3,] <- apply(ancitysim_tot_3[,,,-1],1:3,quantile,0.975, na.rm=T)


# AF USING THE CORRECT DENOMINATOR
af_tot1[,1,,] <- (an_tot1[,1,,]/deaths_tot[1])*100
af_tot1[,2,,] <- (an_tot1[,2,,]/deaths_tot[2])*100
af_tot1[,3,,] <- (an_tot1[,3,,]/deaths_tot[3])*100

af_tot2[,1,,] <- (an_tot2[,1,,]/deaths_tot[1])*100
af_tot2[,2,,] <- (an_tot2[,2,,]/deaths_tot[2])*100
af_tot2[,3,,] <- (an_tot2[,3,,]/deaths_tot[3])*100

af_tot3[,1,,] <- (an_tot3[,1,,]/deaths_tot[1])*100
af_tot3[,2,,] <- (an_tot3[,2,,]/deaths_tot[2])*100
af_tot3[,3,,] <- (an_tot3[,3,,]/deaths_tot[3])*100

# OVERALL
an_tot_over1[,1,] <- apply(ancitysim_tot_1[,,,1],c(1,3),sum)
an_tot_over1[,2,] <- apply(apply(ancitysim_tot_1[,,,-1],c(1,3,4),sum),
                           1:2,quantile,0.025, na.rm=T)
an_tot_over1[,3,] <- apply(apply(ancitysim_tot_1[,,,-1],c(1,3,4),sum),
                           1:2,quantile,0.975, na.rm=T)                           

an_tot_over2[,1,] <- apply(ancitysim_tot_2[,,,1],c(1,3),sum)
an_tot_over2[,2,] <- apply(apply(ancitysim_tot_2[,,,-1],c(1,3,4),sum),
                           1:2,quantile,0.025, na.rm=T)
an_tot_over2[,3,] <- apply(apply(ancitysim_tot_2[,,,-1],c(1,3,4),sum),
                           1:2,quantile,0.975, na.rm=T) 

an_tot_over3[,1,] <- apply(ancitysim_tot_3[,,,1],c(1,3),sum)
an_tot_over3[,2,] <- apply(apply(ancitysim_tot_3[,,,-1],c(1,3,4),sum),
                           1:2,quantile,0.025, na.rm=T)
an_tot_over3[,3,] <- apply(apply(ancitysim_tot_3[,,,-1],c(1,3,4),sum),
                           1:2,quantile,0.975, na.rm=T) 

deaths_tot_over <- sum(deaths_tot)
af_tot_over1 <- (an_tot_over1/deaths_tot_over)*100
af_tot_over2 <- (an_tot_over2/deaths_tot_over)*100
af_tot_over3 <- (an_tot_over3/deaths_tot_over)*100


# ENSEMBLE ESTIMATES OBS / MODELLED / BOTH - ACCOUNT FOR VARIABILITY ACROSS SERIES
# note: to compute the estimate I use mean instead of median because it was more 
# compatible with the single estimates.
# CREATE EMPTY ARRAY TO STORE SUMMARY RESULTS
an_tot_syn1 <- af_tot_syn1 <- an_tot_syn2 <- af_tot_syn2 <- an_tot_syn3 <- 
  af_tot_syn3 <- array(0,dim=c(3,length(d1),3,2)
    ,dimnames=list(c("obs_syn", "model_syn", "obsmod_syn"),
                   paste0("dec", seq(length(d1))),c("est","ci.l","ci.u"),
    c("abs","dif")))

an_tot_over_syn1 <- af_tot_over_syn1 <- an_tot_over_syn2 <- af_tot_over_syn2 <- 
  an_tot_over_syn3 <- af_tot_over_syn3 <-  array(0,dim=c(3,3,2)
    ,dimnames=list(c("obs_syn", "model_syn", "obsmod_syn"),c("est","ci.l","ci.u"),
    c("abs","dif")))

# BY DECADE
an_tot_syn1[1,,1,] <- apply(ancitysim_tot_1[c(9:12),,,-1],c(2,3),mean)
an_tot_syn1[1,,2,] <- apply(ancitysim_tot_1[c(9:12),,,-1],c(2,3),quantile,0.025)
an_tot_syn1[1,,3,] <- apply(ancitysim_tot_1[c(9:12),,,-1],c(2,3),quantile,0.975)

an_tot_syn1[2,,1,] <- apply(ancitysim_tot_1[c(2:8),,,-1],c(2,3),mean)
an_tot_syn1[2,,2,] <- apply(ancitysim_tot_1[c(2:8),,,-1],c(2,3),quantile,0.025)
an_tot_syn1[2,,3,] <- apply(ancitysim_tot_1[c(2:8),,,-1],c(2,3),quantile,0.975)

an_tot_syn1[3,,1,] <- apply(ancitysim_tot_1[c(2:12),,,-1],c(2,3),mean)
an_tot_syn1[3,,2,] <- apply(ancitysim_tot_1[c(2:12),,,-1],c(2,3),quantile,0.025)
an_tot_syn1[3,,3,] <- apply(ancitysim_tot_1[c(2:12),,,-1],c(2,3),quantile,0.975)


an_tot_syn2[1,,1,] <- apply(ancitysim_tot_2[c(9:12),,,-1],c(2,3),mean)
an_tot_syn2[1,,2,] <- apply(ancitysim_tot_2[c(9:12),,,-1],c(2,3),quantile,0.025)
an_tot_syn2[1,,3,] <- apply(ancitysim_tot_2[c(9:12),,,-1],c(2,3),quantile,0.975)

an_tot_syn2[2,,1,] <- apply(ancitysim_tot_2[c(2:8),,,-1],c(2,3),mean)
an_tot_syn2[2,,2,] <- apply(ancitysim_tot_2[c(2:8),,,-1],c(2,3),quantile,0.025)
an_tot_syn2[2,,3,] <- apply(ancitysim_tot_2[c(2:8),,,-1],c(2,3),quantile,0.975)

an_tot_syn2[3,,1,] <- apply(ancitysim_tot_2[c(2:12),,,-1],c(2,3),mean)
an_tot_syn2[3,,2,] <- apply(ancitysim_tot_2[c(2:12),,,-1],c(2,3),quantile,0.025)
an_tot_syn2[3,,3,] <- apply(ancitysim_tot_2[c(2:12),,,-1],c(2,3),quantile,0.975)


an_tot_syn3[1,,1,] <- apply(ancitysim_tot_3[c(9:12),,,-1],c(2,3),mean)
an_tot_syn3[1,,2,] <- apply(ancitysim_tot_3[c(9:12),,,-1],c(2,3),quantile,0.025)
an_tot_syn3[1,,3,] <- apply(ancitysim_tot_3[c(9:12),,,-1],c(2,3),quantile,0.975)

an_tot_syn3[2,,1,] <- apply(ancitysim_tot_3[c(2:8),,,-1],c(2,3),mean)
an_tot_syn3[2,,2,] <- apply(ancitysim_tot_3[c(2:8),,,-1],c(2,3),quantile,0.025)
an_tot_syn3[2,,3,] <- apply(ancitysim_tot_3[c(2:8),,,-1],c(2,3),quantile,0.975)

an_tot_syn3[3,,1,] <- apply(ancitysim_tot_3[c(2:12),,,-1],c(2,3),median)
an_tot_syn3[3,,2,] <- apply(ancitysim_tot_3[c(2:12),,,-1],c(2,3),quantile,0.025)
an_tot_syn3[3,,3,] <- apply(ancitysim_tot_3[c(2:12),,,-1],c(2,3),quantile,0.975)


# AF USING THE CORRECT DENOMINATOR
af_tot_syn1[,1,,] <- (an_tot_syn1[,1,,]/deaths_tot[1])*100
af_tot_syn1[,2,,] <- (an_tot_syn1[,2,,]/deaths_tot[2])*100
af_tot_syn1[,3,,] <- (an_tot_syn1[,3,,]/deaths_tot[3])*100

af_tot_syn2[,1,,] <- (an_tot_syn2[,1,,]/deaths_tot[1])*100
af_tot_syn2[,2,,] <- (an_tot_syn2[,2,,]/deaths_tot[2])*100
af_tot_syn2[,3,,] <- (an_tot_syn2[,3,,]/deaths_tot[3])*100

af_tot_syn3[,1,,] <- (an_tot_syn3[,1,,]/deaths_tot[1])*100
af_tot_syn3[,2,,] <- (an_tot_syn3[,2,,]/deaths_tot[2])*100
af_tot_syn3[,3,,] <- (an_tot_syn3[,3,,]/deaths_tot[3])*100

# OVERALL
an_tot_over_syn1[1,1,] <- apply(apply(ancitysim_tot_1[c(9:12),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn1[1,2,] <- apply(apply(ancitysim_tot_1[c(9:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn1[1,3,] <- apply(apply(ancitysim_tot_1[c(9:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

an_tot_over_syn1[2,1,] <- apply(apply(ancitysim_tot_1[c(2:8),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn1[2,2,] <- apply(apply(ancitysim_tot_1[c(2:8),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn1[2,3,] <- apply(apply(ancitysim_tot_1[c(2:8),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

an_tot_over_syn1[3,1,] <- apply(apply(ancitysim_tot_1[c(2:12),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn1[3,2,] <- apply(apply(ancitysim_tot_1[c(2:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn1[3,3,] <- apply(apply(ancitysim_tot_1[c(2:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)


an_tot_over_syn2[1,1,] <- apply(apply(ancitysim_tot_2[c(9:12),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn2[1,2,] <- apply(apply(ancitysim_tot_2[c(9:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn2[1,3,] <- apply(apply(ancitysim_tot_2[c(9:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

an_tot_over_syn2[2,1,] <- apply(apply(ancitysim_tot_2[c(2:8),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn2[2,2,] <- apply(apply(ancitysim_tot_2[c(2:8),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn2[2,3,] <- apply(apply(ancitysim_tot_2[c(2:8),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

an_tot_over_syn2[3,1,] <- apply(apply(ancitysim_tot_2[c(2:12),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn2[3,2,] <- apply(apply(ancitysim_tot_2[c(2:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn2[3,3,] <- apply(apply(ancitysim_tot_2[c(2:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)


an_tot_over_syn3[1,1,] <- apply(apply(ancitysim_tot_3[c(9:12),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn3[1,2,] <- apply(apply(ancitysim_tot_3[c(9:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn3[1,3,] <- apply(apply(ancitysim_tot_3[c(9:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

an_tot_over_syn3[2,1,] <- apply(apply(ancitysim_tot_3[c(2:8),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn3[2,2,] <- apply(apply(ancitysim_tot_3[c(2:8),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn3[2,3,] <- apply(apply(ancitysim_tot_3[c(2:8),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

an_tot_over_syn3[3,1,] <- apply(apply(ancitysim_tot_3[c(2:12),,,-1],c(1,3,4),sum),
                           2,mean, na.rm=T)
an_tot_over_syn3[3,2,] <- apply(apply(ancitysim_tot_3[c(2:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.025, na.rm=T)
an_tot_over_syn3[3,3,] <- apply(apply(ancitysim_tot_3[c(2:12),,,-1],c(1,3,4),sum),
                           2,quantile,0.975, na.rm=T)

af_tot_over_syn1 <- (an_tot_over_syn1/deaths_tot_over)*100
af_tot_over_syn2 <- (an_tot_over_syn2/deaths_tot_over)*100
af_tot_over_syn3 <- (an_tot_over_syn3/deaths_tot_over)*100


