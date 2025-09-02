################################################################################
################################################################################
##### Refining methods for attributing health 
##### impacts to climate change: a heat-mortality case study in Zürich
#####     by Stuart-Smith et al. 2025 Climatic Change 
################################################################################
################################################################################

################################################################################

# 4. HEAT-MORTATLIY IN SUMMER 2018

################################################################################


################################################################################
# COMPUTATION OF ATTRIBUTABLE HEAT-MORTALITY IN SUMMER 2028

# We use the backward perspective so heat-mortality better follows the daily variation in Tmean.
# Note that the results are very similar between forward and backward perspective. 

source("attrdl.R")

date2018 <- seq(as.Date("2018-06-01"),as.Date("2018-08-31"), by="day")


# CREATE EMPTY ARRAY TO STORE RESULTS
an2018 <- array(NA, dim=list(length(date2018),length(temp_sim),2),
                  dimnames=list(as.character(date2018),temp_sim,c("abs","dif")))

an2018_tot <- array(NA, dim=list(length(temp_sim),3,2,nsim+1),
                  dimnames=list(temp_sim,c("tot","HW","noHW"),c("abs","dif"), 
                                c("est",paste0("sim",seq(nsim)))))


# RUN CODE
for (g in seq(length(temp_sim))){
  
  # SELECT THE 3rd PERIOD
  dta.sub <- subset(dta.attr, year %in% d1[3]:d2[3])
      
  # CREATE CB
  argvar <- list(fun=varfun,knots=quantile(dta.sub$tmean_obs,varper/100,na.rm=T), 
                   Bound=range(dta.sub$tmean_obs,na.rm=T))
  arglag <- list(knots=logknots(lag,lagnk))

  cb <- crossbasis(dta.sub$tmean_obs,lag=lag,
                     argvar=argvar,
                     arglag=arglag, group=as.factor(dta.sub$year)) 
    
  # MODEL FORMULA 
  formula <- deaths ~ cb + dow + ns(yday,df=dfseas):factor(year) + ns(date,df=round(length(unique(year))/dftrend/10))
  
  # RUN MODEL & PREDICT
  model <- glm(formula,dta.sub,family=quasipoisson,na.action="na.exclude")
  cp <- crossreduce(cb,model, by=0.1, at=c(quantile(dta.sub$tmean_obs, 0.20, na.rm=T):max(dta.sub$tmean, na.rm=T)),
                    cen=mmt[3])
    
  set.seed(13041975)

  # STORE DAILY ATT MORTALITY WITH ATTRDL
  an <- attrdl(dta.sub[,temp_sim[g]], cb, dta.sub$death, model, 
                 type="an", dir="back", tot=FALSE, cen=mmt[3], range=c(mmt[3],40))
      
  # SELECT SUMMER 2018
  ind2018 <- year(dta.sub$date)==2018
  
  # STORE DAILY RESULTS
  an2018[,g,1] <- an[ind2018]
    
  # CREATE INDICATOR HW
  indhw <- dta.sub$date %in% seq(as.Date("2018-07-29"),as.Date("2018-08-10"), "day")
  
  # ESTIMATE SUMMARIES WITH ATTRDL
  an2018_tot[g,"tot",1,1] <- attrdl(dta.sub[ind2018,temp_sim[g]], cb, dta.sub$death[ind2018], model, 
                     type="an", dir="back", tot=T
                     , cen=mmt[3], range=c(mmt[3],40))
  an2018_tot[g,"HW",1,1] <- attrdl(dta.sub[ind2018&indhw,temp_sim[g]], cb, dta.sub$death[ind2018&indhw], model, 
                     type="an", dir="back", tot=T
                     , cen=mmt[3], range=c(mmt[3],40))
  an2018_tot[g,"noHW",1,1] <- attrdl(dta.sub[ind2018&!indhw,temp_sim[g]], cb, dta.sub$death[ind2018&!indhw], 
                                       model,type="an", dir="back", tot=T,
                                       cen=mmt[3], range=c(mmt[3],40))

  # COMPUTE SIMULATIONS
  an2018_tot[g,"tot",1,-1] <- attrdl(dta.sub[ind2018,temp_sim[g]], cb, dta.sub$death[ind2018], model, 
                     type="an", dir="back", tot=T
                     , cen=mmt[3], range=c(mmt[3],40), sim=T, nsim=nsim)
  an2018_tot[g,"HW",1,-1] <- attrdl(dta.sub[ind2018&indhw,temp_sim[g]], cb, dta.sub$death[ind2018&indhw], model, 
                     type="an", dir="back", tot=T
                     , cen=mmt[3], range=c(mmt[3],40), sim=T, nsim=nsim)
  an2018_tot[g,"noHW",1,-1] <- attrdl(dta.sub[ind2018&!indhw,temp_sim[g]], cb, 
                                        dta.sub$death[ind2018&!indhw], model, 
                     type="an", dir="back", tot=T
                     , cen=mmt[3], range=c(mmt[3],40), sim=T, nsim=nsim)

}

# COMPUTE DIFFERENCE SUMMARIES
an2018[,,2] <- an2018[,rep(1,length(temp_sim)),1] - an2018[,,1]
an2018_tot[,,2,] <- an2018_tot[rep(1,length(temp_sim)),,1,] - an2018_tot[,,1,]


# COMPUTE SUMMARIES FOR SUMMER 2018
an_2018 <- af_2018 <- array(0,dim=c(2,3,2,3)
    ,dimnames=list(c("fact", "counterfact"),
                   c("tot","HW","noHW"), 
                   c("abs","dif"),
                   c("est","ci.l","ci.u")))

an_2018[1,,,1] <- an2018_tot[1,,,1]
an_2018[1,,,2] <- apply(an2018_tot[1,,,-1],c(1,2), quantile,0.025,na.rm=T)
an_2018[1,,,3] <- apply(an2018_tot[1,,,-1],c(1,2), quantile,0.975,na.rm=T)

an_2018[2,,,1] <- apply(an2018_tot[c(2:12),,,1],c(2,3),mean,na.rm=T)
an_2018[2,,,2] <- apply(apply(an2018_tot[c(2:12),,,-1],c(2:4),mean,na.rm=T), c(1,2),quantile,0.025,na.rm=T)
an_2018[2,,,3] <- apply(apply(an2018_tot[c(2:12),,,-1],c(2:4),mean,na.rm=T), c(1,2),quantile,0.975,na.rm=T)

af_2018[,"tot",,] <- (an_2018[,"tot",,] /sum(dta.sub[ind2018,"deaths"]))*100
af_2018[,"HW",,] <- (an_2018[,"HW",,] /sum(dta.sub[ind2018&indhw,"deaths"]))*100
af_2018[,"noHW",,] <- (an_2018[,"noHW",,] /sum(dta.sub[ind2018&!indhw,"deaths"]))*100



# RESULTS HW2018
table_res_2018_anaf <- matrix(NA, nrow=6, ncol=3)

for(i in c(1:3)){
  table_res_2018_anaf[1,i] <- paste0(round(an_2018["fact",i,"abs","est"],0)," [", 
                                     round(an_2018["fact",i,"abs","ci.l"],0), ";",
                                     round(an_2018["fact",i,"abs","ci.u"],0), "]")
  table_res_2018_anaf[2,i] <- paste0(round(an_2018["counterfact",i,"abs","est"],0)," [", 
                                     round(an_2018["counterfact",i,"abs","ci.l"],0), ";",
                                     round(an_2018["counterfact",i,"abs","ci.u"],0), "]")
  table_res_2018_anaf[3,i] <- paste0(round(an_2018["counterfact",i,"dif","est"],0)," [", 
                                     round(an_2018["counterfact",i,"dif","ci.l"],0), ";",
                                     round(an_2018["counterfact",i,"dif","ci.u"],0), "]")
  table_res_2018_anaf[4,i] <- paste0(round(af_2018["fact",i,"abs","est"],2)," [", 
                                     round(af_2018["fact",i,"abs","ci.l"],2), ";",
                                     round(af_2018["fact",i,"abs","ci.u"],2), "]")
  table_res_2018_anaf[5,i] <- paste0(round(af_2018["counterfact",i,"abs","est"],2)," [", 
                                     round(af_2018["counterfact",i,"abs","ci.l"],2), ";",
                                     round(af_2018["counterfact",i,"abs","ci.u"],2), "]")
  table_res_2018_anaf[6,i] <- paste0(round(af_2018["counterfact",i,"dif","est"],2)," [", 
                                     round(af_2018["counterfact",i,"dif","ci.l"],2), ";",
                                     round(af_2018["counterfact",i,"dif","ci.u"],2), "]")
}

colnames(table_res_2018_anaf) <- c("Total", "HW", "NoHW")

rownames(table_res_2018_anaf) <- c("AN factual", "AN counter", "AN dif",
                                "AF factual", "AF counter", "AF dif")

