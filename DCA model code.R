library(rriskDistributions)
library(ReIns)
library(RColorBrewer)
library(dplyr)
library(tibble)


set.seed(12345)
mycolors <- brewer.pal(n = 9, name = 'Set1')

#### inputs #####
make.params <- function(hosp=T)
{
  params <- list()
  params <- within(params, {
    #weights
    I_T <- 1 # value of preventing transmission from a case
    clinicalVsTransmissionBenefit <- ifelse(hosp, 1, 0.03) # value of preventing mortality/morbidity, in units of T
    # outpatient: 6% of all diagnoses hospitalized (supplement); assume half are already in hospital when tested.
    # Then, ~30x lower risk of severe illness in symptomatic outpatient versus hospitalized patient. Assume outpatient treatment for prevention has similar % morbidity mortality impact as inpatient mgme once sick (Mab trials, reduced for limited uptake).

    falseDiagnosisHarm <- 1/10
    # i_F <- 1/10 # cost of an isolation, relative to value of preventing transmission from a case
    # c_F <- 0.1 #cost of FP vs benefit of TP, clinical
    # s_F <- 0.1 #cost of FP vs benefit of TP, statistical/surveillance
    # n <- 3 #isolated contacts per diagnosis (used for isolation costs, not transmission calculations)
    
    clinicalDiagnosisDiscount <- 0.75 # proportion of TP benefit achievable through clinical diagnosis (e.g. assumes you don't isolate as effectively, don't fully stop loking for another diagnosis, don't get the case notification)
    # hard to base on reduction in PPV, from >90% to prev*sens/(prev*sens+(1-prev)(spec))= 0.09/(0.09+0.9*0.5) = 15% ?

    presx <- 0.2
    incubationPeriod <- 5.5 # incubationPeriod period mean
    inc_sd <- 1.5
    isolationOfKnownCase <- ifelse(hosp, 0.9, 0.7) # transmission averted after diagnosed
    isolationOfCasePendingResult <- ifelse(hosp, 0.7, 0.3) # transmission averted while awaiting diagnostic result (from case, versus asymptomatic or after testing falsely negative)
    isolationOfContacts <- 0.4 #contact notification efficacy for averting future transmission. 
    sxOnsetDay <- 2 # from onset of infectivity, on average
    duration <- 10 # from symptom onset to end of 99% of cumulative infectivity, on average
    
    # viral load and infectivity
    v1=3; v2=8 #95% CI of log10 viral load
    minimumInfectiousLogVirus=3 #minimum logvirus at which infectivity begins
    infectivityScale = 1 # 1 = log-linear -- how skewed is the infectivity distribution, in relation to viral load?
    firstGenerationTransmissionWeight <- 0.25 # how much do you value preventing transmission to immediate contacts,

    # timing
    medianPresentationDay <- 3.5
    turnaroundTimeRDT <- 3/24
    turnaroundTimeNAAT <- ifelse(hosp, 1, 3)
    
    # time sensitivity discount
    dailyDecayClinicalBenefit <- ifelse(hosp, 0.12, 0.18)  # daily discount, clinical result value.
    
    csens_h <- 0.8
    cspec_h <- 0.5 
    csens_o <- csens_h #0.5
    cspec_o <- cspec_h #0.8
    sensitivityClinician <- ifelse(hosp, csens_h, csens_o) # may want to define elsewhere
    specificityClinician <- ifelse(hosp, cspec_h, cspec_o)
    sensitivityRDT_vsNAAT <- 0.85
    specificityRDT <- 0.99
    sensitivityNAAT <- 0.9
    spec_NAAT <- 0.995
    
    prev <- ifelse(hosp, 0.4, 0.1)
    
    weibpar <- get.weibull.par(p=c(0, presx, 0.99), q=c(0,sxOnsetDay,sxOnsetDay+duration), plot = F, show.output = F)
    normpar <- get.norm.par(p = c(0.025, 0.975), q = c(v1,v2), plot = F, show.output = F)
    
  })
  
  return(params)
}


niceParamNames <- names(make.params())
niceParamNames['prev'] <- 
  "Prevalence of COVID-19"
niceParamNames['sensitivityNAAT'] <- 
  "Sensitivity of NAAT"
niceParamNames['sensitivityRDT_vsNAAT'] <- 
  "Sensitivity of RDT (vs NAAT, acute)"
niceParamNames['isolationOfKnownCase'] <- 
  "Isolation after case diagnosis"
niceParamNames['isolationOfCasePendingResult'] <- 
  "Isolation awaiting test result"
niceParamNames['infectivityScale'] <- 
  "Viral burden - infectivity association"
niceParamNames['turnaroundTimeNAAT'] <- 
  "NAAT turnaround time (days)"
niceParamNames['minimumInfectiousLogVirus'] <-
  "Minimum infectious viral burden (log 10)"
niceParamNames['isolationOfContacts'] <- 
  "Isolation of contacts"
niceParamNames['firstGenerationTransmissionWeight'] <- 
  "Importance of 1st-generation transmission"
niceParamNames['sxOnsetDay'] <- 
  "Days from exposure to symptoms"
niceParamNames['turnaroundTimeRDT'] <- 
  "Ag-RDT turnaround time (days)"
niceParamNames['dailyDecayClinicalBenefit'] <- 
  "Time sensivity of clinical intervention"
niceParamNames['specificityRDT'] <- 
  "Specificity of Ag-RDT"
niceParamNames["specificityClinician"] <- 
 "Specificity of clinician judgment"
niceParamNames["sensitivityClinician"] <- 
  "Sensitivity of clinician judgment"
niceParamNames["clinicalDiagnosisDiscount"] <- 
  "Intensity reduction when diagnosed clinically"
niceParamNames["falseDiagnosisHarm"] <- 
  "Relative harm of false positives"
niceParamNames["clinicalVsTransmissionBenefit"] <- 
  "Contribution of clinical versus transmission\neffects to net benefit"
niceParamNames["medianPresentationDay"] <- 
  "Median days from symptoms to presentation"





### utility functions ####

col1 <-  rgb(255,200,200,max = 255, alpha = 70)
col2 <- rgb(173,216,230,max = 255, alpha = 70)


get_diff_vectors <- function(x, y) {
  count_x <- table(x)
  count_y <- table(y)
  same_counts <- match(names(count_y), names(count_x))
  count_x[same_counts] <- count_x[same_counts] - count_y
  as.numeric(rep(names(count_x), count_x))
}


############ detection, infectiousness and transmission #####

#' @description Determines the proportion of transmissions (as a weighted average of transmissions by cases and theri contacts)
  # that can be prevented by using a specified test, accounting for the timing of presentation, the test's turnaround time and sensitivity and 
  # an assumed distribution of infectivity over time and across viral loads relative to tests' limit of detection.
# So this is a proportion of all transmission, not transmission among those detected. 

#' @return proportion of cases detected, proportion of transmission prevented (after weighting by generation).
#' Will assume that full clinical benefit potential remains when cases become severe enogh to present. Because those presenting later may lose some options but may benefit more from the options that remain. 
mapTransmissions <- function(params,
                            npts = 1e5, #arbitrary sim size
                            plots=T, plotclin=F) 
{  
  if(missing(params)) params <- make.params()
  
  with(params, {
  

  # timing of transmissions from index case
  casetransmits <- rweibull(n = npts, shape = weibpar[1], scale = weibpar[2]) - sxOnsetDay

  # incubationPeriod period for secondary case's symptom onset 
  incubation <- rgamma(shape = incubationPeriod^2/inc_sd, rate = incubationPeriod/inc_sd, n = npts)
  contactonset <- casetransmits + incubation
  # using a random resampling of case transmission timing to estimate timing of transmission from contacts, 

  contacttransmits <- contactonset + (rweibull(n = npts, shape = weibpar[1], scale = weibpar[2]) - sxOnsetDay)

  ##  heterogeneous infectivity 
  
  # Will model variation in infectiousness at diagnosis, and assume that corresponds to variation in avertable transmission -- either because the low ones were always low so weren't very infectious and don't ahve many infected contacts, or because you're cathcing them late when they aren't infectious and their contacts' infectious time has mostly passed.)
  logvirus <- rnorm(n = npts, mean=normpar[1], sd = normpar[2]) # at symptom onset (affects detection, correlated with peak)
  
  if(plots)
  {
    ggplot(data.frame(logvirus)) + geom_density(aes(x=logvirus), lwd=0.5) + 
      geom_vline(xintercept=min(logvirus[Pos_NAAT_onset==1]), lty=2, color=mycolors[1], lwd=0.5) + 
      geom_vline(xintercept=min(logvirus[Pos_RDT_onset==1]), lty=3, color=mycolors[2], lwd=0.5) + 
      xlab("Log10 virus at symptom onset") +
      ylab("Proportion of patient population") +
      xlim(1,10)+
      theme_minimal(base_size = 12) +
      scale_y_continuous(breaks = NULL) + 
      coord_flip() + 
      annotate(geom="text", x=min(logvirus[Pos_NAAT_onset==1])+0.1, y=0, 
               label="Limit of NAAT detection", color=mycolors[1], hjust=0, vjust=0, size=4) + 
      annotate(geom="text", x=min(logvirus[Pos_RDT_onset==1])+0.1, y=0, 
               label="Limit of Ag-RDT detection", color=mycolors[2], hjust=0, vjust=0, size=4) 
  #   annotate(geom="text", x=min(logvirus[Pos_NAAT_onset==1])-0.1, y=0, 
  #            label="Limit of NAAT detection", color=mycolors[1], hjust=0, vjust=0) + 
  #     annotate(geom="text", x=min(logvirus[Pos_RDT_onset==1])-0.1, y=0, 
  #              label="Limit of Ag-RDT detection", color=mycolors[2], hjust=0, vjust=0) 
  }

  # infectivity
  
  infectivity <- (pmax(0, pmin(logvirus-minimumInfectiousLogVirus, max(logvirus)-minimumInfectiousLogVirus)/(max(logvirus)-minimumInfectiousLogVirus)))^infectivityScale
  infectivity <- infectivity/mean(infectivity) # this is infectivity relative to average (at any time point, assumes remains proportional)
  
  # detection
  
  sensitivityRDT <- sensitivityRDT_vsNAAT * sensitivityNAAT # at symptom onset
  
  # for each presentation time (in, say, 0.1 day increments), determine the LOD:
  timing <- seq(0,25, by=0.1)
  sensitivityNAAT_delay <- sensitivityNAAT * pmax(0, 1 - pmax(0,(timing-6)/(22-6)))
  sensitivityRDT_delay <- sensitivityRDT * pmax(0, 1 - pmax(0,(timing-6)/(14-6)))
  
  LODs_delay_NAAT <-  quantile(logvirus, 1-sensitivityNAAT_delay)
  LODs_delay_RDT <-  quantile(logvirus, 1-sensitivityRDT_delay)
    
  if(medianPresentationDay==0) presentationDay = 0 else
    presentationDay <- rtweibull(n = npts, shape=0.85, scale=5.5*medianPresentationDay/3.5, endpoint = 25)
  presentationStep <- unlist(lapply(presentationDay, function(x) which.min(timing<x)))
  
  Pos_NAAT_onset <- 1*(logvirus>LODs_delay_NAAT[1])
  Pos_RDT_onset <- 1*(logvirus>LODs_delay_RDT[1])
  
  Pos_NAAT <- 1*(logvirus>LODs_delay_NAAT[presentationStep])
  Pos_RDT <- 1*(logvirus>LODs_delay_RDT[presentationStep])
  Pos_clin <- rbinom(n = npts, size = 1, prob = sensitivityClinician)
  
  if(plots) {

    # detection
    plot(timing, 100*sensitivityNAAT_delay, type='l', 
         xlab="Days since symptom onset", 
         ylab="% symptomatic COVID patients detected",
         col=mycolors[1], lwd=2,
         ylim=c(0,100))
    lines(timing, 100*sensitivityRDT_delay, col=mycolors[2], lwd=2)
    text(x = 0, y=100*sensitivityNAAT+1, "NAAT", col=mycolors[1], adj = c(0,0))
    text(x = 0, y=100*sensitivityRDT+1, "Ag-RDT", col=mycolors[2], adj = c(0,0))
    
    # # old version, infectivity and burden
    # hist(infectivity, col='black', breaks = seq(0,ceiling(max(infectivity)),by=0.1), xlab="Relative infectivity", freq = T,
    #      main="Infectivity and detectability")
    # hist(infectivity[rbinom(npts,1,prob=Pos_NAAT_onset)==1], 
    #      add=T, breaks = seq(0,ceiling(max(infectivity)),by=0.1))
    # hist(infectivity[rbinom(npts,1,prob=Pos_RDT_onset)==1], 
    # # hist(infectivity[rbinom(prob=probPos_RDT)], 
    #      add=T, breaks = seq(0,ceiling(max(infectivity)),by=0.1), col='white')
    # abline(v=mean(infectivity[rbinom(npts,1,prob=Pos_RDT_onset)==1]), lty=2, col=mycolors[2], lwd=2)
    # abline(v=mean(infectivity[rbinom(npts,1,prob=Pos_NAAT_onset - Pos_RDT_onset)==1]), lty=2, col=mycolors[1], lwd=2)
    # legend("topright", c("Not virologically detected", "Detected by NAAT", "Detected by Ag-RDT and NAAT"), fill=c("black","gray","white"))  
    # 
    ratio <- sample(infectivity[Pos_RDT_onset==1], size = 100000, replace = T)/ 
        sample(infectivity[Pos_NAAT_onset & !Pos_RDT_onset], size=100000, replace=T) 
      
    print(mean(ratio)); print(sd(ratio)  )
    
    # plot of infectivity versus peak log viral burden
    v <- seq(0,12,by=0.1)
    i <- (pmax(0, pmin(v-minimumInfectiousLogVirus, max(v)-minimumInfectiousLogVirus)/(max(v)-minimumInfectiousLogVirus)))^infectivityScale
    i <- i/mean(i) # this is infectivity relative to average (at any time point, assumes remains proportional)
    par(mar=c(4,4,1,1))
    # plot(i,v, type='l', xlab="relative infectivity", ylab="peak value of log viral burden")  
    # plot(i,10^v, type='l', xlab="relative infectivity", ylab="peak value of log viral burden")  
    # plot(10^v, i, type='l')
    plot(v,i, type='l', ylab="Relative Infectivity", xlab="Peak value of log viral burden", lwd=2)
    abline(h=mean(infectivity[rbinom(npts,1,prob=Pos_RDT_onset)==1]), lty=2, col=mycolors[2], lwd=2)
    abline(h=mean(infectivity[rbinom(npts,1,prob=Pos_NAAT_onset - Pos_RDT_onset)==1]), lty=2, col=mycolors[1], lwd=2)
    text(x=0, y=0.05+mean(infectivity[rbinom(npts,1,prob=Pos_RDT_onset)==1]), 
         "Mean infectivity, patients detectable\n by Ag-RDT at symptom onset", col=mycolors[2], adj=c(0,0))
    text(x=0, y=0.05+mean(infectivity[rbinom(npts,1,prob=Pos_NAAT_onset - Pos_RDT_onset)==1]), 
         "Mean infectivity, patients detectable\n only by NAAT  at symptom onset", col=mycolors[1], adj=c(0,0))
    
  }
  
  ##transmission events, with weighting and LOD
    
  # To get the transmission directly from cases, weight each case by relative infectivity 
  # Get indexes of cases who could generate transmission
  sourceCases <- sample(1:length(casetransmits), size = npts, replace=T, prob = infectivity)

  caseevents <- enframe(casetransmits[sourceCases])
  caseevents$Pos_RDT <- Pos_RDT[sourceCases]
  caseevents$Pos_NAAT <- Pos_NAAT[sourceCases]
  caseevents$Pos_clin <- Pos_clin[sourceCases]
  
  mean(Pos_RDT[sourceCases]) # proportion of transmission sources detectable, not propotion of patients
  mean(Pos_NAAT[sourceCases])
  
  # to get the probability that each transmission event is averted by a given test, 
  # screen by TAT and a binomial probability that the assay detected them.

  caseevents$Avertable_RDT <- caseevents$value > presentationDay + turnaroundTimeRDT
  caseevents$Avertable_RDTpending <- caseevents$value > presentationDay & caseevents$value < presentationDay + turnaroundTimeRDT
  caseevents$Averted_RDT <- caseevents$Avertable_RDT * rbinom(n=nrow(caseevents), # incorporating test sensitivity and isolation effectiveness
                                                                   size=1, prob = isolationOfKnownCase*caseevents$Pos_RDT) + 
                          caseevents$Avertable_RDTpending * rbinom(n=nrow(caseevents), 
                                    size=1, prob = isolationOfCasePendingResult)

  caseevents$Avertable_NAAT <- caseevents$value > presentationDay + turnaroundTimeNAAT
  caseevents$Avertable_NAATpending <- caseevents$value > presentationDay & caseevents$value < presentationDay + turnaroundTimeNAAT
  caseevents$Averted_NAAT <- caseevents$Avertable_NAAT * rbinom(n=nrow(caseevents), # incorporating test sensitivity and isolation effectiveness
                                                              size=1, prob = isolationOfKnownCase*caseevents$Pos_NAAT) + 
                            caseevents$Avertable_NAATpending * rbinom(n=nrow(caseevents), 
                                                                     size=1, prob = isolationOfCasePendingResult)
                          
  caseevents$Avertable_clin <- caseevents$value > presentationDay + 0
  caseevents$Averted_clin <- caseevents$Avertable_clin * rbinom(n=nrow(caseevents), # incorporating test sensitivity and isolation effectiveness
                                                               size=1, prob = isolationOfKnownCase*caseevents$Pos_clin)
  
  ## Contact trnasmission events timing distribution:
  
  # the transmissions that aren't avertable have contact infections that aren't averted and (have transmission after TAT or are missed by falseDiagnosisHarm)
  
  caseevents$contacttime <- contacttransmits # if contact is infected and transmits one infection, when does that transmission occur relative to index case symptom onset?
  
  caseevents$Transmits_RDT <- 1 - caseevents$Averted_RDT # which contacts become infected despite diagnostic and intevention efforts
  # And which transmissions from infected contacts could be stopped (not including those stopped before contacts get infected)
  caseevents$contactTransmissionAvertable_RDT <- caseevents$Transmits_RDT & caseevents$contacttime >presentationDay + turnaroundTimeRDT
  caseevents$contactTransmissionAverted_RDT <- caseevents$contactTransmissionAvertable_RDT * rbinom(n=nrow(caseevents), 
                                                                            size=1, prob = isolationOfContacts*Pos_RDT)
  caseevents$contactTransmits_RDT <- caseevents$Transmits_RDT & !caseevents$contactTransmissionAverted_RDT

  caseevents$Transmits_NAAT <- 1 - caseevents$Averted_NAAT
  caseevents$contactTransmissionAvertable_NAAT <- caseevents$Transmits_NAAT & caseevents$contacttime >presentationDay + turnaroundTimeNAAT
  caseevents$contactTransmissionAverted_NAAT <- caseevents$contactTransmissionAvertable_NAAT * rbinom(n=nrow(caseevents), 
                                                                                                    size=1, prob = isolationOfContacts*Pos_NAAT)
  caseevents$contactTransmits_NAAT <- caseevents$Transmits_NAAT & !caseevents$contactTransmissionAverted_NAAT
  
  caseevents$Transmits_clin <- 1 - caseevents$Averted_clin
  caseevents$contactTransmissionAvertable_clin <- caseevents$Transmits_clin & caseevents$contacttime >presentationDay + 0
  caseevents$contactTransmissionAverted_clin <- caseevents$contactTransmissionAvertable_clin * rbinom(n=nrow(caseevents), 
                                                                                                      size=1, prob = isolationOfContacts*Pos_clin)
  caseevents$contactTransmits_clin <- caseevents$Transmits_clin & !caseevents$contactTransmissionAverted_clin
  
  # casetransmits are timing, caseTransmissions are indexes fo all events (some repeated), 
   # caseTransmits_X are indexed for events that occur with a given test (subset of caseTransmissions)
  
  caseevents$Transmission <- ifelse(caseevents$Transmits_RDT==0, ifelse(caseevents$Transmits_NAAT==0, "Prevented with either test", "Prevented with Ag-RDT only"), 
                                 ifelse(caseevents$Transmits_NAAT==0, "Prevented with NAAT only", "Occurs regardless of test"))
  caseevents %>% count(Transmits_NAAT, Transmits_RDT, Transmission)
  caseevents$Transmission <- factor(caseevents$Transmission, levels=c("Occurs regardless of test",
    "Prevented with NAAT only",
    "Prevented with either test",
    "Prevented with Ag-RDT only"))

  if(plots){
  
    # to generate this figure in an interpretable way, would need to change presentation day to fixed (e.g 3):
  # gg <- ggplot(data=caseevents, aes(x=value, fill=Transmission)) + 
  #   geom_histogram( alpha=0.6, bins = 60) + 
  #     xlim(-3, 16) + 
  #     xlab("Days after symptom onset") + ylab("Number of transmission events") +
  #   annotate(
  #     geom = "curve", xend = c(1,1.25,4)-0.25, yend = c(0,0,0), x = c(7,8,9), y = c(2500,2000,1500), 
  #     curvature = .3, arrow = arrow(length = unit(2, "mm"))
  #   ) + 
  #   annotate("text", x=c(7,8,9), y=c(2500,2000,1500), hjust="left",
  #        label = c("Presumptive isolation", "Ag-RDT resulted", "NAAT resulted"))
  # 
    
  # by timing of transmission, cases
  arranged <- caseevents %>% 
    arrange(value) %>% 
    mutate(rn = row_number())
  arranged_RDT <- subset(caseevents, Transmits_RDT==1) %>% 
    arrange(value) %>% 
    mutate(rn = row_number())
  arranged_NAAT <- subset(caseevents, Transmits_NAAT==1) %>% 
    arrange(value) %>% 
    mutate(rn = row_number())
  arranged_clin <- subset(caseevents, Transmits_clin==1) %>% 
    arrange(value) %>% 
    mutate(rn = row_number())
  
  # by timing of transmission, contacts
  gen2 <- caseevents %>% 
    arrange(contacttime) %>% 
    mutate(rn = row_number())
  gen2_RDT <- subset(caseevents, contactTransmits_RDT==1) %>% 
    arrange(contacttime) %>% 
    mutate(rn = row_number())
  gen2_NAAT <- subset(caseevents, contactTransmits_NAAT==1) %>% 
    arrange(contacttime) %>% 
    mutate(rn = row_number())
  gen2_clin <- subset(caseevents, contactTransmits_clin==1) %>% 
    arrange(contacttime) %>% 
    mutate(rn = row_number())
  
  # by timing of presentation, for all cases not weighted as source cases**
  allcases <- data.frame(presentationDay)
  allcases$RDTDay <- presentationDay + turnaroundTimeRDT
  allcases$NAATDay <- presentationDay + turnaroundTimeNAAT
  allcases$Pos_RDT <- Pos_RDT
  allcases$Pos_NAAT <- Pos_NAAT
  allcases$Pos_clin <- Pos_clin
  tested <- allcases %>% 
    arrange(presentationDay) %>% 
    mutate(rn = row_number())
  detected_RDT <- subset(allcases, Pos_RDT==1) %>% 
    arrange(RDTDay) %>% 
    mutate(rn = row_number())
  detected_NAAT <- subset(allcases, Pos_NAAT==1) %>% 
    arrange(NAATDay) %>% 
    mutate(rn = row_number())
  detected_clin <- subset(allcases, Pos_clin==1) %>% 
    arrange(presentationDay) %>% 
    mutate(rn = row_number())

  
  colors <- c(
    # "Any/All" = mycolors[5],
              "NAAT" = mycolors[1], "Ag-RDT" = mycolors[2], 
              "Clinical" = mycolors[3])
  linetypes <- c("Diagnostic result available" = "dotted", 
              "COVID-19 diagnosed" = 'solid')
  
  (detection <- ggplot() +
    xlim(-3,20) + 
      geom_step(data=detected_clin[seq(1,nrow(detected_clin), by=100),], aes(x=presentationDay, y=rn, color='Clinical', linetype='COVID-19 diagnosed'),lwd=0.8) +
      geom_step(data=tested[seq(1,nrow(tested), by=100),], aes(x=presentationDay, y=rn, color="Clinical", linetype="Diagnostic result available"),lwd=0.8) +
      geom_step(data=tested[seq(1,nrow(tested), by=100),], aes(x=presentationDay+turnaroundTimeRDT, y=rn, color="Ag-RDT", linetype="Diagnostic result available"),lwd=0.8) +
    geom_step(data=tested[seq(1,nrow(tested), by=100),], aes(x=presentationDay+turnaroundTimeNAAT, y=rn, color="NAAT", linetype="Diagnostic result available"),lwd=0.8) +
    geom_step(data=detected_RDT[seq(1,nrow(detected_RDT), by=100),], aes(x=presentationDay+turnaroundTimeRDT, y=rn, color="Ag-RDT", linetype="COVID-19 diagnosed"),lwd=0.8) +
    geom_step(data=detected_NAAT[seq(1,nrow(detected_NAAT), by=100),], aes(x=presentationDay+turnaroundTimeNAAT, y=rn, color='NAAT', linetype="COVID-19 diagnosed"),lwd=0.8)  +
    labs(x = "Days from symptom onset",
         color = "Diagnostic approach", linetype="Step of diagnosis") +
    # scale_color_manual(values = c(mycolors[5], mycolors[2], mycolors[1], mycolors[3]), breaks = c("Any/All", "Ag-RDT", "NAAT", "Clinical")) + 
      scale_color_manual(values = c(mycolors[2], mycolors[1], mycolors[3]), breaks = c("Ag-RDT", "NAAT", "Clinical")) + 
    scale_linetype_manual(values = c("dotted", "solid"), breaks = c("Diagnostic result available", "COVID-19 diagnosed")) + 
    scale_y_continuous(name = "Cumulative evaluations completed or true-positive results,\namong 100k simulated COVID-19 patients", 
                       breaks = npts*c(0,0.25,0.5,0.75,1), 
                       labels = c("0","25k","50k","75k","100k")) + 
      theme_bw()+
      theme(legend.key.width=unit(1,"cm"))
  )
  
  # plot transmission
  colors <- c("NAAT" = mycolors[1], "Ag-RDT" = mycolors[2], "Clinical diagnosis" = mycolors[3], 
              "No intervention" = mycolors[4])
  linetypes <- c("From index cases" = "dashed", "From contacts" = "dotdash")
  (transmission <- ggplot() +
    xlim(-3,20) + 
    geom_step(data=arranged[seq(1,nrow(arranged), by=100),], aes(x=value, y=rn, color="No intervention", linetype="From index cases"),size=0.8) +
      geom_step(data=arranged_clin[seq(1,nrow(arranged_clin), by=100),], aes(x=value, y=rn, color='Clinical diagnosis', linetype="From index cases"),size=0.8) +
      geom_step(data=gen2_clin[seq(1,nrow(gen2_clin), by=100),], aes(x=contacttime, y=rn, color='Clinical diagnosis', linetype="From contacts"),size=0.8) +
      geom_step(data=arranged_RDT[seq(1,nrow(arranged_RDT), by=100),], aes(x=value, y=rn, color="Ag-RDT", linetype="From index cases"),size=0.8) +
    geom_step(data=arranged_NAAT[seq(1,nrow(arranged_NAAT), by=100),], aes(x=value, y=rn, color='NAAT', linetype="From index cases"),size=0.8)  +
    geom_step(data=gen2[seq(1,nrow(gen2), by=100),], aes(x=contacttime, y=rn, color="No intervention", linetype="From contacts"),size=0.8) +
    geom_step(data=gen2_RDT[seq(1,nrow(gen2_RDT), by=100),], aes(x=contacttime, y=rn, color="Ag-RDT", linetype="From contacts"),size=0.8) +
    geom_step(data=gen2_NAAT[seq(1,nrow(gen2_NAAT), by=100),], aes(x=contacttime, y=rn, color='NAAT', linetype="From contacts"),size=0.8)  +
    labs(x = "Days from symptom onset",
        color = "Diagnostic approach",
         linetype="Type of transmission event") +
    scale_color_manual(values = c(mycolors[4],mycolors[2],mycolors[1],mycolors[3]), breaks=c("No intervention","Ag-RDT", "NAAT", "Clinical diagnosis") )+
    scale_linetype_manual(values = c("dashed","dotdash"), breaks=c("From index cases", "From contacts") ) +
    scale_y_continuous(name = "Cumulative transmission events", 
                         breaks = npts*c(0,0.25,0.5,0.75,1), 
                         labels = c("0%","25%","50%","75%","100%")) + 
      theme_bw() +
      theme(legend.key.width=unit(1,"cm"))
    )
  
  library(patchwork)
  detection / transmission
  
  # tables of frequencies, this plot as curve
    tempallcase <- table(floor(caseevents$value*4)/4)
    allcase <- append(0, tempallcase); names(allcase) <- append(-2.25, names(tempallcase))
    tempRDTcase <- table((floor(caseevents$value*4)/4)[caseevents$Transmits_RDT==1])
    RDTcase <- append(0, tempRDTcase); names(RDTcase) <- append(-2.25, names(tempRDTcase))
    tempNAATcase <- table((floor(caseevents$value*4)/4)[caseevents$Transmits_NAAT==1])
    NAATcase <- append(0, tempNAATcase); names(NAATcase) <- append(-2.25, names(tempNAATcase))
    
    allcontact <- table(floor(caseevents$contacttime))
    RDTcontact <- table(floor(caseevents$contacttime[caseevents$contactTransmits_RDT==1]))
    NAATcontact <- table(floor(caseevents$contacttime[caseevents$contactTransmits_NAAT==1]))
    
    par(mfrow=c(2,1))
    plot(names(allcase), allcase, type='l', xlab="Days from symptom onset", 
         ylab="Transmission from cases", xlim=c(-3,30), yaxt='n',
         main="Timing of transmission events")
    polygon(names(allcase), allcase, col=rgb(0.5,0.5,0.5, alpha=0.3)) 
    polygon(names(RDTcase), RDTcase, col=rgb(0,0,1,0.3))
    polygon(names(NAATcase), NAATcase, col=rgb(1,0,0,0.3))
    lines(names(RDTcase), RDTcase, type='l', col=mycolors[2], lwd=2)
    lines(names(NAATcase), NAATcase, type='l', col=mycolors[1], lwd=2)
    legend('topright', lty=1, lwd=2, col=c('black',mycolors[1],mycolors[2]), legend=c("No intervention", "NAAT", "Ag-RDT"))
    # hist(caseevents$contacttime, col="white", main="Directly from cases", xlab="Day relative to case symptom onset",
    #      freq=T, xlim=c(-2,40))
    # hist(caseevents$contacttime[caseevents$contactTransmits_RDT==1], add=T, col=col1,
    #      freq=T, xlim=c(-2,40))
    # hist(caseevents$contacttime[caseevents$contactTransmits_NAAT==1], add=T, col=col2,
    #      freq=T, xlim=c(-2,40))
    plot(names(allcontact), allcontact, type='l', xlab="Days from index case's symptom onset", 
         ylab="Transmission from contacts", xlim=c(-3,30), yaxt='n')
    polygon(names(allcontact), allcontact, col=rgb(0.5,0.5,0.5, alpha=0.3)) 
    polygon(names(RDTcontact), RDTcontact, col=rgb(0,0,1,0.3))
    polygon(names(NAATcontact), NAATcontact, col=rgb(1,0,0,0.3))
    lines(names(RDTcontact), RDTcontact, type='l', col=mycolors[2], lwd=2)
    lines(names(NAATcontact), NAATcontact, type='l', col=mycolors[1], lwd=2)
    
  } else gg <- NA
    
    # transmission prevented from cases
    print(c(
      1-mean(caseevents$Transmits_NAAT),
      1-mean(caseevents$Transmits_RDT),
      1-mean(caseevents$Transmits_clin),
      (1-mean(caseevents$Transmits_clin))*clinicalDiagnosisDiscount,
    
      # transmission prevented from contacts
      1-mean(caseevents$contactTransmits_NAAT),
      1-mean(caseevents$contactTransmits_RDT),
      1-mean(caseevents$contactTransmits_clin),
      (1-mean(caseevents$contactTransmits_clin))*clinicalDiagnosisDiscount
    ))
    
    # 
    # ## for Overview figure, timing of transmission:
    # 
    par(mfrow=c(1,1), mar=c(3,2,1,1), oma=c(0,0,0,0))
    plot(NA, xlab="Days since symptom onset",
             xlim=c(-3,25), yaxt='n', ylab="",
         ylim=c(0,4.5e3))
    rect(xleft = 3, xright = 30, ybottom = 0, ytop = 4.5e4, col = rgb(0.5,0.5,0.5,1/4), border = NA)
    rect(xleft = 6, xright = 30, ybottom = 0, ytop = 4.5e4, col = rgb(0.5,0.5,0.5,1/4), border=NA)
    lines(names(allcase), allcase, col=mycolors[5], lwd=2)
    lines(names(allcontact), 0.25*allcontact, col=mycolors[4], lwd=2)
    text(5,4.2e4, "From index cases", col=mycolors[5])
    text(16,2.75e4, "From direct contacts", col=mycolors[4])
    segments(x0 = 0, x1=0, y0=0, y1 = allcase['0'], lty=2, lwd=1.5)
    text(-0.8,500, "Symptom onset", srt=90, adj=0)
    segments(x0 = 3, x1=3, y0=0, y1 = allcase['3'], lty=3, lwd=1.5)
    # text(2.2,500, "Early diagnosis", srt=90, adj=0)
    segments(x0 = 6, x1=6, y0=0, y1 = allcase['6'], lty=3, lwd=1.5)
    # text(5.2,500, "Later diagnosis", srt=90, adj = 0)
    mtext("Days since symptom onset", side = 1, line=2)
    # 
    
    
    
    # allv <- table(floor(logvirus*10)/10)
    # plot(seq(2,11,by=0.1), allv[as.character(seq(2,11,by=0.1))], type='l', xaxt='n', yaxt='n')
    # polygon(x = c(c(floor(LODs_delay_NAAT[30]*10)/10),seq(floor(LODs_delay_NAAT[30]*10)/10,12,by=0.1)),c(0, allv[as.character(seq(floor(LODs_delay_NAAT[30]*10)/10,12,by=0.1))]), col=rgb(0.5,0.5,0.5, alpha=0.2)) 
    # polygon(x = c(c(floor(LODs_delay_RDT[30]*10)/10+1),seq(floor(LODs_delay_RDT[30]*10)/10+1,12,by=0.1)),c(0, allv[as.character(seq(floor(LODs_delay_RDT[30]*10)/10+1,12,by=0.1))]), col=rgb(0.5,0.5,0.5, alpha=0.4)) 
    # text(x=LODs_delay_NAAT[30]+0.5, y=300, "Detectable only\nby high-sensitivity assays", srt=90, adj = 0)
    # text(x=LODs_delay_RDT[30]+1.5+0.5, y=300, "Detectable by high- or\nlow-sensitivity\nassays", srt=90, adj=0)
    # mtext("Log viral burden (correlated with infectivity)", side=1)
    # mtext("Proportion of patients", side=2)
    # 
    par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,0,0,0))
    plot(NA, xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n')
    mtext("Certainty of diagnosis", side=1)
    mtext("Intensity of intervention", side=2)
    # 
  # Report transmission outcome, relative to reference of preventing all transmission,  
  # and here, scale contact transmission by gen1.
  
  transmission_benefit_RDT <- 1 - 
   ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_RDT) + 
    (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_RDT) )
  
  transmission_benefit_NAAT <- 1 - 
    ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_NAAT) + 
        (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_NAAT) )
  
  transmission_benefit_clin <- 1 - 
    ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_clin) + 
        (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_clin) ) 
  
  # clinical benefit, exponential decay
  clinical_benefit_RDT <- mean(allcases$Pos_RDT * exp(-(allcases$presentationDay + turnaroundTimeRDT) * dailyDecayClinicalBenefit))
  clinical_benefit_NAAT <- mean(allcases$Pos_NAAT * exp(-(allcases$presentationDay + turnaroundTimeNAAT) * dailyDecayClinicalBenefit))
  clinical_benefit_clin <- mean(allcases$Pos_clin * exp(-(allcases$presentationDay) * dailyDecayClinicalBenefit))
  
  return(list("NAAT"=transmission_benefit_NAAT, "RDT"=transmission_benefit_RDT, "clin"=transmission_benefit_clin,
              "NAATclinical"=clinical_benefit_NAAT, "RDTclinical"=clinical_benefit_RDT, "clinclinical"=clinical_benefit_clin,
              "NAATdetection"=mean(Pos_NAAT), "RDTdetection"=mean(Pos_RDT), "clindetection"=mean(Pos_clin),
              "NAATpeakinf" = sum(infectivity[Pos_NAAT==1])/sum(infectivity),
              "RDTpeakinf" = sum(infectivity[Pos_RDT==1])/sum(infectivity), 
              "clinpeakinf" = sum(infectivity[Pos_clin==1])/sum(infectivity),
              "NAATinf" = sum(infectivity[Pos_NAAT==1])/sum(infectivity),
              "RDTinf" = sum(infectivity[Pos_RDT==1])/sum(infectivity), 
              "clininf" = sum(infectivity[Pos_clin==1])/sum(infectivity),
              "NAATavertable" = sum(subset(caseevents, Pos_NAAT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "RDTavertable" = sum(subset(caseevents, Pos_RDT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "clinavertable" = sum(subset(caseevents, Pos_clin==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "NAATaverted" = sum(subset(caseevents, Averted_NAAT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "RDTaverted" = sum(subset(caseevents, Averted_RDT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "clinaverted" = sum(subset(caseevents, Averted_clin==1)$Avertable_clin)/sum(caseevents$Avertable_clin)
              
              ) ) 
  })
  
}

mapTransmissions(params = make.params(hosp=F))
mapTransmissions(params = make.params(hosp=T))

# plot potential clinical impact,  avertable M&M:
timing <- seq(0,30,by=0.1)
hparams <- make.params(hosp=T)
cparams <- make.params(hosp=F)
hmm <- ( hparams$clinicalVsTransmissionBenefit * exp(-hparams$dailyDecayClinicalBenefit * timing)  )
cmm <- ( cparams$clinicalVsTransmissionBenefit * exp(-cparams$dailyDecayClinicalBenefit * timing)  )
par(mar=c(3,4,1,1))
plot(timing, hmm, type='l', col=mycolors[7], lwd=2,
     xlab="",
     ylab="")  
mtext("Value of avertible morbidity and mortality", side=2, line=3, adj=0)
mtext("benchmarked to value of preventing all transmission", side=2, line=2, cex=0.8)
mtext("Days since symptom onset", side=1, line=2)
lines(timing, cmm, col=mycolors[8], lwd=2)
text(10,0.50,labels = "Hospital", col=mycolors[7])
text(10,0.05,labels = "Outpatient", col=mycolors[8])




########## Simple DCA #######
#  Simple decision curve assuming a fixed clinician sens/spec and assuming all cases are equally important to detect

# plot a DCA for some test with fixed sens and spec:
x <- seq(0,1,by=0.01) # threshold probabilities
weights <- x/(1-x) # number of false positives you'd want to treat per case, at each threshold prob

NB_simple <- function(sens, spec, prev)
{
  TP <- sens*prev
  FP <- (1-spec)*(1-prev)
  NB <- TP - FP*weights
  return(NB)
}


plotDCA <- function(prev=0.4, sensitivityClinician=0.9, specificityClinician=0.5, 
                    sensitivityRDT_vsNAAT=0.85, specificityRDT=0.99, 
                    sensitivityNAAT=0.9, spec_NAAT=0.99, 
                    plotclin=T, ymax, plotlegend=T, title="", ylab="Net benefit", lwd=1) 
{
  if(missing(ymax)) ymax <- prev
  sensitivityRDT <- sensitivityRDT_vsNAAT*sensitivityNAAT
  plot(x,NB_simple(sensitivityNAAT, spec_NAAT, prev), type='l', col=mycolors[1], xlim=c(0,1), ylim=c(0,ymax),
       ylab=ylab, xlab="",
       main=title, cex.main=0.9,
       xaxt='n', lwd=lwd)
  axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
  lines(x, NB_simple(sensitivityRDT, specificityRDT, prev), col=mycolors[2], lwd=lwd)
  if(plotclin)  lines(x, NB_simple(sensitivityClinician, specificityClinician, prev), col=mycolors[3], lwd=lwd)
  lines(x, NB_simple(sens=1, spec=0, prev), col='black', lty=1, lwd=lwd)
  lines(x, NB_simple(sens=0, spec=1, prev), col='black', lty=2, lwd=lwd)
  if(plotlegend)legend(x = "topright", legend = c("NAAT","Ag-RDT","Clinical judgment","Treat all", "Treat none"),
                       col=c(mycolors[1],mycolors[2],mycolors[3],'black','black'), lty=c(1,1,1,1,2), lwd=lwd
  )
}

### Plot conventional DCA ###  
layout(array(c(1,1,1,2,2,2,2), dim=c(7,1)))
par(mar=c(2,4,3,1), oma=c(2,0,0,0))
params <- make.params()
params$sensitivityRDT <- params$sensitivityRDT_vsNAAT * params$sensitivityNAAT
nsens_all <- 0.8
rsens_all <- 0.625
rsens_all_vsNAAT <- 0.78

attach(params)
plotDCA(prev=0.1, sensitivityClinician=csens_o, specificityClinician = cspec_o,  
        sensitivityRDT_vsNAAT=rsens_all_vsNAAT, specificityRDT=specificityRDT, 
        sensitivityNAAT=nsens_all, spec_NAAT=spec_NAAT, ymax=0.1, 
        title="A. Outpatient Setting", ylab="Net benefit (conventional approach)", lwd=1.5)
plotDCA(prev=0.4, sensitivityClinician=csens_h, specificityClinician = cspec_h, 
        sensitivityRDT_vsNAAT=rsens_all_vsNAAT, specificityRDT=specificityRDT, 
        sensitivityNAAT=nsens_all, spec_NAAT=spec_NAAT, ymax=0.4, plotlegend = F, 
        title="B. Hospital Setting", ylab="Net benefit (conventional approach)", lwd=1.5)
mtext(side=1, outer=T, text = "Threshold probability above which one would opt to intervene", line=0.5, cex=0.7)
x[which(NB_simple(rsens_all, specificityRDT, prev) > NB_simple(sensitivityClinician, specificityClinician, prev))]
x[which(NB_simple(rsens_all, specificityRDT, prev) > NB_simple(sens=1, spec=0, prev))]
x[which(NB_simple(nsens_all, spec_NAAT, prev) > NB_simple(sens=1, spec=0, prev))]
x[which(NB_simple(rsens_all, specificityRDT, prev=0.1) > NB_simple(sensitivityClinician, specificityClinician, prev = 0.1))]
x[which(NB_simple(rsens_all, specificityRDT, prev=0.1) > NB_simple(sens=1, spec=0, prev = 0.1))]
detach(params)
rm(params)


  
##### Calculate Net Benefit with modified formula ######

# Net benefit, as a function of test and setting/prevalence and other input parameters incl weights
# Will want to display this across multiple dimensions of variation, perhaps up to two at a time

# get the net benefit for a single set of parameters
# the value calculated is the average across the cases that are detected, or the average across all cases? Need to modify to be those detected. 
NB <- function(params)
{
  #need to base this on the proportion actually detected, after accounting for delays. mapTransmissions will account for the lack of trasmission benefit, but we also need to track who's detected at all and gets clinical benefit.
  
    with(params, 
       {   
         mt <- mapTransmissions(params, plots = F)
         
         trans <- mt[1:3]
         clinical <- mt[4:6]
         TAT <- list(NAAT = turnaroundTimeNAAT, RDT = turnaroundTimeRDT, clin=0)
        sens <- list(NAAT = mt$NAATdetection, RDT=mt$RDTdetection, clin=mt$clindetection)
        spec <- list(NAAT=spec_NAAT, RDT = specificityRDT, clin=specificityClinician)
      
        # note that trans is a proportion of transmission prevented among all patients, not only those detected. 
        # So we can divide by sensitivity to get transmission prevented per detection. 
        transaverted <- as.list(mapply(function(a,b) a/b, trans, sens))
        
        clinaverted <- as.list(mapply(function(a,b) a/b, clinical, sens))
        
        # now combine all sources of TP benefit, and discount action in response to clinical decision
        V_TP <- (mapply(function(x,w,z) I_T *
                          # ( x + clinicalVsTransmissionBenefit * exp(-dailyDecayClinicalBenefit * y)  )/ 
                          ( x + w)  / 
                          (1+clinicalVsTransmissionBenefit) * # normalizing so that one full TP get a 1 on both NB scales
                          z, 
                      x=transaverted, w = clinaverted, z=list(1,1,clinicalDiagnosisDiscount)))
        
        V_FP <-  I_T * falseDiagnosisHarm * c(1,1,clinicalDiagnosisDiscount)
        #( (-c_F*clinicalVsTransmissionBenefit ) - (n*isolationOfContacts+1)*i_F)
        # removed factor (1+clinicalVsTransmissionBenefit), to normalize to NB_simple
        
        NB <- prev*V_TP*unlist(sens) - (1-prev)*V_FP*(1-unlist(spec))
        
        NB <- append(NB, rep(NA,2))
        
  
        # treat all, with clinical discount
        NB[5] <- prev* (I_T * ( (averted$clin + clinicalVsTransmissionBenefit)) * clinicalDiagnosisDiscount ) / (1+clinicalVsTransmissionBenefit )- 
          (1-prev)*(I_T * falseDiagnosisHarm)*clinicalDiagnosisDiscount
        
      return(NB)
  })
}


twowayvary <- function(p1, r1, p2=NA, r2=NA, 
                       resetparams # change if needed, for any that will be fixed at a nondefault value
)
{
  # decide which parameters will vary and over what range
  if(missing(resetparams)) resetparams <- make.params() 
  vary1 <- p1
  range1 <- r1
  vary2 <- p2
  range2 <- r2 
  collateresults <- array(NA, dim=c(5, length(range1), ifelse(is.na(vary2), 1, length(range2))))
  if(!is.na(vary1)) 
    for (set1 in 1:length(range1))
    {
      params <- resetparams
      params[[vary1]] <- range1[set1]
      if(!is.na(vary2)) {
        for (set2 in 1:length(range2))
        {  params[[vary2]] <- range2[set2] 
        collateresults[,set1,set2] <- NB(params)
        }} else collateresults[,set1,1] <- NB(params)
    }
  return(collateresults)
}


###### Decision curves for threshold probability, hospital #######3

vary1 <- "falseDiagnosisHarm"
# if we want probability threshold pt to range from 0 to 1, then 
# pt = falseDiagnosisHarm*(1-pt)
pt <- seq(0.001,0.991,by=0.03)
range1 <- pt/(1-pt)

#hospital
probthreshold <- twowayvary(vary1, range1) # standard param values
altparams <- make.params(); altparams$clinicalDiagnosisDiscount <- 1
probthreshold2 <- twowayvary(vary1, range1, resetparams = altparams)
altparams_sens <- altparams_tat <- altparams
altparams_sens$sensitivityRDT_vsNAAT <- 0.95 # for acute infection, 90% in Kruger, 95% in Igloi
altparams_tat$turnaroundTimeNAAT <- 2 # for acute infection, 90% in Kruger, 95% in Igloi
probthreshold_sens <- twowayvary(vary1, range1, resetparams = altparams_sens) # standard param values
probthreshold_tat <- twowayvary(vary1, range1, resetparams = altparams_tat) # standard param values
#outpatient
altparams <- make.params(hosp = F); 
probthreshold4 <- twowayvary(vary1, range1, resetparams = altparams)
altparams$clinicalDiagnosisDiscount <- 1
probthreshold3 <- twowayvary(vary1, range1, resetparams = altparams)
altparams_sens_outpt <- altparams_tat_outpt <- altparams
altparams_sens_outpt$sensitivityRDT_vsNAAT <- 0.95 # for acute infection, 90% in Kruger, 95% in Igloi
altparams_tat_outpt$turnaroundTimeNAAT <- 2 # for acute infection, 90% in Kruger, 95% in Igloi
probthreshold_sens_outpt <- twowayvary(vary1, range1, resetparams = altparams_sens_outpt) # standard param values
probthreshold_tat_outpt <- twowayvary(vary1, range1, resetparams = altparams_tat_outpt) # standard param values

# plot it (Figure 4)
layout(t(array(c(1,2,1,2,1,2,3,4,3,4,3,4,3,4), dim=c(2,7))))
# par(mfrow=c(2,2))
par(mar=c(3,4,3,2), oma=c(2,0,0,0))

plot(pt, probthreshold3[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability", ylab="Net benefit", 
     main="Outpatient setting",
     ylim=c(0, max(probthreshold3, na.rm=T)),
     xaxt='n')
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0,0.04, "A", font=2, cex=1.5)
lines(pt, probthreshold3[2,,1], col=mycolors[2])
lines(pt, probthreshold_sens_outpt[2,,1], col=mycolors[2], lty=2)
lines(pt, probthreshold_tat_outpt[1,,1], col=mycolors[1], lty=2)
legend(x = "topright", legend = c("NAAT (3 day)", "NAAT (2 day)", 
                                  "Ag-RDT (85% sens)","Ag-RDT (95% sens)"),
       col=c(mycolors[1],mycolors[1],mycolors[2],mycolors[2]), lty=c(1,2,1,2))

plot(pt, probthreshold3[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability", ylab="Net benefit", 
     main="Outpatient setting",
     ylim=c(0, max(probthreshold3, na.rm=T)),
     xaxt='n')
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0,0.04, "B", font=2, cex=1.5)
lines(pt, probthreshold3[2,,1], col=mycolors[2])
lines(pt, probthreshold3[3,,1], col=mycolors[3])
lines(pt, probthreshold3[5,,1], col='black', lty=1)
abline(h = 0, col='black', lty=4)
lines(pt, probthreshold4[3,,1], col=mycolors[3], lty=3)
lines(pt, probthreshold4[5,,1], col='black', lty=3)
legend(x = "topright", legend = c("NAAT (3 day)","Ag-RDT (85% sens)","Clinical, full", "Clinical, reduced",
                                  "Treat all, full", "Treat all, reduced", "Treat none"),
       col=c(mycolors[1],mycolors[2],mycolors[3], mycolors[3],'black','black','black'), lty=c(1,1,1,3,1,3,4))


plot(pt, probthreshold2[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability", ylab="Net benefit", 
     main="Hospital setting",
     ylim=c(0, max(probthreshold2, na.rm=T)),
     xaxt='n')
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0,0.29, "C", font=2, cex=1.5)
lines(pt, probthreshold2[2,,1], col=mycolors[2])
lines(pt, probthreshold_sens[2,,1], col=mycolors[2], lty=2)
lines(pt, probthreshold_tat[1,,1], col=mycolors[1], lty=2)
legend(x = "topright", legend = c("NAAT (1 day)", "NAAT (2 day)", 
                                  "Ag-RDT (85% sens)","Ag-RDT (95% sens)"),
       col=c(mycolors[1],mycolors[1],mycolors[2],mycolors[2]), lty=c(1,2,1,2))

plot(pt, probthreshold2[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability", ylab="Net benefit", 
     main="Hospital setting",
     ylim=c(0, max(probthreshold2, na.rm=T)),
     xaxt='n')
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0,0.29, "D", font=2, cex=1.5)
lines(pt, probthreshold2[2,,1], col=mycolors[2])
lines(pt, probthreshold2[3,,1], col=mycolors[3])
lines(pt, probthreshold2[5,,1], col='black', lty=1)
abline(h = 0, col='black', lty=4)
lines(pt, probthreshold[3,,1], col=mycolors[3], lty=3)
lines(pt, probthreshold[5,,1], col='black', lty=3)
legend(x = "topright", legend = c("NAAT (1 day)","Ag-RDT (85% sens)","Clinical, full", "Clinical, reduced",
                                  "Treat all, full", "Treat all, reduced", "Treat none"),
       col=c(mycolors[1],mycolors[2],mycolors[3], mycolors[3],'black','black','black'), lty=c(1,1,1,3,1,3,4))

mtext(side=1, outer=T, text = "Threshold probability above which early intervention for an average case provides net benefit", line=0.5, cex = 0.8)

pt
probthreshold[5,,1]
pt
probthreshold4[5,,1]


##### One-way sensitivity analysis, NAAT vs RDT. (Hosp setting) #####

gethilo <- function(params, level, hosp=T)
{
  hiloparams <- list()
  hiloparams$clinicalVsTransmissionBenefit <- (if(hosp) c(5,0.1) else c(0.1,0))
  hiloparams$prev <- (if(hosp) c(0.6,0.1) else c(0.2,0.02))
  hiloparams$falseDiagnosisHarm <- c(0.2,0.01)
  hiloparams$clinicalDiagnosisDiscount <- c(1,0.5)
  hiloparams$presx <- c(0.5,0.12)
  hiloparams$Period <- c(8,3)
  hiloparams$isolationOfKnownCase <- (if (hosp) c(1,0.5) else c(0.9,0.4))
  hiloparams$isolationOfCasePendingResult <- (if(hosp) c(0.9, 0.4) else c(0.6, 0))
  hiloparams$isolationOfContacts <- c(0.7,0)
  hiloparams$sxOnsetDay <- c(4,1)
  hiloparams$medianPresentationDay <- c(6,2)
  hiloparams$duration <- c(12,5)
  hiloparams$minimumInfectiousLogVirus <- c(5,2)
  hiloparams$infectivityScale <- c(3,0) # 1 = log-linear -- how skewed is the infectivity distribution, in relation to viral load?
  hiloparams$firstGenerationTransmissionWeight <- c(1,0) # how much do you value preventing transmission to immediate contacts,
  hiloparams$turnaroundTimeRDT <- c(1,0)
  hiloparams$turnaroundTimeNAAT <- (if(hosp) c(3,0.5) else c(7,1))
  hiloparams$dailyDecayClinicalBenefit <- c(0.3, 0.08)
  hiloparams$sensitivityClinician <- c(0.9, 0.7)
  hiloparams$specificityClinician <- c(0.7, 0.3)
  
  hiloparams$sensitivityRDT_vsNAAT <- c(0.95,0.75)
  hiloparams$specificityRDT <- c(0.995, 0.95)
  
  newparams <- params
   if(level=='high') 
   { for (name in names(hiloparams))
     newparams[[name]] <- hiloparams[[name]][1]
      newparams$weibpar <- get.weibull.par(p=c(0, hiloparams$presx[1], 0.99),
                                           q=c(0, hiloparams$sxOnsetDay[1],
                                               hiloparams$sxOnsetDay[1]+hiloparams$duration[1]), 
                                           plot = F, show.output = F)
  } else
  {for (name in names(hiloparams))
     newparams[[name]] <- hiloparams[[name]][2]
     newparams$weibpar <- get.weibull.par(p=c(0, hiloparams$presx[2], 0.99), 
                                        q=c(0, hiloparams$sxOnsetDay[2],
                                            hiloparams$sxOnsetDay[2]+hiloparams$duration[2]), 
                                        plot = F, show.output = F)
  }
  
  return(newparams)
}




####### sensitivty, transmission benefit ##########
hospplot <- T
resetparams <- make.params(hosp = hospplot)
highparams <- gethilo(params = resetparams, level="high", hosp = hospplot)
lowparams <- gethilo(params = resetparams, level="low", hosp = hospplot)
excludeparams <- c("weibpar", "normpar", "csens_h","csens_o","cspec_h","cspec_o")
infecthigh <- infectlow <- array(dim=c(3,length(resetparams)),dimnames = list(c("NAAT","RDT","clin"), names(resetparams)))
for (name in names(resetparams))
{ 
  onehighparam <- onelowparam <- resetparams; 
  onehighparam[[name]] <- highparams[[name]]
  onelowparam[[name]] <- lowparams[[name]]
  ref <- unlist(mapTransmissions(params = resetparams, plots=F))
  infecthigh[,name] <- unlist(mapTransmissions(params = onehighparam, plots = F))[1:3]
  infectlow[,name] <- unlist(mapTransmissions(params = onelowparam,  plots = F))[1:3]
  print(name)
  
}
infecthigh <- rbind(infecthigh, infecthigh["RDT",]/infecthigh["NAAT",])
rownames(infecthigh)[4] <- "assaydiff"
infectlow <- rbind(infectlow, infectlow["RDT",]/infectlow["NAAT",])
rownames(infectlow)[4] <- "assaydiff"
library(reshape2)
dathigh <- melt(infecthigh["assaydiff",]); dathigh$val <- "max"
datlow <- melt(infectlow["assaydiff",]); datlow$val <- "min"
datlow$names <- rownames(datlow); dathigh$names <- rownames(dathigh)
baseline <- ref["RDT"]/ref["NAAT"]
require(dplyr)
dat <- rbind(datlow, dathigh) 
keepnames <- setdiff(dathigh$names[pmax(dathigh$value/baseline, datlow$value/baseline, dathigh$value/datlow$value, datlow$value/dathigh$value, na.rm = T) > 1.1], 
                     excludeparams)
dat <- subset(dat, (dat$names %in% keepnames) )
dat <- dat %>% arrange(names)
# and include a version with extremes of multiple parameters
Nparams <- Rparams <- resetparams
for (name in keepnames) if (dat %>% filter(names==name, val=="max") %>% select(value) > 
                            dat %>% filter(names==name, val=="min") %>% select(value)) 
                        {Rparams[[name]] <- highparams[[name]]; Nparams[[name]] <- lowparams[[name]]} else 
                           {Rparams[[name]] <- lowparams[[name]]; Nparams[[name]] <- highparams[[name]]}
# dat <- rbind(c(unlist(mapTransmissions(params = Nparams,  plots = F))["RDT"]/
#                       unlist(mapTransmissions(params = Nparams,  plots = F))["NAAT"],
#                     "All favor alternative",
#                     "All of the above"),
#              dat)
# dat <- rbind(c(unlist(mapTransmissions(params = Rparams,  plots = F))["RDT"]/
#                       unlist(mapTransmissions(params = Rparams,  plots = F))["NAAT"],
#                     "All favor Ag-RDT",
#                     "All of the above"),
#              dat)

dat$value <- as.numeric(dat$value)
# prep for tornado plotting
class(dat) <- c("tornado", class(dat))
attr(dat, "output_name") <- "value"
dat$names <- paste0(niceParamNames[dat$names], "\n   (",unlist(lowparams)[dat$names], "-",unlist(highparams)[dat$names],")")
# dat$names[3:(nrow(dat))] <- paste0(dat$names[3:(nrow(dat))], " (",unlist(lowparams)[dat$names[3:(nrow(dat))]], "-",unlist(highparams)[dat$names[3:(nrow(dat))]],")")

library(ggplot2)
source("tornadoplot.R")
# label <- data.frame(x = 1, y = 2, label = c(">100"))
if(hospplot) {(gtrans_h <- ggplot_tornado(dat, baseline, title="One-way sensitivity analysis,\nTransmission benefit of Ag-RDT vs NAAT,\nhospital setting") + 
  ylab("Value of Ag-RDT relative to NAAT \n(values >1 favor Ag-RDT)") +
  xlab("Parameter (range explored)")+ 
    geom_hline(yintercept = 1)+
    ylim(c(0,1.6))); baseline_h <- baseline} else
{(gtrans_o <- ggplot_tornado(dat, baseline, title="One-way sensitivity analysis,\nTransmission benefit of Ag-RDT vs NAAT,\noutpatient setting") + 
  # scale_y_log10(limits=c(0.8,1.5), breaks=c(0.5,1,2,4,8,16,32)) + 
  # ylab("Value of Ag-RDT relative to NAAT \n(values >1 favor Ag-RDT)") + 
   ylab("") + 
  xlab("Parameter (range explored)")+
 geom_hline(yintercept = 1)+
   ylim(c(0,1.6))); baseline_o <- baseline}

# For text:
unlist(mapTransmissions(params = Nparams,  plots = F))["RDT"]/
  unlist(mapTransmissions(params = Nparams,   = F))["NAAT"]
unlist(mapTransmissions(params = Rparams,  plots = F))["RDT"]/
  unlist(mapTransmissions(params = Rparams,  plots = F))["NAAT"]


# gtrans_o / gtrans_h
gtrans_o + gtrans_h +   plot_layout(ncol = 1, heights = c(3,2))



###### sensitivity, net benefit  ##########
set.seed(55555)
hospplot <- T
resetparams <- make.params(hosp=hospplot)
resetparams$clinicalDiagnosisDiscount <- 1
highparams <- gethilo(params = resetparams, level="high", hosp = hospplot)
lowparams <- gethilo(params = resetparams, level="low", hosp=hospplot)
excludeparams <- c("weibpar", "normpar", "csens_h","csens_o","cspec_h","cspec_o")
excludeparams_clin <- c("weibpar", "normpar", "csens_h","csens_o","cspec_h","cspec_o", "isolationOfContacts")
NBhigh <- NBlow <- array(dim=c(5,length(resetparams)),dimnames = list(c("NAAT","RDT","clin","flex","all"), names(resetparams)))
for (name in names(resetparams))
{ 
  onehighparam <- onelowparam <- resetparams; 
  onehighparam[[name]] <- highparams[[name]]
  onelowparam[[name]] <- lowparams[[name]]
NBref <- NB(resetparams)
  NBhigh[,name] <- NB(params = onehighparam)
  NBlow[,name] <- NB(params = onelowparam)

}
NBhigh <- rbind(NBhigh, NBhigh["RDT",]/NBhigh["NAAT",], NBhigh["RDT",]/NBhigh["clin",])
NBlow <- rbind(NBlow, NBlow["RDT",]/NBlow["NAAT",], NBlow["RDT",]/NBlow["clin",])
rownames(NBhigh)[6] <- rownames(NBlow)[6] <- "assayratio"
rownames(NBhigh)[7] <- rownames(NBlow)[7] <- "clinratio"
library(reshape2)
dathigh <- melt(NBhigh["assayratio",]); dathigh$val <- "max"
dathighclin <- melt(NBhigh["clinratio",]); dathighclin$val <- "max"
datlow <- melt(NBlow["assayratio",]); datlow$val <- "min"
datlowclin <- melt(NBlow["clinratio",]); datlowclin$val <- "min"
datlow$names <- datlowclin$names <- rownames(datlow); dathigh$names <- dathighclin$names <- rownames(dathigh)
baseline <- NBref["RDT"]/NBref["NAAT"]; baselineclin <- NBref["RDT"]/NBref["clin"]

require(dplyr)
dat <- rbind(datlow, dathigh) 
datclin <- rbind(datlowclin, dathighclin) 
keepnames <- setdiff(dathigh$names[pmax(dathigh$value/baseline, datlow$value/baseline, dathigh$value/datlow$value, datlow$value/dathigh$value) > 1.1], excludeparams)
keepnamesclin <- setdiff(dathighclin$names[pmax(dathighclin$value/baselineclin, datlowclin$value/baselineclin, dathighclin$value/datlowclin$value, datlowclin$value/dathighclin$value) > 1.2], excludeparams_clin)
dat <- subset(dat, (dat$names %in% keepnames) ); datclin <- subset(datclin, (datclin$names %in% keepnamesclin) )
dat <- dat %>% arrange(names); datclin <- datclin %>% arrange(names)

# and include a version with extremes of multiple parameters
Nparams <- Rparams <- Nparamsclin <- Rparamsclin <- Nparamsclin2 <- Rparamsclin2  <- resetparams
for (name in keepnames)
   if (dat %>% filter(names==name, val=="max") %>% select(value) > 
                            dat %>% filter(names==name, val=="min") %>% select(value)) 
    {Rparams[[name]] <- highparams[[name]]; Nparams[[name]] <- lowparams[[name]]} else 
    {Rparams[[name]] <- lowparams[[name]]; Nparams[[name]] <- highparams[[name]]}
for (name in keepnamesclin)
  if (datclin %>% filter(names==name, val=="max") %>% select(value) > 
      datclin %>% filter(names==name, val=="min") %>% select(value)) 
  {Rparamsclin[[name]] <- highparams[[name]]; Nparamsclin[[name]] <- lowparams[[name]];
  Rparamsclin2[[name]] <- highparams[[name]]; Nparamsclin2[[name]] <- lowparams[[name]]} else 
  {Rparamsclin[[name]] <- lowparams[[name]]; Nparamsclin[[name]] <- highparams[[name]];
  Rparamsclin2[[name]] <- lowparams[[name]]; Nparamsclin2[[name]] <- highparams[[name]]}

# dat <- rbind(c(NB(params = Rparams)["RDT"]/ 
#                  NB(params = Rparams)["NAAT"],
#                "All favor Ag-RDT",
#                "All of the above"),
#              dat)
# dat <- rbind(c(NB(params = Nparams)["RDT"]/ 
#                  NB(params = Nparams)["NAAT"],
#                "All favor alternative",
#                "All of the above"),
#              dat)
dat$value <- as.numeric(dat$value)

# datclin <- rbind(c(NB(params = Rparamsclin)["RDT"]/ 
#                  NB(params = Rparamsclin)["clin"],
#                "All favor Ag-RDT",
#                "All of the above"),
#              datclin)
# datclin <- rbind(c(NB(params = Nparamsclin)["RDT"]/ 
#                  NB(params = Nparamsclin)["clin"],
#                "All favor alternative",
#                "All of the above"),
#              datclin)
datclin$value <- as.numeric(datclin$value)

# prep for tornado plotting
class(dat) <- c("tornado", class(dat)); class(datclin) <- c("tornado", class(datclin))
attr(dat, "output_name") <- "value"; attr(datclin, "output_name") <- "value"
dat$names <- paste0(niceParamNames[dat$names], "\n(",unlist(lowparams)[dat$names], "-",unlist(highparams)[dat$names],")")
datclin$names <- paste0(niceParamNames[datclin$names], "\n(",unlist(lowparams)[datclin$names], "-",unlist(highparams)[datclin$names],")")

library(ggplot2)
library(patchwork)
source("tornadoplot.R")
if(hospplot) (rdtnaat_h <- ggplot_tornado(dat, baseline, title="Net Benefit of Ag-RDT\nversus NAAT, hospital setting") +
                ylim(c(0,1.6)+
                       ylab("Value of Ag-RDT relative to NAAT \n(values >1 favor Ag-RDT)") +
                       xlab("Parameter (range explored)") )
    # scale_y_log10(limits=c(0.5,5), breaks = c(0.5,1,2,4,8))
      ) else
      (rdtnaat_o <- ggplot_tornado(dat, baseline, title="Net Benefit of Ag-RDT\nversus NAAT, outpatient setting") +
         ylim(c(0,1.6)+
                ylab("") +
                xlab("Parameter (range explored)") )
         # scale_y_log10(limits=c(0.5,max(dat$value)), breaks = c(0.5,1,2,4,8,16,32))
       ) 
if(hospplot) (rdtclin_h <- ggplot_tornado(datclin, baselineclin, title="Net Benefit of Ag-RDT\nversus clinical judgment, hospital setting") + 
                theme(legend.position = "none")+
                ylab("Value of Ag-RDT relative to clinical judgment \n(values >1 favor Ag-RDT)") +
                 xlab("Parameter (range explored)") 
              ) else
      (rdtclin_o <- ggplot_tornado(datclin, baselineclin, title="Net Benefit of Ag-RDT\nversus clinical judgment, outpatient setting") + theme(legend.position = "none") +
          ylab("") +
          xlab("Parameter (range explored)") 
          )
#to get outpatient setting, modify hospplot and repeat above.
 rdtnaat_o + rdtnaat_h +   plot_layout(ncol = 1, heights = c(3,2))
 rdtclin_o + rdtclin_h +   plot_layout(ncol = 1, heights = c(2,2))
 
if(hospplot) 
  { hi_benefit_clin_h <- NB(params = Nparamsclin)["RDT"]/ 
            NB(params = Nparamsclin)["clin"]
hi_benefit_h <- NB(params = Nparams)["RDT"]/ 
  NB(params = Nparams)["clin"] 
lo_benefit_clin_h <- NB(params = Rparamsclin)["RDT"]/ 
  NB(params = Rparamsclin)["clin"]
lo_benefit_h <- NB(params = Rparams)["RDT"]/ 
  NB(params = Rparams)["clin"] } else
  { hi_benefit_clin_o <- NB(params = Nparamsclin)["RDT"]/ 
    NB(params = Nparamsclin)["clin"]
  hi_benefit_o <- NB(params = Nparams)["RDT"]/ 
    NB(params = Nparams)["clin"] 
  lo_benefit_clin_o <- NB(params = Rparamsclin)["RDT"]/ 
    NB(params = Rparamsclin)["clin"]
  lo_benefit_o <- NB(params = Rparams)["RDT"]/ 
    NB(params = Rparams)["clin"] } 
 
# pdf(file = "DCA fig5.pdf", width = 6.5, height = 9)
# rdtnaat_o /
#   rdtnaat_h / 
#   rdtclin_h
# dev.off()

 
##### find precise values wher Ag-RDT preferred in hospital #####
resetparams <- make.params(hosp=T)
resetparams$turnaroundTimeNAAT <- 2
NB(resetparams)
resetparams <- make.params(hosp=T)
resetparams$sensitivityRDT_vsNAAT <- 0.85 **
NB(resetparams)
exp(-0.3*2.5)
 
 ##### explore other dimensions #######


resetparams <- make.params(); resetparams$clinicalDiagnosisDiscount <- 1
vary1 <- "sensitivityClinician"
vary2 <- "clinicalDiagnosisDiscount"
range1 <- seq(0.5,0.95,by=0.05)
range2 <- seq(0.5,1,by=0.05)
collateresults <- twowayvary(vary1, range1, vary2, range2)
par(mfrow=c(1,2), mar=c(4,4,3,1))
plot(range1, collateresults[1,,which.min(range2<resetparams[[vary2]])], type='l', col=mycolors[1], 
     main=paste0("Varying ",vary1), xlab=vary1, ylab="Net benefit", ylim=c(min(collateresults, na.rm=T), max(collateresults, na.rm=T)))
lines(range1, collateresults[2,,which.min(range2<resetparams[[vary2]])], col=mycolors[2])
lines(range1, collateresults[3,,which.min(range2<resetparams[[vary2]])], col=mycolors[3])
# lines(range1, collateresults[4,,which.min(range2<resetparams[[vary2]])], col='orange')
lines(range1, collateresults[5,,which.min(range2<resetparams[[vary2]])], col='black')
plot(range2, collateresults[1,which.min(range1<resetparams[[vary1]]),], type='l', col=mycolors[1], 
    main=paste0("Varying ",vary2), xlab=vary2, ylab="Net benefit", ylim=c(min(collateresults, na.rm=T), max(collateresults, na.rm=T)))
lines(range2, collateresults[2,which.min(range1<resetparams[[vary1]]),], col=mycolors[2])
lines(range2, collateresults[3,which.min(range1<resetparams[[vary1]]),], col=mycolors[3])
# lines(range2, collateresults[4,which.min(range1<resetparams[[vary1]]),], col='orange')
lines(range2, collateresults[5,which.min(range1<resetparams[[vary1]]),], col='black')
legend(x = "bottomright", legend = c("NAAT","Ag-RDT","Clinical", "Treat all"),
                     col=c(mycolors[1],mycolors[2],mycolors[3], 'black'), lty=1)
                     


# now varying turnaround time, clinical benefit decay
resetparams <- make.params(); resetparams$clinicalDiagnosisDiscount <- 1
vary1 <- "turnaroundTimeNAAT"
range1 <- seq(0.25,7,by=0.25)
collateresults <- twowayvary(p1 = vary1, r1 = range1, resetparams = resetparams)
par(mfrow=c(1,2), mar=c(4,4,3,1))
plot(range1, collateresults[1,,1], type='l', col=mycolors[1], 
     main="", xlab="NAAT turnaround time, days", ylab="Net benefit", ylim=c(min(collateresults, na.rm=T), max(collateresults, na.rm=T)))
lines(range1, collateresults[2,,1], col=mycolors[2])
lines(range1, collateresults[3,,1], col=mycolors[3])
lines(range1, collateresults[5,,1], col='black')

vary1 <- "sensitivityRDT_vsNAAT"
range1 <- seq(0.5,1,by=0.01)
collateresults2 <- twowayvary(p1 = vary1, r1 = range1, resetparams = resetparams)
plot(range1, collateresults2[1,,1], type='l', col=mycolors[1], 
     main="", xlab="Sensitivty of Ag-RDT", ylab="Net benefit", ylim=c(min(collateresults, na.rm=T), max(collateresults, na.rm=T)))
lines(range1, collateresults2[2,,1], col=mycolors[2])
lines(range1, collateresults2[3,,1], col=mycolors[3])
lines(range1, collateresults2[5,,1], col='black')
legend(x = "bottomright", legend = c("NAAT","Ag-RDT","Clinical", "Treat all"),
       col=c(mycolors[1],mycolors[2],mycolors[3], 'black'), lty=1)


# transform to data frame for 3d plotting
vary1 <- "turnaroundTimeNAAT"
range1 <- seq(0.25,5,by=1)
vary2 <- "sensitivityRDT_vsNAAT"
range2 <- seq(0.5,1,by=0.05)
collateresults <- twowayvary(vary1, range1, vary2, range2, resetparams = resetparams)
dimnames(collateresults) <- list(Assay = c("NAAT","RDT","clinical, fixed", "clinical, flexible", "treat all"),
                                 vary1 = range1,
                                 vary2 = range2)
require(reshape2)
dfdata <- melt(collateresults)
library(rgl)
rgl.open()
rgl.bg( sphere = F, fogtype = "none", color = "white",
        back = "lines")
plot3d(x = subset(dfdata, Assay=="NAAT")$vary1, y=subset(dfdata, Assay=="NAAT")$vary2, z=subset(dfdata, Assay=="NAAT")$value, type = 'l', col=mycolors[1],
       xlab=vary1, ylab=vary2, zlab="Net benefit", xlim=range(range1), ylim=range(range2), zlim=range(dfdata$value))
plot3d(x = subset(dfdata, Assay=="RDT")$vary1, y=subset(dfdata, Assay=="RDT")$vary2, z=subset(dfdata, Assay=="RDT")$value, type = 'l', col=mycolors[2], add = T)
plot3d(x = subset(dfdata, Assay=="clinical, fixed")$vary1, y=subset(dfdata, Assay=="clinical, fixed")$vary2, 
       z= subset(dfdata, Assay=="clinical, fixed")$value, type = 'l', col=mycolors[3], add = T)
legend3d(x="bottomright",legend = c("NAAT","Ag-RDT","clinical"), pch=16, col=c(mycolors[1],mycolors[2],mycolors[3]))
rgl.close()

library(plotly)
x <- list(title = vary1)
y <- list(title = vary2)
# fig1 <- plot_ly(x=range1, y=range2, z =collateresults[1,,], type='contour', contours=list(showlabels=T))
# fig1 <- fig1 %>% colorbar(title="Net benefit") %>% layout(xaxis=x, yaxis=y, title="NAAT")
# fig2 <- plot_ly(x=range1, y=range2, z =collateresults[2,,], type='contour', contours=list(showlabels=T))
# fig2 <- fig2 %>% colorbar(title="Net benefit") %>% layout(xaxis=x, yaxis=y, title="RDT")
# fig3 <- plot_ly(x=range1, y=range2, z =collateresults[3,,], type='contour', contours=list(showlabels=T))
# fig3 <- fig3 %>% colorbar(title="Net benefit") %>% layout(xaxis=x, yaxis=y, title="Clinical judgment")

# incremental NB contour plots:
#NAAT vs RDT
figa <- plot_ly(x=range2, y=range1, z =collateresults[1,,]- collateresults[2,,], type='contour', contours=list(showlabels=T))
figa <- figa %>% colorbar(title="Net benefit") %>% layout(xaxis=y, yaxis=x, title="NAAT versus RDT")
figa
# RDT vs clinical (flexible)
figb <- plot_ly(x=range2, y=range1, z =collateresults[1,,]- collateresults[4,,], type='contour', contours=list(showlabels=T))
figb <- figb %>% colorbar(title="Net benefit") %>% layout(xaxis=y, yaxis=x, title="NAAT versus clinical")
figb






#   
#   [How to account for error/imprecision in the clinician's calibration?]
# 
# [Or for RDT, account for there being some the clinician would judge so probable that they wouldn't do the RDT??]
#   
#   [For next piece: Perspective with multiple issues in decision curves (wo too much depth on any one), OR methods piece just on the baseline of DCA that is better than treat all/treat none and that some get treated regardless]
#   
  
#   
# #example for a specific threshold:
# tprob <- 0.2
# i <- which(x==tprob)
#   cutoffNBs <- numeric(length(scores$score))
#   for (j in 1:length(cutoffNBs))
#   {
#     cutoffNBs[j] <- NBdiscrete(scores$score[j], x[i])
#   }
#   scores$prevhalfcutoffs <- cutoffNBs
#   library(ggplot2)
#   ggplot(data=scores, aes(x=score, fill=factor(covid))) + 
#     geom_histogram(alpha=0.5, position="identity", binwidth = 0.1) + 
#     geom_line(aes(x=score, y=prevhalfcutoffs))
#   
#   
#   choosescore <- numeric(length(x))
#   considercutoffs <- seq(floor(min(scores$score)), ceiling(max(scores$score)), by = 0.1)
#   for (i in 1:length(x))
#   {
#     cutoffNBs <- numeric(length(considercutoffs))
#     for (j in 1:length(cutoffNBs))
#     {
#       cutoffNBs[j] <- NBdiscrete(considercutoffs[j], x[i])
#     }
#     choosescore[i] <- considercutoffs[which.max(cutoffNBs)] # account for uncertainty here somewhere?
#   }
#   
#   # what score cutoff should be used for each threshold probability, to maximize net benefit?
#   plot(x, choosescore)
#   NBs <- numeric(length(x))
#   for (i in 1:length(NBs)) NBs[i] <- NBdiscrete(cutoff = choosescore[i], thresholdprob = x[i])
#   plot(x, NBs)
#   
#   
#   plotDCA(prev=prev)
#   plotDCA(prev=prev, plotclin = F)
#   lines(x, NBs/nsim, col=mycolors[3])


###### Model of clinician judgment #########

# If a clinician knows what the above plot looks like, they can choose the point which 
# has the PPV that they want.
## -- or, for more complicated NB model, can choose the cutoff that maximizes NB.

# point sensitivityClinician and specificityClinician -> binormal model of ROC curve -> 
# For each cutoff and prev, can calcuate number of TPs and number of FPs ->
# within NB model, can choose the cutoff that maximizes NB for a given set of params, and 
# use the corresponding sensitivityClinician and specificityClinician in the overall NB calculations.


# clinician_options <- function(params, histplot=F, cdfplot=F, s = 100)
# {
#   # Assume clinicians evaluate patients based on clinicial features (incl symptoms, exposure risks, and potential for alternative diagnoses) and could hypothetically rank patients on some arbitary "liklihood of COVID" scale. This ranking may depend on individual-level exposure risk but will be independent of the general prevalence of COVID in the population. In keeping with common approaches to ROC analysis, assume that patients with COVID have clinical-judgment liklihoods normally distributed as $N(m_p,\inc_p)$, and patients without COVID follow another normal distribution $N(m_n, \inc_n)$. 
#   m_p <- 0 # arbitrary mean likelihood of positives
#   sd_p <- 1 # also arbitrary
#   sensquantile <- qnorm(p = params$sensitivityClinician, mean = 0, sd=1, lower.tail = F)
#   
#   # so we want sensquantile of the first distribution, to equal specquantile of the second (lowertail)
#   sd_n <- sd_p # assuming equal for now, but can adjust.
#   # Now find the mean for the negatives:  
#   # m_n + sd_n * qnorm(p = 0.5, mean=0, sd=1) = sensquantile
#   m_n = sensquantile - sd_n * qnorm(p = params$specificityClinician, mean=0, sd=1)
#   
#   if(histplot) {
#     # Plot the distributions, at two levels of prevalence:
#     plothist <- function(prev, nsim, title)
#     { L_p <- rnorm(n = nsim*prev, mean = m_p, sd = sd_p)
#     L_n <- rnorm(n = nsim*(1-prev), mean = m_n, sd = sd_n)
#     hist(L_n, freq = T, breaks=seq(-6,6,by=0.1), main=title, xlab="arbitrary ranking scale")
#     hist(L_p, freq = T, add=T, col=mycolors[3], breaks=seq(-6,6,by=0.1) )
#     }
#     par(mfrow=c(1,2))
#     plothist(prev = 0.1, nsim = 100000, title="10% prevalence")
#     plothist(prev = 0.4, nsim = 100000, title="40% prevalence")
#   }
#   prev <- 0.1
#   nsim <- 100000
#   L_p <- rnorm(n = nsim*prev, mean = m_p, sd = sd_p)
#   L_n <- rnorm(n = nsim*(1-prev), mean = m_n, sd = sd_n)
#   if(cdfplot){
#     # cdf plot illustrates overlap  
#     par(mfrow=c(1,1))
#     plot(ecdf(x = L_p), xlab="Diagnostic threshold on arbitrary scale", ylab="Proportion classified as COVID-negative")
#     lines(ecdf(x = L_n), col=mycolors[3])
#     abline(v=0, lty=2)
#     # abline(v=-1, lty=3, col='gray')
#     # abline(v=0.5, lty=3, col='gray')
#     legend("bottomright", c("Negatives","Positives","Threshold used in fixed-cutoff model"), lty=c(1,1,2), col=c(mycolors[1],"black","black"))
#   }  
#   
#   cutoffs <- seq(floor(min(L_n)*s)/s,ceiling(max(L_p)*s)/s, by=1/s)
#   sns <- sps <- numeric(length(cutoffs))
#   for (i in 1:length(cutoffs))
#   {
#     sns[i] <- mean(L_p>cutoffs[i])
#     sps[i] <- mean(L_n<cutoffs[i])
#   }
#   
#   return(rbind(sns, sps))
# }
# 
# clinician_options(make.params(), s=10, histplot = T, cdfplot = T)  
# 



