library(rriskDistributions)
library(ReIns)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(ggplot2)


set.seed(12345)
mycolors <- brewer.pal(n = 9, name = 'Set1')

#### inputs #####
make.params <- function(hosp=T)
{
  params <- list()
  params <- within(params, {
    #weights
    I_T <- 1 # value of preventing transmission from a case
    clinicalVsTransmissionBenefit <- ifelse(hosp, 2, 0.06) # value of preventing mortality/morbidity, in units of T
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
    firstGenerationTransmissionWeight <- 0.2 # how much do you value preventing transmission to immediate contacts,

    # timing
    medianPresentationDay <- ifelse(hosp, 5, 3.5)
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

mapTransmissions <- function(params,
                            npts = 1e6, #arbitrary sim size
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
  
  # if(plots)
  # {
    # print(
    #   ggplot(data.frame(logvirus)) + geom_density(aes(x=logvirus), lwd=0.5) + 
    #   geom_vline(xintercept=min(logvirus[Pos_NAAT_onset==1]), lty=2, color=mycolors[1], lwd=0.5) + 
    #   geom_vline(xintercept=min(logvirus[Pos_RDT_onset==1]), lty=3, color=mycolors[2], lwd=0.5) + 
    #   xlab("Log10 virus at symptom onset") +
    #   ylab("Proportion of patient population") +
    #   xlim(1,10)+
    #   theme_minimal(base_size = 12) +
    #   scale_y_continuous(breaks = NULL) + 
    #   coord_flip() + 
    #   annotate(geom="text", x=min(logvirus[Pos_NAAT_onset==1])+0.1, y=0, 
    #            label="Limit of NAAT detection", color=mycolors[1], hjust=0, vjust=0, size=4) + 
    #   annotate(geom="text", x=min(logvirus[Pos_RDT_onset==1])+0.1, y=0, 
    #            label="Limit of Ag-RDT detection", color=mycolors[2], hjust=0, vjust=0, size=4) 
    # )
  # }

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
    
    ratio <- sample(infectivity[Pos_RDT_onset==1], size = 100000, replace = T)/ 
        sample(infectivity[Pos_NAAT_onset & !Pos_RDT_onset], size=100000, replace=T) 
      
    print(mean(ratio)); print(sd(ratio)  )
    
    # plot of infectivity versus peak log viral burden
    v <- seq(0,12,by=0.1)
    i <- (pmax(0, pmin(v-minimumInfectiousLogVirus, max(v)-minimumInfectiousLogVirus)/(max(v)-minimumInfectiousLogVirus)))^infectivityScale
    i <- i/mean(i) # this is infectivity relative to average (at any time point, assumes remains proportional)
    par(mar=c(4,4,1,1))
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
                                                               size=1, prob = isolationOfKnownCase*caseevents$Pos_clin*clinicalDiagnosisDiscount)
  caseevents$Averted_treatall <- caseevents$Avertable_clin * rbinom(n=nrow(caseevents), # incorporating test sensitivity and isolation effectiveness
                                                                size=1, prob = isolationOfKnownCase*clinicalDiagnosisDiscount)
  
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
  caseevents$Transmits_treatall <- 1 - caseevents$Averted_treatall
  caseevents$contactTransmissionAvertable_clin <- caseevents$Transmits_clin & caseevents$contacttime >presentationDay + 0
  caseevents$contactTransmissionAverted_clin <- caseevents$contactTransmissionAvertable_clin * rbinom(n=nrow(caseevents), 
                                                                                                      size=1, prob = isolationOfContacts*Pos_clin*clinicalDiagnosisDiscount)
  caseevents$contactTransmissionAverted_treatall <- caseevents$contactTransmissionAvertable_clin * rbinom(n=nrow(caseevents), 
                                                                                                      size=1, prob = isolationOfContacts*clinicalDiagnosisDiscount)
  caseevents$contactTransmits_clin <- caseevents$Transmits_clin & !caseevents$contactTransmissionAverted_clin
  caseevents$contactTransmits_treatall <- caseevents$Transmits_treatall & !caseevents$contactTransmissionAverted_treatall
  
  # casetransmits are timing, caseTransmissions are indexes fo all events (some repeated), 
   # caseTransmits_X are indexed for events that occur with a given test (subset of caseTransmissions)
  
  caseevents$Transmission <- ifelse(caseevents$Transmits_RDT==0, ifelse(caseevents$Transmits_NAAT==0, "Prevented with either test", "Prevented with Ag-RDT only"), 
                                 ifelse(caseevents$Transmits_NAAT==0, "Prevented with NAAT only", "Occurs regardless of test"))
  caseevents %>% dplyr::count(Transmits_NAAT, Transmits_RDT, Transmission)
  caseevents$Transmission <- factor(caseevents$Transmission, levels=c("Occurs regardless of test",
    "Prevented with NAAT only",
    "Prevented with either test",
    "Prevented with Ag-RDT only"))

    
    # by timing of transmission, cases
    arranged <- caseevents %>% 
      arrange(value) %>% 
      dplyr::mutate(rn = row_number())
    arranged_RDT <- subset(caseevents, Transmits_RDT==1) %>% 
      arrange(value) %>% 
      dplyr::mutate(rn = row_number())
    arranged_NAAT <- subset(caseevents, Transmits_NAAT==1) %>% 
      arrange(value) %>% 
      dplyr::mutate(rn = row_number())
    arranged_clin <- subset(caseevents, Transmits_clin==1) %>% 
      arrange(value) %>% 
      dplyr::mutate(rn = row_number())
    
    # by timing of transmission, contacts
    gen2 <- caseevents %>% 
      arrange(contacttime) %>% 
      dplyr::mutate(rn = row_number())
    gen2_RDT <- subset(caseevents, contactTransmits_RDT==1) %>% 
      arrange(contacttime) %>% 
      dplyr::mutate(rn = row_number())
    gen2_NAAT <- subset(caseevents, contactTransmits_NAAT==1) %>% 
      arrange(contacttime) %>% 
      dplyr::mutate(rn = row_number())
    gen2_clin <- subset(caseevents, contactTransmits_clin==1) %>% 
      arrange(contacttime) %>% 
      dplyr::mutate(rn = row_number())
    
    # by timing of presentation, for all cases not weighted as source cases**
    allcases <- data.frame(presentationDay)
    allcases$RDTDay <- presentationDay + turnaroundTimeRDT
    allcases$NAATDay <- presentationDay + turnaroundTimeNAAT
    allcases$Pos_RDT <- Pos_RDT
    allcases$Pos_NAAT <- Pos_NAAT
    allcases$Pos_clin <- Pos_clin
    tested <- allcases %>% 
      arrange(presentationDay) %>% 
      dplyr::mutate(rn = row_number())
    detected_RDT <- subset(allcases, Pos_RDT==1) %>% 
      arrange(RDTDay) %>% 
      dplyr::mutate(rn = row_number())
    detected_NAAT <- subset(allcases, Pos_NAAT==1) %>% 
      arrange(NAATDay) %>% 
      dplyr::mutate(rn = row_number())
    detected_clin <- subset(allcases, Pos_clin==1) %>% 
      arrange(presentationDay) %>% 
      dplyr::mutate(rn = row_number())
  
    
    if(plots){
      
      colors <- c(
      # "Any/All" = mycolors[5],
                "NAAT" = mycolors[1], "Ag-RDT" = mycolors[2], 
                "Clinical" = mycolors[3])
    linetypes <- c("All results" = "dotted", 
                "Positive results" = 'solid')
    
    xmax <- 25
    (detection <- ggplot() +
      xlim(-3,xmax) + 
        geom_step(data=detected_clin[seq(1,nrow(detected_clin), by=npts/1000),], aes(x=presentationDay, y=rn, color='Clinical', linetype='Positive results'),lwd=0.8) +
        geom_step(data=tested[seq(1,nrow(tested), by=npts/1000),], aes(x=presentationDay, y=rn, color="Clinical", linetype="All results"),lwd=0.8) +
        geom_step(data=tested[seq(1,nrow(tested), by=npts/1000),], aes(x=presentationDay+turnaroundTimeRDT, y=rn, color="Ag-RDT", linetype="All results"),lwd=0.8) +
      geom_step(data=tested[seq(1,nrow(tested), by=npts/1000),], aes(x=presentationDay+turnaroundTimeNAAT, y=rn, color="NAAT", linetype="All results"),lwd=0.8) +
      geom_step(data=detected_RDT[seq(1,nrow(detected_RDT), by=npts/1000),], aes(x=presentationDay+turnaroundTimeRDT, y=rn, color="Ag-RDT", linetype="Positive results"),lwd=0.8) +
      geom_step(data=detected_NAAT[seq(1,nrow(detected_NAAT), by=npts/1000),], aes(x=presentationDay+turnaroundTimeNAAT, y=rn, color='NAAT', linetype="Positive results"),lwd=0.8)  +
      labs(x = "Days from symptom onset",
           color = "Diagnostic approach", linetype="Cumulative results\nfrom patients w/COVID") +
        scale_color_manual(values = c(mycolors[2], mycolors[1], mycolors[3]), breaks = c("Ag-RDT", "NAAT", "Clinical")) + 
      scale_linetype_manual(values = c("dotted", "solid"), breaks = c("All results", "Positive results")) + 
      scale_y_continuous(name = "Cumulative test results received,\namong 1,000,000 simulated patients with COVID-19", 
                         breaks = npts*c(0,0.25,0.5,0.75,1), 
                         labels = c("0","250k","500k","750k","1000k")) + 
        theme_bw()+
        theme(legend.key.width=unit(1,"cm"))+ 
        theme(legend.position = c(0.75, 0.3), 
              legend.title = element_text(size = 10), legend.text = element_text(size = 9))
    )
    print(
      detection + geom_label(aes(x=c(15,xmax,22), 
                                 y=c(max(detected_RDT$rn), max(detected_clin$rn)+npts/20, max(detected_NAAT$rn)-npts/20)),
                    label = paste0(round(c(max(detected_RDT$rn), max(detected_clin$rn), max(detected_NAAT$rn))/(npts)*100,1),"%"),
                      label.size=0.2, 
                    col=colors[c(2,3,1)])
    )
    
    
    # # ## for Overview figure, timing of transmission:
    # # 
    # par(mfrow=c(1,1), mar=c(3,2,1,1), oma=c(0,0,0,0))
    # plot(NA, xlab="Days since symptom onset",
    #          xlim=c(-3,25), yaxt='n', ylab="",
    #      ylim=c(0,4.5e3))
    # rect(xleft = 3, xright = 30, ybottom = 0, ytop = 4.5e4, col = rgb(0.5,0.5,0.5,1/4), border = NA)
    # rect(xleft = 6, xright = 30, ybottom = 0, ytop = 4.5e4, col = rgb(0.5,0.5,0.5,1/4), border=NA)
    # lines(names(allcase), allcase, col=mycolors[5], lwd=2)
    # lines(names(allcontact), 0.25*allcontact, col=mycolors[4], lwd=2)
    # text(5,4.2e4, "From index cases", col=mycolors[5])
    # text(16,2.75e4, "From direct contacts", col=mycolors[4])
    # segments(x0 = 0, x1=0, y0=0, y1 = allcase['0'], lty=2, lwd=1.5)
    # text(-0.8,500, "Symptom onset", srt=90, adj=0)
    # segments(x0 = 3, x1=3, y0=0, y1 = allcase['3'], lty=3, lwd=1.5)
    # # text(2.2,500, "Early diagnosis", srt=90, adj=0)
    # segments(x0 = 6, x1=6, y0=0, y1 = allcase['6'], lty=3, lwd=1.5)
    # # text(5.2,500, "Later diagnosis", srt=90, adj = 0)
    # mtext("Days since symptom onset", side = 1, line=2)
  }
  
  transmission_benefit_RDT <- 1 - 
   ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_RDT) + 
    (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_RDT) )
  
  transmission_benefit_NAAT <- 1 - 
    ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_NAAT) + 
        (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_NAAT) )
  
  transmission_benefit_clin <- 1 - 
    ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_clin) + 
        (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_clin) ) 
  
  transmission_benefit_treatall <- 1 - 
    ( firstGenerationTransmissionWeight*mean(caseevents$Transmits_treatall) + 
        (1-firstGenerationTransmissionWeight)*mean(caseevents$contactTransmits_treatall) ) 
  
  # clinical benefit, exponential decay
  clinical_benefit_RDT <- mean(allcases$Pos_RDT * exp(-(allcases$presentationDay + turnaroundTimeRDT) * dailyDecayClinicalBenefit))
  clinical_benefit_NAAT <- mean(allcases$Pos_NAAT * exp(-(allcases$presentationDay + turnaroundTimeNAAT) * dailyDecayClinicalBenefit))
  clinical_benefit_clin <- mean(allcases$Pos_clin * exp(-(allcases$presentationDay) * dailyDecayClinicalBenefit)) *clinicalDiagnosisDiscount
  clinical_benefit_treatall <- mean(exp(-(allcases$presentationDay) * dailyDecayClinicalBenefit)) *clinicalDiagnosisDiscount
  
  return(list(# transmission reduction after g1 weighting:
              "NAAT"=transmission_benefit_NAAT, "RDT"=transmission_benefit_RDT, "clin"=transmission_benefit_clin,
                "treatall"=transmission_benefit_treatall,
              # reduction in avertable M&M 
              "NAATclinical"=clinical_benefit_NAAT, "RDTclinical"=clinical_benefit_RDT, "clinclinical"=clinical_benefit_clin,
                "treatallclinical"=clinical_benefit_treatall,
              # proportion of cases detected
              "NAATdetection"=mean(Pos_NAAT), "RDTdetection"=mean(Pos_RDT), "clindetection"=mean(Pos_clin),
              # proportion of transmission generated by those who are detected
              "NAATpeakinf" = sum(infectivity[Pos_NAAT==1])/sum(infectivity),
              "RDTpeakinf" = sum(infectivity[Pos_RDT==1])/sum(infectivity), 
              "clinpeakinf" = sum(infectivity[Pos_clin==1])/sum(infectivity),
              # proportion of future transmission (after clinical presentation) generated by those who are detected
              "NAATfuturetrans" = sum(subset(caseevents, Pos_NAAT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "RDTfuturetrans" = sum(subset(caseevents, Pos_RDT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
              "clinfuturetrans" = sum(subset(caseevents, Pos_clin==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
    #               "NAATavertable" = sum(subset(caseevents, Avertable_NAAT==1 & Pos_NAAT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
    #               "RDTavertable" = sum(subset(caseevents, Avertable_RDT==1 & Pos_RDT==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
    #               "clinavertable" = sum(subset(caseevents, Avertable_clin==1 & Pos_clin==1)$Avertable_clin)/sum(caseevents$Avertable_clin),
    #               "NAATcontactavertable" = sum(subset(caseevents, contactTransmissionAvertable_NAAT==1 & Pos_NAAT==1)$contactTransmissionAvertable_clin)/sum(caseevents$contactTransmissionAvertable_clin),
    #               "RDTcontactavertable" = sum(subset(caseevents, contactTransmissionAvertable_RDT==1 & Pos_RDT==1)$contactTransmissionAvertable_clin)/sum(caseevents$contactTransmissionAvertable_clin),
    #               "clincontactavertable" = sum(subset(caseevents, contactTransmissionAvertable_clin==1 & Pos_clin==1)$contactTransmissionAvertable_clin)/sum(caseevents$contactTransmissionAvertable_clin),
              # allavertable: of all (case) transmission events, which occurred after a positive test resulted
              "allNAATavertable" = sum(subset(caseevents, Avertable_NAAT==1 & Pos_NAAT==1)$Avertable_clin)/nrow(caseevents),
              "allRDTavertable" = sum(subset(caseevents, Avertable_RDT==1 & Pos_RDT==1)$Avertable_clin)/nrow(caseevents),
              "allclinavertable" = sum(subset(caseevents, Avertable_clin==1 & Pos_clin==1)$Avertable_clin)/nrow(caseevents),
              # of all (contact) transmission events, which occurred after a positive test resulted
              "allNAATcontactavertable" = sum(subset(caseevents, contactTransmissionAvertable_NAAT==1 & Pos_NAAT==1)$contactTransmissionAvertable_clin)/nrow(caseevents),
              "allRDTcontactavertable" = sum(subset(caseevents, contactTransmissionAvertable_RDT==1 & Pos_RDT==1)$contactTransmissionAvertable_clin)/nrow(caseevents),
              "allclincontactavertable" = sum(subset(caseevents, contactTransmissionAvertable_clin==1 & Pos_clin==1)$contactTransmissionAvertable_clin)/nrow(caseevents), 
              # split out T1j and T2j separately:    
              "casetransRDT" = mean(caseevents$Transmits_RDT), 
              "casetransNAAT" = mean(caseevents$Transmits_NAAT),
              "contacttransRDT" = mean(caseevents$contactTransmits_RDT),
              "contacttransNAAT" = mean(caseevents$contactTransmits_NAAT)
              
              
              ) ) 
  })
  
}


(mo <- mapTransmissions(params = make.params(hosp=F), plots=T, npts=1e6))
(mh <- mapTransmissions(params = make.params(hosp=T), plots=T, npts=1e6))

mo$RDTdetection
mo$RDTfuturetrans

mo$RDT
mo$NAAT
mo$clin

# # Supplemental results, individual components
# (1-mh$casetransRDT)*1000
# (1-mo$casetransRDT)*1000
# (1-mh$casetransNAAT)*1000
# (1-mo$casetransNAAT)*1000
# (1-params$firstGenerationTransmissionWeight)/params$firstGenerationTransmissionWeight *(1-mh$contacttransRDT)*1000
# (1-params$firstGenerationTransmissionWeight)/params$firstGenerationTransmissionWeight *(1-mo$contacttransRDT)*1000
# (1-params$firstGenerationTransmissionWeight)/params$firstGenerationTransmissionWeight *(1-mh$contacttransNAAT)*1000
# (1-params$firstGenerationTransmissionWeight)/params$firstGenerationTransmissionWeight *(1-mo$contacttransNAAT)*1000
# 
# (1-0.1)*(1-spec$NAAT) * 1000
# (1-0.1)*(1-spec$RDT) * 1000
# (1-0.4)*(1-spec$NAAT) * 1000
# (1-0.4)*(1-spec$RDT) * 1000
# 
# 1000 * mh$NAATclinical * 0.04
# 1000 * mo$NAATclinical * 0.04/30
# 1000 * mh$RDTclinical * 0.04
# 1000 * mo$RDTclinical * 0.04/30
# 
# # stochastic variability and sample size
# size <- 1e6
# o1 <- o2 <- o3 <- o4 <- numeric(10)
# for (i in 1:10)
#   o1[i] <- mapTransmissions(params = make.params(hosp=T), plots=F, plotclin = F, npts=size)$RDT
# summary(o1)
# for (i in 1:10)
#   o2[i] <- mapTransmissions(params = make.params(hosp=T), plots=F, plotclin = F, npts=size)$NAAT
# summary(o2)
# for (i in 1:10)
#   o3[i] <- mapTransmissions(params = make.params(hosp=T), plots=F, plotclin = F, npts=size)$RDTclinical
# summary(o3)
# for (i in 1:10)
#   o4[i] <- mapTransmissions(params = make.params(hosp=T), plots=F, plotclin = F, npts=size)$NAATclinical
# summary(o4)


##### Calculate Net Benefit with modified formula account for time-sensitivty, diagnostic-associated infectivity  ######

# Net benefit, as a function of test and setting/prevalence and other input parameters incl weights
# Will want to display this across multiple dimensions of variation, perhaps up to two at a time

# get the net benefit for a single set of parameters
# the value calculated is the average across the cases that are detected, or the average across all cases? Need to modify to be those detected. 
NB <- function(params, npts=1e6)
{
  #need to base this on the proportion actually detected, after accounting for delays. mapTransmissions will account for the lack of trasmission benefit, but we also need to track who's detected at all and gets clinical benefit.
  
    with(params, 
       {   
         mt <- mapTransmissions(params, plots = F, plotclin = F, npts = npts)
         
         trans <- mt[c("NAAT","RDT","clin","treatall")] # transmission benefit, proportion of all transmission prevented, weighted by g1 and 1-g1
         clinical <- mt[c("NAATclinical" ,"RDTclinical", "clinclinical","treatall")] # mean of (detected)*exponential decay, not yet weighted by c
        sens <- list(NAAT = mt$NAATdetection, RDT=mt$RDTdetection, clin=mt$clindetection, treatall=1) # proportion detected at time of clinical presentation
        spec <- list(NAAT=spec_NAAT, RDT = specificityRDT, clin=specificityClinician, treatall=0)
      
        # # already accounted for clinical diagnosis discount effects on both prevented transmission and clinical benefits, in maptransmission

        # now combine all sources of TP benefit, and discount action in response to clinical decision
        # 
        V_TP <- (mapply(function(x,w,z,v) I_T *( # normalizing constant
                          ( x + clinicalVsTransmissionBenefit*w)  / 
                          ((1-presx)*(isolationOfKnownCase)) - 
                          z * falseDiagnosisHarm* v), # normalizing so that one full TP get a 1 on both NB scales, and so that we're in units of avertible transmission
                      x=trans, w = clinical, z = sens, v=c(1,1,clinicalDiagnosisDiscount, clinicalDiagnosisDiscount))) #already accounted for clinical diagnosis discount in clin detection
        
        V_FP <-  I_T * falseDiagnosisHarm * c(1,1,clinicalDiagnosisDiscount, clinicalDiagnosisDiscount) * (1-unlist(spec))
        #( (-c_F*clinicalVsTransmissionBenefit ) - (n*isolationOfContacts+1)*i_F)
        # removed factor (1+clinicalVsTransmissionBenefit), to normalize to NB_simple
        
        NB <- prev*V_TP - (1-prev)*V_FP
        
      return(NB)
  })
}

NB(params = make.params(hosp=F))
NB(params = make.params(hosp=T))

twowayvary <- function(p1, r1, p2=NA, r2=NA, 
                       resetparams,
                       npts=1e5# change if needed, for any that will be fixed at a nondefault value
)
{
  # decide which parameters will vary and over what range
  if(missing(resetparams)) resetparams <- make.params() 
  vary1 <- p1
  range1 <- r1
  vary2 <- p2
  range2 <- r2 
  collateresults <- array(NA, dim=c(4, length(range1), ifelse(is.na(vary2), 1, length(range2))))
  if(!is.na(vary1)) 
    for (set1 in 1:length(range1))
    {
      params <- resetparams
      params[[vary1]] <- range1[set1]
      if(!is.na(vary2)) {
        for (set2 in 1:length(range2))
        {  params[[vary2]] <- range2[set2] 
        collateresults[,set1,set2] <- NB(params, npts)
        }} else collateresults[,set1,1] <- NB(params, npts)
    }
  return(collateresults)
}


###### Decision curves for threshold probability #######

vary1 <- "falseDiagnosisHarm"
# if we want probability threshold pt to range from 0 to 1, then 
# pt = falseDiagnosisHarm*(1-pt)
pt <- seq(0.001,0.991,by=0.03)
range1 <- pt#/(1-pt)

#hospital
probthreshold <- twowayvary(vary1, range1) # standard param values
altparams <- make.params(hosp=T); altparams$clinicalDiagnosisDiscount <- 1
probthreshold2 <- twowayvary(vary1, range1, resetparams = altparams)
altparams_sens <- altparams_tat <- altparams_tat3 <- altparams
altparams_sens$sensitivityRDT_vsNAAT <- 0.95 # for acute infection, 90% in Kruger, 95% in Igloi
altparams_tat$turnaroundTimeNAAT <- 2
altparams_tat3$turnaroundTimeNAAT <- 3
probthreshold_sens <- twowayvary(vary1, range1, resetparams = altparams_sens) # standard param values
probthreshold_tat <- twowayvary(vary1, range1, resetparams = altparams_tat) 
probthreshold_tat3 <- twowayvary(vary1, range1, resetparams = altparams_tat3)
#outpatient
altparams <- make.params(hosp = F); 
probthreshold4 <- twowayvary(vary1, range1, resetparams = altparams)
altparams$clinicalDiagnosisDiscount <- 1
probthreshold3 <- twowayvary(vary1, range1, resetparams = altparams)
altparams_sens_outpt <- altparams_tat_outpt <- altparams_tat1_outpt <- altparams
altparams_sens_outpt$sensitivityRDT_vsNAAT <- 0.95 # for acute infection, 90% in Kruger, 95% in Igloi
altparams_tat_outpt$turnaroundTimeNAAT <- 2 # for acute infection, 90% in Kruger, 95% in Igloi
altparams_tat1_outpt$turnaroundTimeNAAT <- 1
probthreshold_sens_outpt <- twowayvary(vary1, range1, resetparams = altparams_sens_outpt) # standard param values
probthreshold_tat_outpt <- twowayvary(vary1, range1, resetparams = altparams_tat_outpt) # standard param values
probthreshold_tat1_outpt <- twowayvary(vary1, range1, resetparams = altparams_tat1_outpt) # standard param values

# Figure 3
layout(t(array(c(1,1,1,2,2,2,2), dim=c(1,7))))
par(mar=c(3,4,3,2), oma=c(2,0,0,0))

plot(pt, probthreshold3[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability for intervention", ylab="Net benefit", 
     main="Outpatient setting",
     ylim=c(0, max(probthreshold3, na.rm=T)*0.9),
     xaxt='n', lty=5, lwd=2)
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0,0.08, "A", font=2, cex=1.5)
lines(pt, probthreshold3[2,,1], col=mycolors[2], lwd=2)
lines(pt, probthreshold_sens_outpt[2,,1], col=mycolors[2], lty=4, lwd=1)
lines(pt, probthreshold_tat_outpt[1,,1], col=mycolors[1], lty=3, lwd=1)
lines(pt, probthreshold_tat1_outpt[1,,1], col=mycolors[1], lty=1, lwd=1)
legend(x = "bottomleft", legend = c("NAAT (1 day)", "NAAT (2 day)", "NAAT (3 day)", 
                                  "Ag-RDT (95% acute sens)","Ag-RDT (85% acute sens)"),
       col=c(mycolors[1],mycolors[1],mycolors[1],mycolors[2],mycolors[2]), 
       lty=c(1,3,5,4,1), 
       lwd=c(1,1,2,1,2), seg.len = 3)

plot(pt, probthreshold2[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability for intervention", ylab="Net benefit", 
     main="Hospital setting",
     ylim=c(0, max(probthreshold2, na.rm=T)),
     xaxt='n', lwd=2)
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0,0.195, "B", font=2, cex=1.5)
lines(pt, probthreshold2[2,,1], col=mycolors[2], lwd=2)
lines(pt, probthreshold_sens[2,,1], col=mycolors[2], lty=4, lwd=1)
lines(pt, probthreshold_tat[1,,1], col=mycolors[1], lty=3, lwd=1)
lines(pt, probthreshold_tat3[1,,1], col=mycolors[1], lty=5, lwd=1)
legend(x = "bottomleft", legend = c("NAAT (1 day)", "NAAT (2 day)", "NAAT (3 day)", 
                                    "Ag-RDT (95% sens)", "Ag-RDT (85% sens)"),
       col=c(mycolors[1],mycolors[1],mycolors[1],mycolors[2],mycolors[2]), lty=c(1,3,5,4,1),
       lwd = c(2,1,1,1,2), seg.len = 3)



# Figure 4

layout(t(array(c(1,1,1,2,2,2,2), dim=c(1,7))))
par(mar=c(3,4,3,2), oma=c(2,0,0,0))

plot(pt, probthreshold3[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability for intervention", ylab="Net benefit", 
     main="Outpatient setting",
     ylim=c(0, max(probthreshold3, na.rm=T)*0.9),
     xlim=c(0,0.6),
     xaxt='n', lty=5, lwd=2)
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0.05,0.08, "A", font=2, cex=1.5)
lines(pt, probthreshold3[2,,1], col=mycolors[2], lwd=2)
lines(pt, probthreshold3[3,,1], col=mycolors[3], lwd=1.5, lty=3)
lines(pt, probthreshold3[5,,1], col='black', lty=3, lwd=1.5)
abline(h = 0, col='black', lty=4)
lines(pt, probthreshold4[3,,1], col=mycolors[3], lty=1, lwd=1.5)
lines(pt, probthreshold4[5,,1], col='black', lty=1, lwd=1.5)
# legend(x = "bottomright", legend = c(
#   # "NAAT (3 day)","Ag-RDT (85% sens)",
#   "Clinical, full", "Clinical, reduced",
#   "Treat all, full", "Treat all, reduced", "Treat none"),
#   # col=c(mycolors[1],mycolors[2],mycolors[3], mycolors[3],'black','black','black'), lty=c(5,1,1,3,1,3,4),
#   col=c(mycolors[3], mycolors[3],'black','black','black'), lty=c(3,1,3,1,4),
#   lwd = c(1.5,1.5,1.5,1.5,1), seg.len=2.5)

plot(pt, probthreshold2[1,,1], type='l', col=mycolors[1],
     xlab="Threshold probability for intervention", ylab="Net benefit", 
     main="Hospital setting",
     ylim=c(0, max(probthreshold2, na.rm=T)), xlim=c(0,0.6),
     xaxt='n', lwd=2)
axis(side=1, at=seq(0,1,by=0.2), labels = c("0%", "20%", "40%", "60%", '80%', "100%"))
text(0.05,0.195, "B", font=2, cex=1.5)
lines(pt, probthreshold2[2,,1], col=mycolors[2], lwd=2)
lines(pt, probthreshold2[3,,1], col=mycolors[3], lwd=1.5, lty=3)
lines(pt, probthreshold2[5,,1], col='black', lty=3, lwd=1.5)
abline(h = 0, col='black', lty=4)
lines(pt, probthreshold[3,,1], col=mycolors[3], lty=1, lwd=1.5)
lines(pt, probthreshold[5,,1], col='black', lty=1, lwd=1.5)
legend(x = "bottomright", legend = c(
  "NAAT","Ag-RDT",
  "Clinical, full", "Clinical, reduced",
                                  "Treat all, full", "Treat all, reduced", "Treat none"),
       col=c(mycolors[1], mycolors[2], mycolors[3], mycolors[3],'black','black','black'), lty=c(5,1,3,1,3,1,4),
  lwd=c(2,2,1.5,1.5,1.5,1.5,1), seg.len = 2)

mtext(side=1, outer=T, text = "Threshold probability for intervention", line=0.5, cex = 0.8)


pt
probthreshold2[5,,1]
probthreshold2[1,,1]
pt
probthreshold3[5,,1]
probthreshold3[2,,1]


##### One-way sensitivity analysis, NAAT vs RDT. (Hosp setting) #####

gethilo <- function(params, level, hosp=T)
{
  hiloparams <- list()
  hiloparams$clinicalVsTransmissionBenefit <- (if(hosp) c(5,0.2) else c(0.1,0))
  hiloparams$prev <- (if(hosp) c(0.6,0.1) else c(0.2,0.02))
  hiloparams$falseDiagnosisHarm <- c(0.2,0.01)
  hiloparams$clinicalDiagnosisDiscount <- c(1,0.5)
  hiloparams$presx <- c(0.5,0.12)
  hiloparams$Period <- c(8,3)
  hiloparams$isolationOfKnownCase <- (if (hosp) c(1,0.5) else c(0.9,0.4))
  hiloparams$isolationOfCasePendingResult <- (if(hosp) c(0.9, 0.4) else c(0.8, 0))
  hiloparams$isolationOfContacts <- c(0.7,0)
  hiloparams$sxOnsetDay <- c(4,1)
  hiloparams$medianPresentationDay <- (if(hosp) c(8, 3) else c(5,2))
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




###### sensitivity, net benefit  ##########
set.seed(55555)
hospplot <- F
resetparams <- make.params(hosp=hospplot)
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
dathigh <- melt(NBhigh["assayratio",]); dathigh$val <- "High parameter value"
dathighclin <- melt(NBhigh["clinratio",]); dathighclin$val <- "High parameter value"
datlow <- melt(NBlow["assayratio",]); datlow$val <- "Low parameter value"
datlowclin <- melt(NBlow["clinratio",]); datlowclin$val <- "Low parameter value"
datlow$names <- datlowclin$names <- rownames(datlow); dathigh$names <- dathighclin$names <- rownames(dathigh)
baseline <- NBref["RDT"]/NBref["NAAT"]; baselineclin <- NBref["RDT"]/NBref["clin"]

require(dplyr)
dat <- rbind(datlow, dathigh) 
datclin <- rbind(datlowclin, dathighclin) 
keepnames <- setdiff(dathigh$names[pmax(dathigh$value/baseline, datlow$value/baseline, dathigh$value/datlow$value, datlow$value/dathigh$value) >
                                     1.05], excludeparams)
keepnamesclin <- setdiff(dathighclin$names[pmax(dathighclin$value/baselineclin, datlowclin$value/baselineclin, dathighclin$value/datlowclin$value, datlowclin$value/dathighclin$value) > 
                                             1.1], excludeparams_clin)
dat <- subset(dat, (dat$names %in% keepnames) ); datclin <- subset(datclin, (datclin$names %in% keepnamesclin) )
dat <- dat %>% arrange(names); datclin <- datclin %>% arrange(names)

# and include a version with extremes of multiple parameters
Nparams <- Rparams <- Nparamsclin <- Rparamsclin <- Nparamsclin2 <- Rparamsclin2  <- resetparams
for (name in keepnames)
   if (dat %>% filter(names==name, val=="High parameter value") %>% select(value) > 
                            dat %>% filter(names==name, val=="Low parameter value") %>% select(value)) 
    {Rparams[[name]] <- highparams[[name]]; Nparams[[name]] <- lowparams[[name]]} else 
    {Rparams[[name]] <- lowparams[[name]]; Nparams[[name]] <- highparams[[name]]}
for (name in keepnamesclin)
  if (datclin %>% filter(names==name, val=="High parameter value") %>% select(value) > 
      datclin %>% filter(names==name, val=="Low parameter value") %>% select(value)) 
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

# rename params for display
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

# prep for tornado plotting
class(dat) <- c("tornado", class(dat)); class(datclin) <- c("tornado", class(datclin))
attr(dat, "output_name") <- "value"; attr(datclin, "output_name") <- "value"
dat$names <- paste0(niceParamNames[dat$names], "\n(",unlist(lowparams)[dat$names], "-",unlist(highparams)[dat$names],")")
datclin$names <- paste0(niceParamNames[datclin$names], "\n(",unlist(lowparams)[datclin$names], "-",unlist(highparams)[datclin$names],")")

library(ggplot2)
library(patchwork)
source("tornadoplot.R")
if(hospplot) (rdtnaat_h <- ggplot_tornado(dat, baseline, title="Net Benefit of Ag-RDT\nversus NAAT, hospital setting") +
                ylim(c(0,1.6))+
                       ylab("Value of Ag-RDT relative to NAAT \n(values >1 favor Ag-RDT)") +
                       xlab("Parameter (range explored)") + 
                geom_hline(yintercept = 1, lty=3)
    # scale_y_log10(limits=c(0.5,5), breaks = c(0.5,1,2,4,8))
      ) else
      (rdtnaat_o <- ggplot_tornado(dat, baseline, title="Net Benefit of Ag-RDT\nversus NAAT, outpatient setting") +
         ylim(c(0,2.2))+
                ylab("") +
                xlab("Parameter (range explored)") + 
         geom_hline(yintercept = 1, lty=3)
         # scale_y_log10(limits=c(0.5,max(dat$value)), breaks = c(0.5,1,2,4,8,16,32))
       ) 
if(hospplot) (rdtclin_h <- ggplot_tornado(datclin, baselineclin, title="Net Benefit of Ag-RDT\nversus clinical judgment, hospital setting") + 
                theme(legend.position = "none")+
                ylab("Value of Ag-RDT relative to clinical judgment \n(values >1 favor Ag-RDT)") +
                 xlab("Parameter (range explored)") +
                geom_hline(yintercept = 1, lty=3)
              ) else
      (rdtclin_o <- ggplot_tornado(datclin, baselineclin, title="Net Benefit of Ag-RDT\nversus clinical judgment, outpatient setting") + theme(legend.position = "none") +
          ylab("") +
          xlab("Parameter (range explored)") +
         geom_hline(yintercept = 1, lty=3) 
          )
  #to get outpatient setting, modify hospplot and repeat above.
 rdtnaat_o + theme(legend.position = c(0.25,0.15)) +
 rdtnaat_h +   theme(legend.position = 'none') + 
           plot_layout(ncol = 1, heights = c(5,2))
 pdf(file = "DCA Fig4.pdf", width=9, height=9)
 rdtnaat_o + theme(legend.position = c(0.25,0.15)) +
   rdtnaat_h +   theme(legend.position = 'none') + 
   plot_layout(ncol = 1, heights = c(5,2))
    dev.off()

      rdtclin_o + rdtclin_h + plot_layout(ncol = 1, heights = c(2,2))
      pdf(file = "DCA Fig S6.pdf", width=9, height=9)
        rdtclin_o + theme(legend.position = c(0.7,0.15)) +
          rdtclin_h + theme(legend.position = 'none') +
          plot_layout(ncol = 1, heights = c(4,2))
        dev.off()

# for caption:
                
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
 



