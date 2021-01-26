#Modified version of maptransmissions function, for plotting transmission under different assumptions (for Fig 3)

plotTransmissions <- function(params,
                             npts = 1e5 )#arbitrary sim size) 
{  
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
    
    if(medianPresentationDay==0) presentationDay = rep(0, npts) else
        presentationDay <- rtweibull(n = npts, shape=0.85, scale=5.5*medianPresentationDay/3.5, endpoint = 25)
    presentationStep <- unlist(lapply(presentationDay, function(x) which.min(timing<x)))
    
    Pos_NAAT_onset <- 1*(logvirus>LODs_delay_NAAT[1])
    Pos_RDT_onset <- 1*(logvirus>LODs_delay_RDT[1])
    
    Pos_NAAT <- 1*(logvirus>LODs_delay_NAAT[presentationStep])
    Pos_RDT <- 1*(logvirus>LODs_delay_RDT[presentationStep])
    Pos_clin <- rbinom(n = npts, size = 1, prob = sensitivityClinician)
    
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
    
    ## Contact transmission events timing distribution:
    
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
                                                                                                        size=1, prob = isolationOfContacts*Pos_clin*clinicalDiagnosisDiscount)
    caseevents$contactTransmits_clin <- caseevents$Transmits_clin & !caseevents$contactTransmissionAverted_clin
    
    # casetransmits are timing, caseTransmissions are indexes fo all events (some repeated), 
    # caseTransmits_X are indexed for events that occur with a given test (subset of caseTransmissions)
    
    caseevents$Transmission <- ifelse(caseevents$Transmits_RDT==0, ifelse(caseevents$Transmits_NAAT==0, "Prevented with either test", "Prevented with Ag-RDT only"), 
                                      ifelse(caseevents$Transmits_NAAT==0, "Prevented with NAAT only", "Occurs regardless of test"))
    caseevents %>% dplyr::count(Transmits_NAAT, Transmits_RDT, Transmission)
    caseevents$Transmission <- factor(caseevents$Transmission, levels=c("Occurs regardless of test",
                                                                        "Prevented with NAAT only",
                                                                        "Prevented with either test",
                                                                        "Prevented with Ag-RDT only"))
    
   
  library(reshape2)
  times <- melt(data = caseevents, id.vars = c("name", "Transmits_RDT", "contactTransmits_RDT"), measure.vars = c("value", "contacttime"))
  times_rdt <- melt(data = subset(caseevents, Transmits_RDT==1 | contactTransmits_RDT==1), id.vars = c("name", "Transmits_RDT", "contactTransmits_RDT"), measure.vars = c("value", "contacttime"))
  times_naat <- melt(data = subset(caseevents, Transmits_NAAT==1 | contactTransmits_NAAT==1), id.vars = c("name", "Transmits_NAAT", "contactTransmits_NAAT"), measure.vars = c("value", "contacttime"))
  times_clin <- melt(data = subset(caseevents, Transmits_clin==1 | contactTransmits_clin==1), id.vars = c("name", "Transmits_clin", "contactTransmits_clin"), measure.vars = c("value", "contacttime"))
  
  times_rdt <- subset(times_rdt, (variable=="value" & Transmits_RDT==1)|(variable=="contacttime" & contactTransmits_RDT))
  times_naat <- subset(times_naat, (variable=="value" & Transmits_NAAT==1)|(variable=="contacttime" & contactTransmits_NAAT))
  times_clin <- subset(times_clin, (variable=="value" & Transmits_clin==1)|(variable=="contacttime" & contactTransmits_clin))
  

    # return(dailytrans_case)
    return(list("noint" = times, "rdt" = times_rdt, "naat" = times_naat, "clin" = times_clin))
  })
}


colors <- c(
  "NAAT-guided" = mycolors[1], "Ag-RDT-guided" = mycolors[2], 
  "Clinical judgment-based" = mycolors[3],
  "No intervention" = 'gray')

Transmission_source <- c(
  value = "Transmission from index cases",
  contacttime = "Transmission from their direct contacts" )

dailytrans_case <- list()
for (hosp in c(T, F))
{
  params <- make.params(hosp=hosp)
  base <- plotTransmissions(params = params)
  width <- 0.7
  (dailytrans_case[[paste0("hosp",hosp)]] <- ggplot(base$noint, aes(x=value, fill=variable, col="No intervention")) +
      geom_freqpoly(binwidth = 1, center=0, aes(y=..count../100, 
                                                col='No intervention'), lty='33', lwd=width) +
      geom_freqpoly(data=base$rdt, binwidth = 1, center=0, aes(y=..count../100, 
                                                               col="Ag-RDT-guided"), lty='33', lwd=width) +
      geom_freqpoly(data=base$naat, binwidth = 1,center=0,  aes(y=..count../100, 
                                                                col="NAAT-guided"), lty='34', lwd=width) +
      geom_freqpoly(data=base$clin, binwidth = 1,center=0,  aes(y=..count../100, 
                                                                col="Clinical judgment-based"), lty='22', lwd=width) +
      labs(x = "Days since symptom onset",
           y = "Transmission events per 1000 patients per day") + 
      theme_bw() +
      scale_color_manual(name= "", values = colors, guide = guide_legend(reverse = TRUE)) + 
      facet_wrap("variable", nrow=2, labeller = labeller(
    variable = Transmission_source
      )) + 
    xlim(c(-3, 25))  + 
   theme(legend.position = c(0.7,0.8)))
}

dailytrans_case$hospFALSE
dailytrans_case$hospTRUE

#### Then show cumulative sensitivty and transmission as bar graph ###
# params <- make.params(hosp=F)
params <- make.params(hosp=T)
sensparam <- c(params$sensitivityNAAT, 
               params$sensitivityNAAT*params$sensitivityRDT_vsNAAT,
               params$sensitivityClinician)
senstable <- data.frame(sensparam, row.names = c("NAAT", "Ag-RDT", "Clinical judgment"))
senstable$detected <- 
  unlist((mapTransmissions(params, plots=F))[c("NAATdetection", "RDTdetection", "clindetection")])
senstable$futureinfweighted <- 
  unlist((mapTransmissions(params, plots=F))[c("NAATfuturetrans", "RDTfuturetrans", "clinfuturetrans")])
senstable$avertable <- 
  unlist((mapTransmissions(params, plots=F))[c("NAATavertable", "RDTavertable", "clinavertable")])

unlist((mapTransmissions(params, plots=F))[c("NAATpeakinf", "RDTpeakinf", "clinpeakinf")])
unlist((mapTransmissions(params, plots=F))[c("NAATcontactavertable", "RDTcontactavertable", "clincontactavertable")])
unlist((mapTransmissions(params, plots=F))[c("allNAATavertable", "allRDTavertable", "allclinavertable")])
unlist((mapTransmissions(params, plots=F))[c("allNAATcontactavertable", "allRDTcontactavertable", "allclincontactavertable")])

senstable$test <- rownames(senstable)

longsenstable <- melt(senstable, id.vars=c('test'))
require(plyr)
longsenstable$variable <- revalue(longsenstable$variable, c(
  "sensparam"="Sensitivity during acute illness\n(input parameter)", 
  "detected" = "Sensitivity at clinical presentation,\nunweighted (simulated)",
  "futureinfweighted" = "Sensitivity at clinical presentation,\nweighted by future transmission\n(simulated, accounting for infectivity)",
  "avertable" = "Avertable future transmission\n(simulated, accounting for\ninfectivity and diagnostic delay)"))
longsenstable$variable <- factor(longsenstable$variable, levels = rev(levels(longsenstable$variable)))
  
longsenstable %>% 
  filter(test %in% c("NAAT", "Ag-RDT")) %>%
  mutate(test = factor(test, levels = c("NAAT", "Ag-RDT"))) %>%
  ggplot() + geom_col(aes(x=variable, y=value*100, fill=test), position='dodge') + 
 ylab("% detected or avertable") + xlab("")  + 
  ylim(c(0,100)) +
  geom_text(aes( x = variable,y = value*100-5, fill=test, label = round(value*100)), 
            position = position_dodge(width=1), col='white', fontface='bold') + 
  scale_fill_manual(values = mycolors[1:2], 
                    name = element_blank(),
                    guide = guide_legend(reverse = TRUE) ) + 
  coord_flip()  + 
  theme(axis.text=element_text(size=9.5)) 
