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
    
    #   # by timing of transmission, cases
    #   arranged <- caseevents %>% 
    #     arrange(value) %>% 
    #     mutate(rn = row_number())
    #   arranged_RDT <- subset(caseevents, Transmits_RDT==1) %>% 
    #     arrange(value) %>% 
    #     mutate(rn = row_number())
    #   arranged_NAAT <- subset(caseevents, Transmits_NAAT==1) %>% 
    #     arrange(value) %>% 
    #     mutate(rn = row_number())
    #   arranged_clin <- subset(caseevents, Transmits_clin==1) %>% 
    #     arrange(value) %>% 
    #     mutate(rn = row_number())
    #   
    #   # by timing of transmission, contacts
    #   gen2 <- caseevents %>% 
    #     arrange(contacttime) %>% 
    #     mutate(rn = row_number())
    #   gen2_RDT <- subset(caseevents, contactTransmits_RDT==1) %>% 
    #     arrange(contacttime) %>% 
    #     mutate(rn = row_number())
    #   gen2_NAAT <- subset(caseevents, contactTransmits_NAAT==1) %>% 
    #     arrange(contacttime) %>% 
    #     mutate(rn = row_number())
    #   gen2_clin <- subset(caseevents, contactTransmits_clin==1) %>% 
    #     arrange(contacttime) %>% 
    #     mutate(rn = row_number())
    #   
    #   # by timing of presentation, for all cases not weighted as source cases**
    #   allcases <- data.frame(presentationDay)
    #   allcases$RDTDay <- presentationDay + turnaroundTimeRDT
    #   allcases$NAATDay <- presentationDay + turnaroundTimeNAAT
    #   allcases$Pos_RDT <- Pos_RDT
    #   allcases$Pos_NAAT <- Pos_NAAT
    #   allcases$Pos_clin <- Pos_clin
    #   tested <- allcases %>% 
    #     arrange(presentationDay) %>% 
    #     mutate(rn = row_number())
    #   detected_RDT <- subset(allcases, Pos_RDT==1) %>% 
    #     arrange(RDTDay) %>% 
    #     mutate(rn = row_number())
    #   detected_NAAT <- subset(allcases, Pos_NAAT==1) %>% 
    #     arrange(NAATDay) %>% 
    #     mutate(rn = row_number())
    #   detected_clin <- subset(allcases, Pos_clin==1) %>% 
    #     arrange(presentationDay) %>% 
    #     mutate(rn = row_number())
    #   
    #   
    #   
    #   # plot transmission
    #   colors <- c("NAAT" = mycolors[1], "Ag-RDT" = mycolors[2], "Clinical diagnosis" = mycolors[3], 
    #               "No intervention" = mycolors[4])
    #   linetypes <- c("From index cases" = "dashed", "From contacts" = "dotdash")
    #   (transmission <- ggplot() +
    #       xlim(-3,20) + 
    #       geom_step(data=arranged[seq(1,nrow(arranged), by=100),], aes(x=value, y=rn, color="No intervention", linetype="From index cases"),size=0.8) +
    #       geom_step(data=arranged_clin[seq(1,nrow(arranged_clin), by=100),], aes(x=value, y=rn, color='Clinical diagnosis', linetype="From index cases"),size=0.8) +
    #       geom_step(data=gen2_clin[seq(1,nrow(gen2_clin), by=100),], aes(x=contacttime, y=rn, color='Clinical diagnosis', linetype="From contacts"),size=0.8) +
    #       geom_step(data=arranged_RDT[seq(1,nrow(arranged_RDT), by=100),], aes(x=value, y=rn, color="Ag-RDT", linetype="From index cases"),size=0.8) +
    #       geom_step(data=arranged_NAAT[seq(1,nrow(arranged_NAAT), by=100),], aes(x=value, y=rn, color='NAAT', linetype="From index cases"),size=0.8)  +
    #       geom_step(data=gen2[seq(1,nrow(gen2), by=100),], aes(x=contacttime, y=rn, color="No intervention", linetype="From contacts"),size=0.8) +
    #       geom_step(data=gen2_RDT[seq(1,nrow(gen2_RDT), by=100),], aes(x=contacttime, y=rn, color="Ag-RDT", linetype="From contacts"),size=0.8) +
    #       geom_step(data=gen2_NAAT[seq(1,nrow(gen2_NAAT), by=100),], aes(x=contacttime, y=rn, color='NAAT', linetype="From contacts"),size=0.8)  +
    #       labs(x = "Days from symptom onset",
    #            color = "Diagnostic approach",
    #            linetype="Type of transmission event") +
    #       scale_color_manual(values = c(mycolors[4],mycolors[2],mycolors[1],mycolors[3]), breaks=c("No intervention","Ag-RDT", "NAAT", "Clinical diagnosis") )+
    #       scale_linetype_manual(values = c("dashed","dotdash"), breaks=c("From index cases", "From contacts") ) +
    #       scale_y_continuous(name = "Cumulative transmission events", 
    #                          breaks = npts*c(0,0.25,0.5,0.75,1), 
    #                          labels = c("0%","25%","50%","75%","100%")) + 
    #       theme_bw() +
    #       theme(legend.key.width=unit(1,"cm"))
    #   )

      
  library(reshape2)
  times <- melt(data = caseevents, id.vars = c("name", "Transmits_RDT", "contactTransmits_RDT"), measure.vars = c("value", "contacttime"))
  times_rdt <- melt(data = subset(caseevents, Transmits_RDT==1 | contactTransmits_RDT==1), id.vars = c("name", "Transmits_RDT", "contactTransmits_RDT"), measure.vars = c("value", "contacttime"))
  times_naat <- melt(data = subset(caseevents, Transmits_NAAT==1 | contactTransmits_NAAT==1), id.vars = c("name", "Transmits_NAAT", "contactTransmits_NAAT"), measure.vars = c("value", "contacttime"))
  
  times_rdt <- subset(times_rdt, (variable=="value" & Transmits_RDT==1)|(variable=="contacttime" & contactTransmits_RDT))
  times_naat <- subset(times_naat, (variable=="value" & Transmits_NAAT==1)|(variable=="contacttime" & contactTransmits_NAAT))
  

  # # plot daily (not comulative) transmission
  # (dailytrans_case <- ggplot(times, aes(x=value, fill=variable, col=variable)) +
  #     geom_freqpoly(binwidth = 1, aes(y=..count../100)) +
  #     geom_freqpoly(data=times_rdt, binwidth = 1, aes(y=..count../100), lty=2) +
  #     geom_freqpoly(data=times_naat, binwidth = 1, aes(y=..count../100), lty=3) +
  #     ylab("Transmission eventss per 1000 patients per day") +
  #     xlab("Days since symptom onset")
  # )


    # return(dailytrans_case)
    return(list("noint" = times, "rdt" = times_rdt, "naat" = times_naat))
  })
}



colors <- c("No intervention" = "gray", 
            "Immediate presentation, fully test-driven management" = mycolors[1], 
            "Delayed presentation, fully test-driven management" = mycolors[2],
            "Delayed presentation, some empiric isolation (final model)" = mycolors[4])

linetypes <- c("None" = 1,
               "Ag-RDT" = 2,
               "NAAT" = 3)

params <- make.params(hosp=F)
base <- plotTransmissions(params = params)
(dailytrans_case <- ggplot(base$noint, aes(x=value, fill=variable, lty="None")) +
    geom_freqpoly(binwidth = 1, aes(y=..count../100, 
                  col='No intervention', lty="None")) +
    geom_freqpoly(data=base$rdt, binwidth = 1, aes(y=..count../100, 
                  col='Delayed presentation, some empiric isolation (final model)', lty="Ag-RDT")) +
    geom_freqpoly(data=base$naat, binwidth = 1, aes(y=..count../100, 
                  col='Delayed presentation, some empiric isolation (final model)', lty="NAAT")) +
    labs(x = "Days since symptom onset",
         y = "Transmission events per 1000 patients per day") + 
    theme_bw() +
    scale_color_manual(name = "Intervention assumptions", values = colors) + 
    scale_linetype_manual(name= "Diagnostic test", values = linetypes)
)

# rdt vs naat with: 

## if presents immediately and treated fully basd on test, but tests have diagnostic delay (NAAT>RDT)
sxparams <- within(params, {isolationOfContacts=1; isolationOfKnownCase=1; 
  isolationOfCasePendingResult=0
    medianPresentationDay=0})
fulltestsx <- plotTransmissions(params = sxparams)
(dailytrans_case2 <- dailytrans_case + 
    geom_freqpoly(data=fulltestsx$rdt, binwidth = 1, aes(y=..count../100, 
                  col='Immediate presentation, fully test-driven management', lty='Ag-RDT')) +
    geom_freqpoly(data=fulltestsx$naat, binwidth = 1, aes(y=..count../100, 
                  col='Immediate presentation, fully test-driven management', lty='NAAT'))
)


# delay in presentation, but treat fully based on test
tparams <- within(params, {isolationOfCasePendingResult=0; isolationOfContacts=1; isolationOfKnownCase=1})
fulltest <- plotTransmissions(params = tparams)
(dailytrans_case3 <- dailytrans_case2 + 
    geom_freqpoly(data=fulltest$rdt, binwidth = 1, aes(y=..count../100, 
                  col='Delayed presentation, fully test-driven management', lty='Ag-RDT')) +
    geom_freqpoly(data=fulltest$naat, binwidth = 1, aes(y=..count../100, 
                  col='Delayed presentation, fully test-driven management', lty='NAAT'))
)

Transmission_source <- c(
    value = "From index cases",
    contacttime = "From immediate contacts"
  )

dailytrans_case3 + facet_wrap("variable", nrow=2, labeller = labeller(
  variable = Transmission_source
))


# Then show cumulative sensitivty and transmission as bar graph
params <- make.params(hosp=F)
sensparam <- c(params$sensitivityNAAT, 
               params$sensitivityNAAT*params$sensitivityRDT_vsNAAT,
               params$sensitivityClinician)
senstable <- data.frame(sensparam, row.names = c("NAAT", "RDT", "clin"))
senstable$detected <- 
  (mapTransmissions(params, plots=F))[c("NAATdetection", "RDTdetection", "clindetection")]
senstable$infweighted <- 
  (mapTransmissions(params, plots=F))[c("NAATpeakinf", "RDTpeakinf", "clinpeakinf")]
senstable$futureinfweighted <- 
  (mapTransmissions(params, plots=F))[c("NAATavertable", "RDTavertable", "clinavertable")]
senstable$prevented <- 
  (mapTransmissions(params, plots=F))[c("NAATaverted", "RDTaverted", "clinaverted")]
senstable$test <- rownames(senstable)
longsenstable <- melt(senstable, id.vars='test')

ggplot(senstable) + geom_bar()
