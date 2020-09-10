#### inputs #####
presx <- 0.2
si <- 6
si_sd <- 4
rt <- 1
f1 <- 0.9 # transmission averted after diagnosed
f2 <- 0.5 #contact notification efficacy for averting future transmission
sxonsetday <- 2 # from onset of infectivity, on average
duration <- 10 # from symptom onset to end of 99% of cumulative infectivity, on average

# timing
presentationDay <- 0
t_RDT <- 0.25
t_NAAT <- 2

# time sensitivity discount
c <- 0.15
s <- 0.05

csens_h <- 0.9
cspec_h <- 0.5 
csens_o <- 0.8
cspec_o <- 0.5
rsens <- 0.7
rspec <- 0.995
nsens <- 0.9
nspec <- 0.995

### utility functions ####

get_diff_vectors <- function(x, y) {
  count_x <- table(x)
  count_y <- table(y)
  same_counts <- match(names(count_y), names(count_x))
  count_x[same_counts] <- count_x[same_counts] - count_y
  as.numeric(rep(names(count_x), count_x))
}

# time sensitivity
discount <- function(base, discountrate, t)
{
  return(base * exp(-discountrate*t))
}


########## Simple DCA #######
#  Simple decision curve assuming a fixed clinician sens/spec and assuming all cases are equally important to detect

  # plot a DCA for some test with fixed sens and spec:
  x <- seq(0,1,by=0.01) # threshold probabilities
  weights <- x/(1-x) # number of false positives you'd want to treat per case, at each threshold prob
  
  NB <- function(sens, spec, prev)
  {
    TP <- sens*prev
    FP <- (1-spec)*(1-prev)
    NB <- TP - FP*weights
    return(NB)
  }
  
  
  plotDCA <- function(prev=0.4, csens=0.9, cspec=0.5, rsens=0.7, rspec=0.99, nsens=0.9, nspec=0.99, plotclin=T, ymax, plotlegend=T, title="") 
  {
    if(missing(ymax)) ymax <- prev
    plot(x,NB(nsens, nspec, prev), type='l', col='green', xlim=c(0,1), ylim=c(0,ymax),
         xlab="Threshold probability above which one opts to intervene",
         ylab="Net benefit",
         main=title)
    lines(x, NB(rsens, rspec, prev), col='blue')
    if(plotclin)  lines(x, NB(csens, cspec, prev), col='red')
    lines(x, NB(sens=1, spec=0, prev), col='black', lty=1)
    lines(x, NB(sens=0, spec=1, prev), col='black', lty=2)
    if(plotlegend)legend(x = "topright", legend = c("NAAT","VAT","Clinician","Treat all", "Treat none"),
                         col=c('green','blue','red','black','black'), lty=c(1,1,1,1,2)
    )
  }
  
  par(mfrow=c(1,2))
  plotDCA(prev=0.4, csens=csens_h, cspec = cspec_h, rsens=rsens, rspec=rspec, nsens=nsens, nspec=nspec, ymax=0.5, plotlegend = F, title="Hospitalized patients")
  plotDCA(prev=0.1, csens=csens_o, cspec = cspec_o, rsens=rsens, rspec=rspec, nsens=nsens, nspec=nspec, ymax=0.5, title="Mildly symptomatic patients")
  
############ infectiousness and transmission #####
npts<-1e5 #arbitrary simulation size for transmission averted


# time distribution of transmission events
# can model multiple potential transmissions per case, and then assume only a subset actually occur
library(rriskDistributions)
weibpar <- get.weibull.par(p=c(0, presx, 0.99), q=c(0,sxonsetday,sxonsetday+duration))
# # assuming equal infectivity of all cases and assuming no overdispersion. (1-p)/p = rt --> p = 1/(rt+1)
# transmissions_case <- rgeom(n = npts, prob = 1/(rt+1))
# contactindex <- rep(1:npts, times=transmissions_case)

# timing of transmissions from index case, assuming rt=1 (can scale results later)
casetransmits <- rweibull(n = npts, shape = weibpar[1], scale = weibpar[2]) - sxonsetday

  # serial interval for secondary case's symptom onset 
  serial <- rgamma(shape = si^2/si_sd, rate = si/si_sd, n = npts)
  contactonset <- casetransmits + serial
  # using a random resampling of case transmission timing to estimate timing of transmission from contacts, 
  
  contacttransmits <- contactonset + sample(casetransmits, npts, T)
  
  bins <- 200
  x <- density(casetransmits, n=33*bins, from=-3, to=30)$x
  
  #timing of transmission events
  par(mfrow=c(1,1))
  plot(density(casetransmits, n=33*bins, from=-3, to=30), xlim=c(-3,30), col="blue", xlab="Day relative to symptom onset in case", 
       main="t=Timing of transmission events")
  lines(density(contacttransmits,n=33*bins, from=-3, to=30), col="red")
  lines(x, 1/(1+rt)*density(casetransmits, n=33*bins, from=-3, to=30)$y + rt/(1+rt)*density(contacttransmits, n=33*bins, from=-3, to=30)$y, col="purple")
  legend(x="topright", legend = c("Transmission from case", "Transmission from contacts", "Combined"), lty =1, col=c("blue","red", "purple"), bty = "n")
  
  
  # timing of avertible transmission
  ## need a scaled density plot, integrating to less than 1, reflecting avertible proportion of all transmission
  ## Procedure: Identify potentially avertible at a given time, select subset (proportion f) that are averted, 
  ## and scale resulting density plot by the number averted / total number of events.
  d0_avertible_cases <- sample(which(casetransmits>0), size = round(sum(casetransmits>0)), replace=F)
  d0_averted_cases <- sample(which(casetransmits>0), size = round(f1*sum(casetransmits>0)), replace=F)
  # averted contact transmissions includes those from contacts who were never infected, and those from contact notification
  d0_averted_contacts <- union(d0_averted_cases, 
                               sample(which(contacttransmits>0), size = round(f2*sum(contacttransmits>0)), replace=F))
  d2_averted_cases <- sample(which(casetransmits>2), size = round(f1*sum(casetransmits>2)), replace=F)
  d2_averted_contacts <- union(d2_averted_cases, 
                               sample(which(contacttransmits>2), size = round(f2*sum(contacttransmits>2)), replace=F))
  
  plot(density(casetransmits, n=33*bins, from=-3, to=30), xlim=c(-3,30), col="blue", xlab="Day relative to symptom onset in case", 
       main="Avertible transmission")
  lines(density(contacttransmits,n=33*bins, from=-3, to=30), col="red")
  
  lines(x, 1/(1+rt)*length( d0_averted_cases)/npts * density(casetransmits[d0_averted_cases], n=33*bins, from=-3, to=30)$y,  col='blue', lty=2)
  lines(x, 1/(1+rt)*length( d2_averted_cases)/npts * density(casetransmits[d2_averted_cases], n=33*bins, from=-3, to=30)$y,  col='blue', lty=3)
  
  # lines(x, 1/(1+rt)*density(casetransmits, n=33*bins, from=-3, to=30)$y + rt/(1+rt)*density(contacttransmits, n=33*bins, from=-3, to=30)$y, col="purple")
  
  lines(x, 1/(1+rt)*length( d0_averted_contacts)/npts * density(contacttransmits[d0_averted_contacts], n=33*bins, from=-3, to=30)$y, col='red', lty=2)

  lines(x, rt/(1+rt)*length( d2_averted_contacts)/npts * density(contacttransmits[d2_averted_contacts], n=33*bins, from=-3, to=30)$y, 
        col='red', lty=3)
  # abline(v=0, col='gray')  
  # abline(v=2, col='gray')  
  legend(x="topright", legend = c("All", "Avertable at day 0", "Avertable at day 2",
                                  "From cases", "From their contacts"), 
         lty =c(1:3,1,1), col=c(rep("black",3),"blue","red"), bty = "n")
 
  
  ####### detection and heterogeneous infectivity ######
  
  #( Will model variation in infectiousness at diagnosis, and assume that corresponds to variation in avertable transmission -- either because the low ones were always low so weren't very infectious and don't ahve many infected contacts, or because you're cathcing them late when they aren't infectious and their contacts' infectious time has mostly passed.)
  # Dca infxn conservative assumption infectiousness ~ prob(CX+). More likely grows w VL, e.g.~Ct. 
  # Fit log-gnormal distribution of viral qty, between 10^2 and 10^8 copies. 
  # (https://www.bmj.com/content/bmj/369/bmj.m1443.full.pdf)
  normpar <- get.norm.par(p = c(0.025, 0.975), q = c(3,8), plot = F)
  logvirus <- rnorm(n = npts, mean=normpar[1], sd = normpar[2])

  sens_RDT <- 0.7
  sens_NAAT <- 0.9
  LODs <-  quantile(logvirus, 1-c(sens_NAAT, sens_RDT)) 
  
  probPos_NAAT <- 1 - (1/10)^(pmax(0, logvirus - LODs[1]))# ** add to methods
  probPos_RDT <- 1 - (1/10)^(pmax(0, logvirus - LODs[2]))
  
  col1 <-  rgb(255,200,200,max = 255, alpha = 80)
  col2 <- rgb(173,216,230,max = 255, alpha = 80)
  
  
  hist(probPos_NAAT, col=col1)
  hist(probPos_RDT, add=T, col=col2)
  
  # and model infectiousness as proportional to probability of positive culture: 
  # from ~0.1 at 10^4 to ~0.9 at 10^9 (https://wwwnc.cdc.gov/eid/article/26/10/20-2403-f1)
  # similarly from 0 to 100 over ~20 ct = log(2^20)/log(10) = 6 log10s (https://link.springer.com/article/10.1007%2Fs10096-020-03913-9)
  infectivityScale <- 2 # 1 = log-linear -- how skewed is the infectivity distribution, in relation to viral load?
  infectivity <- (pmax(0, pmin(logvirus-3, 9-3)/(9-3)))^infectivityScale# ** add to methods
  infectivity <- infectivity/mean(infectivity)

  hist(infectivity, col='black', breaks = seq(0,ceiling(max(infectivity)),by=0.1))
  hist(infectivity[rbinom(npts,1,prob=probPos_NAAT)==1], 
       add=T, breaks = seq(0,ceiling(max(infectivity)),by=0.1))
  hist(infectivity[rbinom(npts,1,prob=probPos_RDT)==1], 
  # hist(infectivity[rbinom(prob=probPos_RDT)], 
       add=T, breaks = seq(0,ceiling(max(infectivity)),by=0.1), col='white')
    

  #### transmission events, with weighting and LOD #######
    
  # To get the transmission directly from cases, weight each case by relative infectivity 
  # Get indexes of cases who could generate transmission
  caseTransmissions <- sample(1:length(casetransmits), size = 1e5, replace=T, prob = infectivity)
  
  # to get the probability that each transmission event is averted by a given test, 
  # screen by TAT and a binomial probability that the assay detected them.
  caseTransmissionAvertable_RDT <- caseTransmissions[which(casetransmits[caseTransmissions]>presentationDay + t_RDT)] 
  caseTransmissionAverted_RDT <- caseTransmissionAvertable_RDT[rbinom(n=length(caseTransmissionAvertable_RDT), 
                                              size=1, prob = f1*probPos_RDT[caseTransmissionAvertable_RDT]) == 1]
  length(caseTransmissionAverted_RDT)/length(caseTransmissionAvertable_RDT) # should be >sens*f1

  caseTransmissionAvertable_NAAT <- caseTransmissions[which(casetransmits[caseTransmissions]>presentationDay + t_NAAT)] 
  caseTransmissionAverted_NAAT <- caseTransmissionAvertable_NAAT[rbinom(n=length(caseTransmissionAvertable_NAAT), 
                                                              size=1, prob = f1*probPos_NAAT[caseTransmissionAvertable_NAAT]) == 1]
  length(caseTransmissionAverted_NAAT)/length(caseTransmissionAvertable_NAAT) # should be >sens*f1
  
  par(mfrow=c(1,2))
  hist(casetransmits[caseTransmissions], col="white")
  # hist(casetransmits[caseTransmissionAverted_NAAT], add=T, col='lightblue', alpha=0.5)
  hist(casetransmits[get_diff_vectors(caseTransmissions, caseTransmissionAverted_NAAT)], add=T, col=col1)
  hist(casetransmits[get_diff_vectors(caseTransmissions, caseTransmissionAverted_RDT)], add=T, col=col2)
  # hist(casetransmits[caseTransmissionAverted_RDT], add=T, col=NULL)
  
  # for contact transmissions averted, first determine which contacts still get infected 
   # (based on case transmission averted, weighted y case infecrivity), 
  # and then determine whether transmission from the infected contacts is averted 
   # based on the sensitivity and TAT of the case's test
  ## Contact trnasmission events timing distribution:
  hist(contacttransmits[caseTransmissions], col="white")
  
  # the transmissions that aren't avertable have cases that aren't averted and (have transmission after TAT or are missed by f)
  caseTransmits_RDT <- get_diff_vectors(caseTransmissions, caseTransmissionAverted_RDT)
  contactTransmissionAvertable_RDT <- caseTransmits_RDT[contacttransmits[caseTransmits_RDT]>presentationDay + t_RDT]
  
  contactTransmissionAverted_RDT <- contactTransmissionAvertable_RDT[rbinom(n=length(contactTransmissionAvertable_RDT), 
                                       size=1, prob = f2*probPos_RDT[contactTransmissionAvertable_RDT]) == 1]

  caseTransmits_NAAT <- get_diff_vectors(caseTransmissions, caseTransmissionAverted_NAAT)
  contactTransmissionAvertable_NAAT <- caseTransmits_NAAT[contacttransmits[caseTransmits_NAAT]>presentationDay + t_NAAT]
  
  contactTransmissionAverted_NAAT <- contactTransmissionAvertable_NAAT[rbinom(n=length(contactTransmissionAvertable_NAAT), 
                                                                            size=1, prob = f2*probPos_NAAT[contactTransmissionAvertable_NAAT]) == 1]
  hist(contacttransmits[get_diff_vectors(caseTransmits_NAAT, contactTransmissionAverted_NAAT)], add=T, col=col1)
  hist(contacttransmits[get_diff_vectors(caseTransmits_RDT, contactTransmissionAverted_RDT)], add=T, col=col2)
  
  
  
  # ** make this a function
  

##### Calculate Net Benefit ######
    
    
  Then make separate plots for each, and the whole thing becomes 3D
  
  
  For both, downweight by how much of that will have passed by the time you get a result.
  
  
  
  Finally, discuss the net benefit as something you can weigh relative to cost. Or can apply to RDT vs clinical after you've maxed your PCR capacity. 




###### Model of clinician judgment #########

Assume clinicians evaluate patients based on clinicial features (incl symptoms, exposure risks, and potential for alternative diagnoses) and could hypothetically rank patients on some arbitary "liklihood of COVID" scale. This ranking may depend on individual-level exposure risk but will be independent of the general prevalence of COVID in the population. In keeping with common approaches to ROC analysis, assume that patients with COVID have clinical-judgment liklihoods normally distributed as $N(m_p,\sigma_p)$, and patients without COVID follow another normal distribution $N(m_n, \sigma_n)$. 

We'll start by assuming equal variances for the positive and negative populations. Take one observed clinician sens and spec (sens 0.95, spec 0.5, in a population with prev 0.1). Set arbitrary means, and fit the variance so that such a point exists. 
  
  ```{r}
  # These are on an arbitrary scale, so set an arbitrary mean and sd for the first
  m_p <- 0
  sd_p <- 1
  sensquantile <- qnorm(p = sens, mean = 0, sd=1, lower.tail = F)
  
  # so we want sensquantile of the first distribution, to equal specquantile of the second (lowertail)
  sd_n <- sd_p # assuming equal for now, but can adjust.
  # Now find the mean for the negatives:
  # m_n + sd_n * qnorm(p = 0.5, mean=0, sd=1) = sensquantile
  m_n = sensquantile - sd_n * qnorm(p = spec, mean=0, sd=1)
  
  # Plot the distributions, at two levels of prevalence:
  plothist <- function(prev, nsim, title)
  { L_p <- rnorm(n = nsim*prev, mean = m_p, sd = sd_p)
  L_n <- rnorm(n = nsim*(1-prev), mean = m_n, sd = sd_n)
  hist(L_n, freq = T, breaks=seq(-6,6,by=0.1), main=title)
  hist(L_p, freq = T, add=T, col='red', breaks=seq(-6,6,by=0.1) )
  }
  par(mfrow=c(1,2))
  plothist(prev = 0.1, nsim = 10000, title="10% prevalence")
  plothist(prev = 0.5, nsim = 10000, title="50% prevalence")
  
  prev <- 0.1
  nsim <- 10000
  L_p <- rnorm(n = nsim*prev, mean = m_p, sd = sd_p)
  L_n <- rnorm(n = nsim*(1-prev), mean = m_n, sd = sd_n)
  plot(ecdf(x = L_p))
  lines(ecdf(x = L_n), col='red')
  abline(v=0, lty=2)
  abline(v=-1, lty=3, col='gray')
  abline(v=0.5, lty=3, col='gray')
  
  ```
  
  
  So how do clinicians adjust their judgment depending on threshold prob and prev?
    
    If a clinician knows what the above plot looks like, they can choose the point which has the PPV that they want. 
  
  ```{r}
  nsim <- 10000
  prev <- 0.5
  scores <-data.frame(c(rnorm(n = nsim*prev, mean = m_p, sd = sd_p),
                        rnorm(n = nsim*(1-prev), mean = m_n, sd = sd_n)))
  colnames(scores) <- "score"
  scores$covid <- rep(c(1,0), times=nsim*c(prev, 1-prev))
  
  x
  weights
  
  # for each x, find the cutoff that maximizes NB. 
  
  NBdiscrete <- function(cutoff, thresholdprob)
  {
    weight <- thresholdprob/(1-thresholdprob)
    TP <- sum(scores$score>cutoff & scores$covid==1)
    FP <- sum(scores$score>cutoff & scores$covid==0)
    NB <- TP - weight*FP
    
    return(NB)
  }
  
  #example for a specific threshold:
  tprob <- 0.2
  i <- which(x==tprob)
  cutoffNBs <- numeric(length(scores$score))
  for (j in 1:length(cutoffNBs))
  {
    cutoffNBs[j] <- NBdiscrete(scores$score[j], x[i])
  }
  scores$prevhalfcutoffs <- cutoffNBs
  library(ggplot2)
  ggplot(data=scores, aes(x=score, fill=factor(covid))) + 
    geom_histogram(alpha=0.5, position="identity", binwidth = 0.1) + 
    geom_line(aes(x=score, y=prevhalfcutoffs))
  
  
  choosescore <- numeric(length(x))
  considercutoffs <- seq(floor(min(scores$score)), ceiling(max(scores$score)), by = 0.1)
  for (i in 1:length(x))
  {
    cutoffNBs <- numeric(length(considercutoffs))
    for (j in 1:length(cutoffNBs))
    {
      cutoffNBs[j] <- NBdiscrete(considercutoffs[j], x[i])
    }
    choosescore[i] <- considercutoffs[which.max(cutoffNBs)] # account for uncertainty here somewhere?
  }
  
  # what score cutoff should be used for each threshold probability, to maximize net benefit?
  plot(x, choosescore)
  NBs <- numeric(length(x))
  for (i in 1:length(NBs)) NBs[i] <- NBdiscrete(cutoff = choosescore[i], thresholdprob = x[i])
  plot(x, NBs)
  
  
  plotDCA(prev=prev)
  plotDCA(prev=prev, plotclin = F)
  lines(x, NBs/nsim, col='red')
  
  
  ```
  
  [How to account for error/imprecision in the clinician's calibration?]

[Or for RDT, account for there being some the clinician would judge so probable that they wouldn't do the RDT??]
  
  [For next piece: Perspective with multiple issues in decision curves (wo too much depth on any one), OR methods piece just on the baseline of DCA that is better than treat all/treat none and that some get treated regardless]
  
  
  
  