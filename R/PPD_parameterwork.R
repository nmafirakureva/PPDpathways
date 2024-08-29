## Thinking through parameters for the PPD Model
library(HEdtree)

# TB symptom screening
# 90% (80% - 100%) based on London prisons data
(p <- getAB(0.90,((80 - 100)/392)^2))
curve(dbeta(x,p$a,p$b),from=0,to=1,n=200)
png(here("plots/phc.tbprev.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Prevalence of true TB at PHC', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Prison GP assessment
# 90% (87% - 94%) based on crude assumption
(p <- getAB(0.90,((100 - 80)/392)^2))
curve(dbeta(x,p$a,p$b),from=0,to=1,n=200)
png(here("plots/dh.tbprev.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Prevalence of true TB at DH', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Proportion of patients getting a chest x-ray
# Assuming 77% (56-88%) based on prison piolot data
(p <- getAB(0.77, ((56-88)^2)/392^2))
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Initial care-seeking at PHC', 
      ylab = 'Density', 
      main = '')
summary(rbeta(1000, p$a, p$b))

# Attending NHS referral
# Using Expected versus observed outpatient appointment rates for prisoners, 2017/18
# from Davies, Miranda, et al. "Locked out? Prisoners’ use of hospital care." Research summary, Nuffield Trust (2020).
# Expected = 100388, Observed = 76673
# Observed/Expected
# curve(dbeta(x,76673,100388-76673),n=1e3)
76673/10000
(100388-76673)/10000
fac <- 10000 #how much weaker (same mean)? I want a>1 so there is a peak within interval
curve(dbeta(x,76673/10000,(100388-76673)/10000),n=1e3)
(p <- getLNparms(2.1,(10.1-0.1)^2/3.92^2,med=FALSE))
curve(dlnorm(x,meanlog=p$mu,sdlog=p$sig),from=0,to=200, 
      xlab = 'OR for initial care-seeking at DH-level given true TB (0-4 years)', 
      ylab = 'Density', 
      main = '')
summary(rbeta(1000, 76673/10000,(100388-76673)/10000))

# Contact tracing
# Using 6 (4-10) based on xxx
(p <- getLNparms(6.5,(3.25-10)^2/3.92^2,med=FALSE))
curve(dlnorm(x,meanlog=p$mu,sdlog=p$sig),from=0,to=200, 
      xlab = 'OR for initial care-seeking at DH-level given true TB (5-14 years)', 
      ylab = 'Density', 
      main = '')
summary(rlnorm(1000,meanlog=1.83788909983218,sdlog=0.1))

# TB presumption in children with active TB disease
# Using sensitivity:89% (95% CI 52% to 98%) based on https://doi.org/10.1002/14651858.CD013693.pub2
(p <- getAB(0.89, ((98-60)^2)/392^2))  
png(here("plots/TB presumption.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'TB presumption in children with active TB disease', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# TB presumption in children without active TB disease
# Using specificity: 69% (95% CI 51% to 83%) based on https://doi.org/10.1002/14651858.CD013693.pub2
(p <- getAB(0.31, ((49-17)^2)/392^2))
png(here("plots/TB presumption no TB.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'TB presumption in children with no TB disease', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Factor to scale down presumption@ PHC
PD2 <- setDT(PD2)
summary(PD2[,dh.presumed/phc.presumed]) # Presumptive TB in IHVN data

## Bacteriological testing: 
# facilities with access to bacteriological testing * Xpert/Truenat availability * sample possibility

## Bacteriological test availability (defined as the proportion of facilities 
# with access to bacteriological testing)
# using Odume et al., 2023 to approximate this
# Plot in Odume et al., 2023 suggests 
# 0-10%? availability of Xpert Ultra in PHC
# ~20-100% availability of Xpert Ultra in DH

# DH
(p <- getAB(0.70, ((50-100)^2)/392^2))
png(here("plots/DH_testing.png.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Facilities with bacteriological testing (DH)', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# PHC, SOC without Truenat
# Could assume no test availability or use Odume et al., 2023 data
# to approximate a low coverage of Xpert Ultra @ PHC e.g. 5% (0-10%)
(p <- getAB(0.05, ((0-10)^2)/392^2))
png(here("plots/PHC_testing.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Facilities with bacteriological testing (DH)', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# PHC, SOC with Truenat introduction
# increased coverage with Truenat introduction
# Rough assumption of 60% (20-100%) based on Xpert availability in Odume et al., 2023
# But could just assume 100% coverage
(p <- getAB(0.60, ((20-100)^2)/392^2))
png(here("plots/PHC_Truenat_testing.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Truenat coverage at PHC (Intervention)', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Type of test available (in facilities with access to bacteriological testing, 
# the proportion of tests done on Xpert/Truenat)
# Currently set to depend on Xpert availability
# TODO:think of a better way to model this
# Assuming 100% Xpert Ultra in DH (no Truenat)
# Assuming 0% Xpert Ultra in PHC (SOC): maybe just assume 1 and let this be dictated by test availability & sample possibility?

# Bacteriological test possibility
# Based on possibility of obtaining a suitable (currently spontaneous sputum) sample for testing

## spontaneous sputum
# DH 5-14 years: 0.291 (0.215 - 0.375) based on TB SPEED Decent https://doi.org/10.1016/j.eclinm.2024.102528
(p <- getAB(0.291,((21.5 - 37.5)/392)^2))
png(here("plots/DH_SPO.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Sputum possibility (5-14 year) at DH',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# DH 0-4 years: 0.024 (0.020 - 0.027) based on Xpert Ultra on stool paper 10.1136/bmjopen-2021-058388
(p <- getAB(0.024,((2.0 - 2.7)/392)^2))
png(here("plots/DH_SPY.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Sputum possibility (0-4 year) at DH',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# PHC 5-14 years: 0.092 (0.063 - 0.129 based on TB SPEED Decent https://doi.org/10.1016/j.eclinm.2024.102528
(p <- getAB(0.092,((6.3 - 12.9)/392)^2))
png(here("plots/PHC_SPO.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Sputum possibility (5-14 year) at PHC',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# PHC 0-4 years: 0.024 (0.020 - 0.027) based on Xpert Ultra on stool paper 10.1136/bmjopen-2021-058388
(p <- getAB(0.024,((2.0 - 2.7)/392)^2))
png(here("plots/PHC_SPY.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Sputum possibility (0-4 year) at PHC',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

## Bacteriological test accuracy
# Sensitivity of Xpert on respiratory samples
# Using 75.3% (95% CI 64.3 to 83.8) based on Kay et al., 2022 (PMID: 36065889)
(p <- getAB(0.753, ((64.3-83.8)^2)/392^2))
png(here("plots/Xpert_sens.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Sensitivity (Xpert Ultra on respiratory samples)',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Specificity of Xpert on respiratory samples
# Using 97.1% (95% CI 94.7 to 98.5) based on Kay et al., 2020 (PMID: 36065889)
(p <- getAB(0.971, ((94.7-98.5)^2)/392^2))
curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
png(here("plots/Xpert_speci.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Specificity (Xpert Ultra on respiratory samples)',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Sensitivity of Truenat on stool
# Using 0.571 (0.48,0.659) based on Singh et al. 2023
(p <- getAB(0.571, ((48.0-65.9)^2)/392^2))
curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
png(here("plots/Truenat_sens.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Sensitivity (Truenat on sputum samples)',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Specificity of Truenat on stool
# Using 0.92 (0.892, 0.942) based on Singh et al. 2023
(p <- getAB(0.92, ((89.2-94.2)^2)/392^2))
curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
png(here("plots/Truenat_speci.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Specificity (Truenat on sputum samples)',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Sensitivity of clinical diagnosis
# Using 0.627 (0.615 - 0.639) based on Marais et al. 2006
(p <- getAB(0.627, ((61.5-65.9)^2)/392^2))
curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
png(here("plots/Clinical_sens.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Specificity (Clinical diagnosis)',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Specificity of clinical diagnosis
# Using 0.901 (0.894 - 0.908) based on Marais et al. 2006
(p <- getAB(0.901, ((89.4-90.8)^2)/392^2))
curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
png(here("plots/Clinical_speci.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Specificity (Clinical diagnosis)',
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))


# hospital referral loss to follow-up
# Assuming ~30% based on previous analyses
(p <- getAB(0.3, ((0-50)^2)/392^2))
png(here("plots/rltfu.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Hospital referral loss to follow-up', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Reassessment at 7 days loss to follow-up
# Assuming ~20% based on previous analyses
(p <- getAB(0.2, ((0-50)^2)/392^2))
png(here("plots/ltfu7d.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Reassessment at 7 days loss to follow-up', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

## Bacteriological test results
# Proportion of bacteriologically confirmable TB cases
# Using 17% (95% CI 10% to 25%) based on https://doi.org/10.1002/14651858.CD013693.pub2
## bacteriologically confirmable
curve(dbeta(x,17,40-17),from=0,to=1,n=200)

(p <- getAB(0.20,((50 - 0)/392)^2))
curve(dbeta(x,p$a,p$b),from=0,to=1,n=200)

# pre-treatment loss to follow-up
# 2-7% based on IHVN ATT initiation data
(p <- getAB(0.05,((0 - 8)/392)^2))
curve(dbeta(x,p$a,p$b),from=0,to=1,n=200)
png(here("plots/ptltfu.png"), width = 800, height = 600)
curve(dbeta(x, p$a, p$b), from = 0, to = 1, n = 200, 
      xlab = 'Pre-treatment loss to follow-up', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rbeta(1000, p$a, p$b))

# Case fatality ratios
# Based on Jenkins et al., 2017

# 0-4 years on treatment: 0.019 (0.012 - 0.029)
(p <- getLNparms(0.019,((2.9 - 1.2)/392)^2))
curve(dlnorm(x,p$mu,p$sig),from=0,to=1,n=200)
png(here("plots/cfrontxY.png"), width = 800, height = 600)
curve(dlnorm(x, -3.96, 0.64), from = 0, to = 1, n = 200, 
      xlab = 'CFR children <5 on TB treatment', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rlnorm(1000, p$mu, p$sig))

# 0-4 years without treatment: 0.436 (0.413 - 0.460)
(p <- getLNparms(0.436,((46 - 41.3)/392)^2))
curve(dlnorm(x,p$mu,p$sig),from=0,to=1,n=200)
png(here("plots/cfrnotxY.png"), width = 800, height = 600)
curve(dlnorm(x, -0.83, 0.08), from = 0, to = 1, n = 200, 
      xlab = 'CFR children <5 without TB treatment', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rlnorm(1000, p$mu, p$sig))

# 5-14 years on treatment: 0.008 (0.006 - 0.011)
(p <- getLNparms(0.008,((1.1 - 0.6)/392)^2))
curve(dlnorm(x,p$mu,p$sig),from=0,to=1,n=200)
png(here("plots/cfrontxO.png"), width = 800, height = 600)
curve(dlnorm(x, -4.82, 0.48), from = 0, to = 1, n = 200, 
      xlab = 'CFR children 5-14 on TB treatment', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rlnorm(1000, p$mu, p$sig))

# 5-14 years without treatment: 0.149 (0.137 - 0.162) 
(p <- getLNparms(0.149,((16.2 - 13.7)/392)^2))
curve(dlnorm(x,p$mu,p$sig),from=0,to=1,n=200)
png(here("plots/cfrnotxO.png"), width = 800, height = 600)
curve(dlnorm(x, -1.90, 0.12), from = 0, to = 1, n = 200, 
      xlab = 'CFR children 5-14 without TB treatment', 
      # xlab = expression(atop("TB presumption in children", "with active TB disease")),
      ylab = 'Density', 
      main = '')
dev.off()
summary(rlnorm(1000, p$mu, p$sig))

## things based on WHO data
E <- fread(here('who_data/data2023/TB_burden_age_sex_2024-07-23.csv'))

## NGA
NGA <- rbind(E[iso3=='NGA' & age_group=='0-4'],
             E[iso3=='NGA' & age_group=='5-14'])
NGA[,V:=((hi-lo)/3.92)^2]
NGA <- NGA[,.(mid=sum(best),V=sum(V)),by=age_group]

N <- fread(here('who_data/data2023/TB_notifications_2024-07-23.csv'))
N[iso3=='NGA' & year==max(N$year),.(newrel_f04,newrel_f1014,newrel_f514,newrel_f014)]

## NGA
NGA2 <- N[iso3=='NGA' & year==max(N$year),
          .(N.04=newrel_f04+newrel_m04,N.514=newrel_f514+newrel_m514)]

## NGA join
NGA2 <- melt(NGA2)
NGA2[,age_group:=c('0-4','5-14')]
NGA <- merge(NGA,NGA2,by='age_group')

## NGA calc
NGA[,cdr:=value/mid]
NGA[,cdr.hi:=value/(mid-sqrt(V))] #using 1 SD
NGA[,cdr.lo:=value/(mid+sqrt(V))]
NGA[,cdr.V:=(cdr.hi-cdr.lo)^2]
NGA[,c('a','b'):=getAB(cdr,cdr.V)] #TODO too low

curve(dbeta(x,NGA$a[1],NGA$b[1]),from=0,to=1,n=200)
curve(dbeta(x,NGA$a[2],NGA$b[2]),from=0,to=1,n=200)

NGA[,paste(a,b,sep=',')]

## fractions under 5
NGA[,frac:=mid/sum(mid)]

## simulate

## NGA
P <- getLNparms(NGA$mid,NGA$V)
yz <- rlnorm(1e4,P$mu[1],P$sig[1])
nyz <- rlnorm(1e4,P$mu[2],P$sig[2])
totz <- yz+nyz
qplot(yz/totz)
F <- getAB(mean(yz/totz),var(yz/totz))
paste(F$a,F$b,sep=',')

## read common data
PD <- read.csv(here('indata/parms.csv'))
P <- parse.parmtable(PD,testdir = here('test'),
                     outfile = here('test/00parms.csv'))
save(P,file=here(here('test/P.Rdata')))

