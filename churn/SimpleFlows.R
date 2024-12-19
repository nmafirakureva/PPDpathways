## see https://docs.google.com/presentation/d/1JFFa0SykXvwsf5IN1iXySOUKQyCo59DcTdyg0mwPMxw/edit?pli=1#slide=id.p

library(posterior)
library(corplot) #github: petedodd/corplot
library(rstan)
library(data.table)
library(here)

## === fraction open:
E <- fread(here("churn/indata/The_Prison_Estate_and_Probation_Service_Region_register_June_23.xlsx - Prison Estate.csv"),
           skip=1)
E[,unique(`Predominant Function`)]
E[,unique(`Cohort of Prisoners Held`)]
E[,unique(Designation)]
E[,unique(Notes)]

exclude <- c( "Reception","Female","Predominant Function","YJB","Predominant Function")
E <- E[!`Predominant Function` %in% exclude]
E <- E[Telephone!=''] #ditch extra lines
E[,Open:=ifelse(grepl('Open',`Predominant Function`),TRUE,FALSE)]
E <- E[,.(Prison,Open)]
E[,table(Open)]


P <- fread(here("churn/indata/Population_data_tool_31Dec2023.csv"))
P <- P[view=='a Establishment*Sex*Age Group']
P <- merge(P,E,by.x='Establishment',by.y='Prison')
P <- P[Date=="2023-12"]
P <- P[,.(pop=sum(Population,na.rm=TRUE)),by=Open]
P

## release data
load(file=here('churn/indata/RLM.rda'))
RLM <- RLM[year==2022]
## TODO release includes remand?
RLM[,md:=c(0.25,0.75,1.5,3,4.5,6,7.5,12,16)]
RLM[,weighted.mean(md,w=value)]


## -- parameters --
## E = remand inflow
## d R = mean durn remand
## d L = mean durn Long
## d S = mean durn Short
## d O = mean durn Open
## 𝛕 L = transfer rate from Long to Short
## 𝛕 S = transfer rate from Short to Open
## p S = fraction of remand not released
## p L = fraction of sentenced remand to Long



## -- evidence --
## Remand inflow (remand first acceptance) [E]
## Sentence inflow* [→p S ]
## Total population
## Remand population
## Fraction in open (shared key)
## Total outflow (releases*)
## Mean duration at release


## https://www.gov.uk/government/statistics/offender-management-statistics-quarterly-july-to-september-2023/offender-management-statistics-quarterly-july-to-september-2023

## totals
totpop <-  87489
totsentence <- 71042
totremand <- 16005

## receptions
fstrecep_remand <- 12747
fstrecep_sentence <- 6193

## releases
releases <- 12351

## fraction open
fracopen <- P[Open==TRUE,pop]/P[,sum(pop)]

## mean durn TODO
RLM[,weighted.mean(md,w=value)] #2.2


## mean durations:
## Long:   1 / (1/dL + 𝛕L)
## Short:  1 / (1/dS + 𝛕S)
## Open:   dO



## -- eqns --
## dR/dt = E - R/d R                                         --[1]
## dL/dt = p S p L R/d R - 𝛕 L L - L/d L                     --[2]
## dS/dt = p S (1-p L )R/d R + 𝛕 L L - S/d S - 𝛕 S S         --[3]
## dΩ/dt = 𝛕 S S - Ω/d O                                     --[4]


## getRemandDur <- function(E,R) R/E
## (dR <- getRemandDur(fstrecep_remand,totremand)) #1.3 yrs


## getPs <- function(E,sentenced) sentenced / E
## (Ps <- getPs(fstrecep_remand,fstrecep_sentence)) #~50%


## eqn2 <- function(pS,pL,dR,tL,dL,
##                  R,L)  pS*pL*R/dR-tL*L-L/dL

## eqn3 <- function(pS,pL,dR,tL,tS,dS,
##                  S,R,L)  pS*(1-pL)*R/dR+tL*L-S/dS-tS/S

## eqn4 <- function(tS,dO,
##                  S,Op)  tS/S - Op/dO



## ================================================
mdl <- stan_model(here('churn/stan/eqm1.stan'))

## totals
totpop <-  87489
totsentence <- 71042
totremand <- 16005

## receptions
fstrecep_remand <- 12747
fstrecep_sentence <- 6193

## releases
releases <- 12351

## fraction open
fracopen <- P[Open==TRUE,pop]/P[,sum(pop)]

## TODO outflow including remand or not?

## test
sdat <- list(
  ## data points: mean and sd
  E_m = fstrecep_remand,  E_s = fstrecep_remand/100,   #remand inflow
  R_m = totremand,  R_s = totremand/100,   #remand stock
  SI_m = fstrecep_sentence,  SI_s = fstrecep_sentence/100, #sentence inflow
  N_m = totpop,  N_s = totpop/100,   #total population
  FO_a = 100*fracopen,  FO_b = 100*(1-fracopen), #fraction open NOTE arbitrary precision
  TO_m = releases,  TO_s = releases/100, #total outflow
  ## mean duration at release 
  DR_m = RLM[,weighted.mean(md,w=value)], #2.2
  DR_s = RLM[,weighted.mean(md,w=value)]/20, 
  ## other priors
  dR_m = 0.5,  dR_s = .1, #duration of remand NOTE guess
  dS_m = 1, dS_s = 0.5, #short duration NOTE defn
  dL_m = 10,            #long NOTE defn
  dom_m = 1,
  ts_m = 0.1,  tl_m = 0.1,
  ps_a = 2,   ps_b = 2, #fraction sentenced
  pl_a = 2,   pl_b = 18, #fraction to long
  ## tolerance for errors
  tol = rep(0.05,9)
  ## tol = c(fstrecep_remand/100,
  ##         totremand/10,
  ##         totremand/10,
  ##         totremand/10,
  ##         totremand/10,
  ##         totpop/100,
  ##         releases/100,
  ##         releases/10, #mean duration * releases
  ##         totremand/10 #total open?
  ##         )
)
## TODO add remand pop as 10? or remove total pop?

##
samps <- sampling(mdl,data=sdat,iter=2000,chains=4,cores=4)

## pairs(samps)

## TODO improve:



## library(bayesplot)

smps <- as_draws_df(samps)
names(smps)
(dtmp <- summarise_draws(smps))
print(dtmp,n=Inf)

smpsd <- as.data.table(smps)
smpsd[,summary(N)]; totpop
smpsd[,summary(R)]; totremand

save(smpsd,file=here('churn/outdata/smpsd.Rdata'))



## make tables:
ttmp <- as.data.table(dtmp)

tab1 <- ttmp[variable %in% c('N','S','R','L','Omega'),
             .(variable,mean=signif(mean,2),sd=signif(sd,2))]
tab1[,name:=c('Total Pop','Remand','Short','Long','Open')]
tab2 <- ttmp[variable %in% c('E','SI','TO','DR','FO'),
             .(variable,mean=signif(mean,2),sd=signif(sd,2))]
tab2[,name:=c('inflow','sentenced','outflow','duration at release','fraction open')]
tab3 <- ttmp[variable %in% c('dR','dS','dL','dom','ts','tl','ps','pl'),
             .(variable,mean=signif(mean,2),sd=signif(sd,2))]
tab3[,name:=c('remand timescale','short timescale','long timescale','open timescale',
              'short transfer rate','long transfer rate','fraction sentenced','sentenced long')]


cls <- c('name','variable','mean','sd')
setcolorder(tab1,cls)
setcolorder(tab2,cls)
setcolorder(tab3,cls)

fwrite(tab1,file=here('churn/outdata/tab1.csv'))
fwrite(tab2,file=here('churn/outdata/tab2.csv'))
fwrite(tab3,file=here('churn/outdata/tab3.csv'))

## TODO neater

## rounding tables
## always round 0.5 up
round2 <- function(x, digits = 0) sign(x) * trunc(abs(x) * 10^digits + 0.5) / 10^digits
ft <- Vectorize(function(x){
  smallpos <- x > 0 & x < 0.01
  one2ten <- x >= 1 & x < 10
  zero2one <- x >= 0.1 & x < 1
  dg <- ifelse(abs(x) > 0.01 & abs(x) < 100, 2, 3)
  x2 <- signif(x, dg)
  trailing.0 <- x2 == round2(x) & one2ten == TRUE
  trailing0 <- x2 * 10 == round2(x * 10) & zero2one == TRUE & x2 < 1
  format(
    x2,
    digits = dg,
    nsmall = 0L,
    big.mark = " ",
    justify = 'right',
    drop0trailing = TRUE,
    scientific = FALSE
  )
})

brkt <- function(x,y,z) paste0(ft(x),' (',
                                        ft(y),' to ',
                                        ft(z),')')

ttmp <- ttmp[!variable %in% c(paste0('eqn[',1:9,']'),'LP__')]
ttmp <- ttmp[,.(variable,mean,q5,q95,rhat,ess_bulk)]

popnames <- c('N','E','R','S','L','Omega','SI','TO')
pops <- ttmp[variable %in% popnames]
nonpops <- ttmp[!variable %in% popnames]


## pops
emps <- subset_draws(as_draws_df(samps),variable=popnames)
emps[,c('.chain','.iteration','.draw')] <- NULL
corplot(emps)
corplot(emps,file=here('churn/outdata/p_pops.png'))


pops <- pops[,.(variable,value=brkt(mean,q5,q95),Rhat=round2(rhat,2),ESS=round2(ess_bulk))]

pops[variable=='TO',variable:='Released']
pops[variable=='SI',variable:='Sentenced']

pops
fwrite(pops,file=here('churn/outdata/pops.out.csv'))

## non pops
nonpops <- nonpops[,.(variable,value=brkt(mean,q5,q95),Rhat=round2(rhat,2),ESS=round2(ess_bulk))]

emps <- subset_draws(as_draws_df(samps),variable=nonpops$variable)
emps[,c('.chain','.iteration','.draw')] <- NULL
corplot(emps)
corplot(emps,file=here('churn/outdata/p_nonpops.png'))


nonpops
fwrite(nonpops,file=here('churn/outdata/nonpops.out.csv'))


## ## TODO check whether frac open includes remand

## 1

## ## smpsf <- smps #with the chain/iteration/draw included
## ## smpsf[,c('.chain','.iteration','.draw')] <- NULL

## ## head(smpsf)

## ## corplot(as.matrix(smpsf))


## emps <- subset_draws(as_draws_df(samps),variable=c('E','dR','R','N'))
## emps[,c('.chain','.iteration','.draw')] <- NULL
## corplot(emps)

## corplot(emps,file='p_eg.png')

## ## pops
## emps <- subset_draws(as_draws_df(samps),variable=c('E','R','L','S','Omega','N'))
## emps[,c('.chain','.iteration','.draw')] <- NULL
## corplot(emps)

## corplot(emps,file='p_pops.png')


## emps <- subset_draws(as_draws_df(samps),variable=c('dR','dS','dL','dom','ts','tl','ps','pl'))
## emps[,c('.chain','.iteration','.draw')] <- NULL
## corplot(emps,file='p_parms.png')

## ## corplot(emps,points=TRUE)

## nrow(emps)

## ## fn <- here('exploratory/figures/post.log.pdf')
## ## corplot(emps,file=fn)


## https://github.com/nmafirakureva/PPDpathways
