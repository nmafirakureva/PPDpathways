## Pete version
## quick run through the PPD pathways tree

## flags for sensitivity analyses
shell <- FALSE # whether running from shell script or not
if(shell){
  ## running from shell
  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  SA <- args[1]                  # none,base/lo/hi,tptru,hicoprev
  if(SA == 'none'){
    SA <- ''
  } 
} else { #set by hand
  rm(list=ls()) #clear all
  shell <- FALSE #whether running from shell script or not
  ##sensitivity analyses (mostly for PT):
  ## '' = basecase
  ## 'discr'='base'/'lo'/'hi'
  ## 'cdr' = making cdr higher for incidence
  ## 'txd' = making the completion influence tx/pt outcome
  sacases <- c('','lo','tptru','hicoprev', 'ctryeff','ugaattcsts', 'cdr', 'hivprev')
  SA <- sacases[1]
}

# rm(list=ls())
library(here)
library(tidyverse)

## load other scripts
source(here('R/ppd_active_tb_pathway1.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/ppd_pathways_functions2.R'))      #functions for tree parameters

## number of reps
nreps <- 1e3
set.seed(1234)

## attributes to use
tblevels <- c('TBD','TBI', 'noTB') # TB disease, TB infection, no TB
# ltbilevels <- c('LTBI','noLTBI') #LTBI, not LTBI
# agelevels <- c('15-64','65-100')
agelevels <- c('15-100')
isoz <- c('GBR') #relevant countries

## --- life years and other outputs NOTE needs to be set FALSE on first run thru
LYSdone <- TRUE
if(!LYSdone){
  ## make discounted life-years if they haven't been done
  LYKc <- GetLifeYears(isolist=isoz,discount.rate=0.03,yearfrom=2021)
  LYKc0 <- GetLifeYears(isolist=isoz,discount.rate=0.00,yearfrom=2021)
  LYKc5 <- GetLifeYears(isolist=isoz,discount.rate=0.05,yearfrom=2021)
  LYKc <- merge(LYKc,LYKc0[,.(iso3,age,LYS0=LYS)],by=c('iso3','age'))
  LYKc <- merge(LYKc,LYKc5[,.(iso3,age,LYS5=LYS)],by=c('iso3','age'))
  LYK <- LYKc[,.(LYS=mean(LYS),LYS0=mean(LYS0),LYS5=mean(LYS5)),by=.(age)] #averaged life-years 4 generic tests
  save(LYKc,file=here('indata/LYKc.Rdata'))
  save(LYK,file=here('indata/LYK.Rdata'))
} else {
  load(file=here('indata/LYKc.Rdata'))
  load(file=here('indata/LYK.Rdata'))
}

# Sensitivity analysis: 0% & 5% discount rates
if(SA %in% c('hi','lo')){
  LYKc[,LYS:=ifelse(SA=='lo', LYS0, 
                    ifelse(SA=='hi',LYS5, LYS))]  
}

## prior parameters
PD <- read.csv(here('indata/CSV/ProbParms2.csv')) #read in probability parameters
AD <- read.csv(here('indata/DiagnosticAccuracy.csv')) #read in accuracy parameters
RD <- fread(gh('indata/RUParms.csv'))    #read resource use data
CD <- fread(gh('indata/CostParms.csv'))    #read cost data

names(PD)
names(RD)
names(CD)

names <- c("ParameterName", "Mean", "Range","Year","Description","Source")
names(PD) <- names(RD) <- names(CD) <- names

PD1 <- PD |> 
  filter(ParameterName != 'tb.prev') |>
  mutate(ParameterName = paste0('soc.', ParameterName))

PD2 <- PD |> 
  filter(ParameterName != 'tb.prev') |>
  mutate(ParameterName = paste0('int.', ParameterName))

PD3 <- PD |> 
  filter(ParameterName %in% c('tb.prev', 'ltbi.prev', 'prog.tb', 'progInf', 'progNInf')) 

PD0 <- rbind(PD1, PD2, PD3, RD)

PD0$NAME <- PD0$ParameterName

PD0 <- PD0 |> 
  mutate(Mean = as.numeric(Mean), Range = as.character(Range),
         Median = ifelse(Range!='', paste0(Mean, ' (', Range, ')'), Range),
         Median = ifelse(!is.na(as.numeric(Range)), '', Median)) 

PD1 <- PD0 |> 
  filter(Median != '') # parameters to be sampled

PD2 <- PD0 |> 
  filter(Median == '') # parameters to be fixed (no uncertainty estimates)

tmp <- PD1 %>%
  extract(Median, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") %>%
  separate(lo,c("lo","hi"),"-") %>%
  mutate_at(c("mid", "lo","hi"), as.numeric)

tmp <- setDT(tmp)
tmp1 <- getLNparms(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2,med=FALSE)
tmp[,DISTRIBUTION:=paste0("LN(",tmp1$mu,",",tmp1$sig,")")] #LN distributions

# tmp1 <- getAB(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2)
# tmp[,DISTRIBUTION:=paste0("B(",tmp1$a,",",tmp1$b,")")] # Beta distributions
# tmp[,DISTRIBUTION:=paste0("G(",tmp[,mid]^2/((tmp[,hi]-tmp[,lo])^2/3.92^2),",",((tmp[,hi]-tmp[,lo])^2/3.92^2)/tmp[,mid],")")] # Gamma distributions
PD1 <- PD1 |> 
  mutate(NAME=ParameterName,
         DISTRIBUTION=tmp$DISTRIBUTION,
         # DISTRIBUTION=ifelse(DISTRIBUTION=='LN(NA,NA)', NA, DISTRIBUTION))
         DISTRIBUTION=ifelse(DISTRIBUTION=='B(NA,NA)', NA, DISTRIBUTION))

## 
PD1 <- PD1 |> 
  filter(DISTRIBUTION!="")

# Fixed parameters to wide format
PD3 <- PD2 |> 
  select(ParameterName, Mean) |>
  distinct() |> 
  pivot_wider(names_from = ParameterName, values_from = Mean) 

# Diagnostic accuracy 
tmp <- AD %>%
  extract(mqrng, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") %>%
  separate(lo,c("lo","hi"),"to") %>%
  mutate_at(c("mid", "lo","hi"), as.numeric)

tmp <- setDT(tmp)
tmp1 <- getAB(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2)
tmp[,DISTRIBUTION:=paste0("B(",tmp1$a,",",tmp1$b,")")] # Beta distributions

AD1 <- tmp |> 
  filter(grepl('symptom|any.abn.xray|xpert', NAME)) |>
  select(NAME, DISTRIBUTION) |> 
  as.data.frame() 


# convert into parameter object
P <- rbind(
  PD1 |> 
  select(NAME, DISTRIBUTION),
  AD1) |> 
  parse.parmtable(outfile='out.csv',
                  testdir = here('plots/test'))             

names(P)

## make base PSA dataset
set.seed(1234) #random number seed

D0 <- makePSA(nreps,P)

# some checks
prob_vrz <- names(P)
prob_vrz <- prob_vrz[!grepl('verbal_screen_time', prob_vrz)]
summary(D0[,..prob_vrz])

## some probabilities are > 1, `quick fix` set to 1
# TODO: check for better distribution assumptioms
D0[, (prob_vrz) := lapply(.SD, function(x) ifelse(x > 1, 1, x)), .SDcols = prob_vrz]

# Filter variables in prob_vrz that have values > 1
impossible_values <- sapply(D0[, ..prob_vrz], function(x) any(x > 1))
filtered_cols <- prob_vrz[impossible_values]

# Summary of filtered columns
summary(D0[, ..filtered_cols])

summary(D0[,.(ltbi.prev, prog.tb, progInf, progNInf)])

## use these parameters to construct input data by attribute
D0 <- makeAttributes(D0)
D0[,sum(value),by=.(isoz,id)] #CHECK
D0[,sum(value),by=.(id, tb)] #CHECK

# merge in fixed parameters
D0 <- cbind(D0, PD3)
  
## read and make cost data
rcsts <- CD

names(CD)

rcsts <- rcsts |> 
  mutate(Mean = as.numeric(Mean), Range = as.character(Range),
         Median = ifelse(Range!='', paste0(Mean, ' (', Range, ')'), Range))

rcsts <- rcsts %>%
  extract(Median, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") %>%
  separate(lo,c("lo","hi"),"-") %>%
  mutate_at(c("mid", "lo","hi"), as.numeric)

rcsts <- rcsts |> 
  mutate(cost.m = Mean,
         cost.sd = (hi-lo)/3.92)

rcsts <- setDT(rcsts)

## turn cost data into PSA
rcsts[is.na(rcsts)] <- 0 # some quick fix >> setting NA to 0
rcsts[cost.sd==0,cost.sd:=cost.m/40]        #SD such that 95% UI ~ 10% of mean

allcosts <- rcsts[,.(iso3=isoz, cost=ParameterName, cost.m, cost.sd)]

C <- MakeCostData(allcosts[iso3=='GBR'],nreps)               # make cost PSA NOTE using CMR cost data

## NOTE can re-run from here to implement changes to MakeTreeParms
## add cost data
D <- merge(D0,C,by=c('id'),all.x=TRUE)       # merge into PSA (differentiated D and D0 to facilitate rerunning)

## compute other parameters (adds by side-effect)
MakeTreeParms(D,P)


names(D)[grepl('cost', names(D))]

## checks
D[,sum(value),by=.(isoz,id)] #CHECK
# D[,sum(value),by=.(isoz,id,age)] #CHECK

D[,.(isoz,age,soc.prop.tb.sympt.screen, int.prop.tb.sympt.screen)]

## check for leaks
head(SOC_ATB.F$checkfun(D)) #SOC arm
head(INT_ATB.F$checkfun(D)) #INT arm

names(SOC_ATB.F)

## === RUN MODEL
arms <- c('SOC_ATB','INT_ATB')
D <- runallfuns(D,arm=arms)                      #appends anwers

## restricted trees:
D[['soc_att_check']] <- SOC.att.F$checkfun(D)
D[['soc_att_cost']] <- SOC.att.F$costfun(D)
# D[['soc_tpt_check']] <- SOC.tpt.F$checkfun(D)
# D[['soc_tpt_cost']] <- SOC.tpt.F$costfun(D)
D[['soc_notx_check']] <- SOC.notx.F$checkfun(D)
D[['soc_notx_cost']] <- SOC.notx.F$costfun(D)

D[['int_att_check']] <- INT.att.F$checkfun(D)
D[['int_att_cost']] <- INT.att.F$costfun(D)
# D[['int_tpt_check']] <- INT.tpt.F$checkfun(D)
# D[['int_tpt_cost']] <- INT.tpt.F$costfun(D)
D[['int_notx_check']] <- INT.notx.F$checkfun(D)
D[['int_notx_cost']] <- INT.notx.F$costfun(D)

## cross checks: compute probability of endpoints from whole tree vs pruned trees
head(D[,int_att_check])
head(D[, attend.int])

## NOTE OK
all(D[, attend.int] == D[, int_att_check])
all(D[, attend.soc] == D[, soc_att_check])
# all(D[, tptend.int] == D[, int_tpt_check])
# all(D[, tptend.soc] == D[, soc_tpt_check])
## NOTE confusingly there is also a notxend variable, which is different
all(D[, notx.int] == D[, int_notx_check])
all(D[, notx.soc] == D[, soc_notx_check])


D[,table(tb)]

## create restricted PSA
DRPCF <- D[,.(id,tb,
           soc_att_check,
           soc_att_cost,
           # soc_tpt_check,
           # soc_tpt_cost,
           soc_notx_check,
           soc_notx_cost,
           int_att_check,
           int_att_cost,
           # int_tpt_check,
           # int_tpt_cost,
           int_notx_check,
           int_notx_cost
           )]

## condition costs on outcome
DRPCF[,c('soc_att_cost',
      # 'soc_tpt_cost',
      'soc_notx_cost',
      'int_att_cost',
      # 'int_tpt_cost',
      'int_notx_cost'):=.(
        soc_att_cost/soc_att_check,
        # soc_tpt_cost/soc_tpt_check,
        soc_notx_cost/soc_notx_check,
        int_att_cost/int_att_check,
        # int_tpt_cost/int_tpt_check,
        int_notx_cost/int_notx_check
      )]

save(DRPCF,file=here('outdata/DRPCF.Rdata'))

## summary
DRSPCF <- DRPCF[,lapply(.SD,mean),.SDcols=names(DRPCF)[-c(1,2)],by=tb]
DRSPCF <- melt(DRSPCF,id='tb')
DRSPCF[,c('arm','outcome','quantity'):=tstrsplit(variable,split='_')]
(DRSPCF <- dcast(data=DRSPCF,formula=arm+quantity+outcome~tb,value.var='value'))
fwrite(DRSPCF,file=here('outdata/DRSPCF.csv'))

# D[tb=='noTB',.(soc.prop.prev.tb.dx,
#     soc.prop.no.prev.tb.dx.symp,
#     soc.prop.no.prev.tb.dx.symp.gp.assess,
#     soc.prop.no.prev.tb.dx.symp.tb.suspicion,
#     soc.prop.no.prev.tb.dx.symp.nhs.referral,
#     soc.prop.no.prev.tb.dx.symp.tb.dx,
#     soc.prop.starting.att,
#     soc.prop.completing.att,
#     sens.any.abn.xray,spec.any.abn.xray)]

## soc.prop.no.prev.tb.dx.symp.tb.dx??

## cost of getting ATT (from CSV output)
D[,mean((pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + (1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visit + cost.prison.escort)) + (pDSTB*DurDSTB*cost.dsatt.drugs + (1-pDSTB)*DurMDRTB*cost.mdratt.drugs) + cost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB)),by=tb]
