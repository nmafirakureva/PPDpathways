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
  ##sensitivity analyses flags:
  ## '' = basecase
  ## 'AllattendNHS' = Assumes 100% NHS referral attendance. The rest of the SA are built on this assumption
  ## 'noltfu'= 100% GP assessment, clinical suspicion, NHS attendance, starting & completing ATT,  starting & completing TPT
  ## 'FUVisitsCost' = No follow-up visit costs for both ATT and TPT
  ## 'DOTsCost' = No DOTs costs for both ATT and TPT
  ## 'XrayCost' = No X-ray costs
  ## 'PrisonEscort' = No prison escort costs
  ## 'ContactTracing' = No contact tracing costs
  ## 'InpatientCost' = No inpatient costs
  sacases <- c('', 'AllattendNHS','noltfu', 'FUVisitsCost', 'DOTsCost', 'XrayCost', 'PrisonEscort', 'ContactTracing', 'InpatientCost')
  SA <- sacases[1]
}

# rm(list=ls())
library(here)
library(tidyverse)

## load other scripts
source(here('R/ppd_pathways_tree.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/ppd_pathways_functions.R'))      #functions for tree parameters

# SA<-'noltfu' # TOPDO: remober after testing quick way of switching on sensitivity analysis
## number of reps
nreps <- 1e3
set.seed(1234)

## attributes to use
tblevels <- c('TBD','TBI', 'noTB') # TB disease, TB infection, no TB
agelevels <- c('15-100')
isoz <- c('GBR') #relevant countries

## prior parameters
PD <- read.csv(here('indata/ProbParms.csv')) # read in probability parameters
PS <- read.csv(here('indata/ProbParmsFixed.csv')) #read in fixed parameters
RD <- fread(gh('indata/RUParms.csv'))    # read resource use data
CD <- fread(gh('indata/CostParms.csv'))    # read cost data

# pre-process parameters
names(PD)
names(PS)
names(RD)
names(CD)

unique(PD$ParameterName)
names <- c("ParameterName", "Mean", "Range","Description","Source","SourceFull")
names(PS) <- names(RD) <- names(CD) <- names

PD.SOC <- PD |> 
  mutate(NAME = case_when(
    !grepl('tb.presum|started.att|ltbi|verbal|prog.tb|nContacts|sens.|spec.|pIsolation|DurMDRTB|Incomp|DurDSTB|pDSTB|smear', NAME) ~ paste0('soc.', NAME),
    .default = NAME
  ))

PD.INT <- PD |> 
  mutate(NAME = case_when(
    !grepl('tb.presum|started.att|ltbi|verbal|prog.tb|nContacts|sens.|spec.|pIsolation|DurMDRTB|Incomp|DurDSTB|pDSTB|smear', NAME) ~ paste0('int.', NAME),
    .default = NAME
  ))

unique(PD.INT$NAME)
PD0 <- rbind(PD.SOC, PD.INT)

# Fixed parameters to wide format
PS.SOC <- PS |> 
  mutate(ParameterName = paste0('soc.', ParameterName))

PS.INT <- PS |> 
  mutate(ParameterName = paste0('int.', ParameterName))

unique(PD.INT$NAME)
PS <- rbind(PS.SOC, PS.INT)

PD1 <- PS |> 
  select(ParameterName, Mean) |>
  distinct() |> 
  pivot_wider(names_from = ParameterName, values_from = Mean) 

PD1 <- rbind(
  PS |> 
    select(ParameterName, Mean),
  RD |> 
    select(ParameterName, Mean) |> 
    filter(!ParameterName %in% unique(PD0$NAME))
) |>
  distinct() |> 
  pivot_wider(names_from = ParameterName, values_from = Mean) 


P <- PD0 |> 
  select(NAME, DISTRIBUTION) |> 
  parse.parmtable()

# names(soc)

## make base PSA dataset
set.seed(1234) #random number seed

D0 <- makePSA(nreps,P)

# some checks
summary(D0)

## check if any probabilities are > 1, 
# Filter variables in prob_vrz that have values > 1
prob_vrz <- names(P)
impossible_values <- sapply(D0[, ..prob_vrz], function(x) any(x > 1))
filtered_cols <- prob_vrz[impossible_values] # Note: these are Okay to be > 1

# Summary of filtered columns
# These are okay
summary(D0[, ..filtered_cols])

## some probabilities are > 1, `quick fix` set to 1
# TODO: check for better distribution assumptioms
# D0[, (prob_vrz) := lapply(.SD, function(x) ifelse(x > 1, 1, x)), .SDcols = prob_vrz] # TODO: looks like a bug

## use these parameters to construct input data by attribute
D0 <- makeAttributes(D0)
D0[,sum(value),by=.(isoz,id)] #CHECK
D0[,sum(value),by=.(id, tb)] #CHECK

# merge in fixed parameters
D0 <- cbind(D0, PD1)

## read and make cost data
rcsts <- CD

names(CD)

rcsts <- rcsts |> 
  mutate(ParameterName = paste0('u', ParameterName),
         Mean = as.numeric(Mean), Range = as.character(Range),
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

D[,.(ucost.att.dots)]


names(D)[grepl('cost', names(D))]

# defining a range of scenarios
# TODO: check if these should be before makeTreeParms or after
if(SA == 'FUVisitsCost'){
  # No followup visit costs
  D[,int.prop.attend.nhs.referral:=1]
  D[,dstb.visits:=0]
  D[,mdrtb.visits:=0]
  D[,TPT.visits:=0]
}

if(SA == 'DOTsCost'){
  # No DOTs costs
  D[,int.prop.attend.nhs.referral:=1]
  D[,ucost.att.dots:=0]
}


if(SA == 'PrisonEscort'){
  # No followup visit costs
  D[,int.prop.attend.nhs.referral:=1]
  D[,ucost.prison.escort:=0]
  D[,mdrtb.visits:=0]
}

if(SA == 'XrayCost'){
  # No DOTs costs
  D[,int.prop.attend.nhs.referral:=1]
  D[,pxray:=0]
}

if(SA == 'InpatientCost'){
  # No inpatient costs
  D[,int.prop.attend.nhs.referral:=1]
  D[,DurDSTBIsolation:=0]
  D[,DurMDRTBIsolation:=0]
  D[,DurMDRTB2Isolation:=0]
}

## compute other parameters (adds by side-effect)
MakeTreeParms(D,P)

if(SA == 'noltfu'){
  # TB symptom screening
  # D[,int.prop.tb.sympt.screen:=1]
  # TB symptoms at screening
  # D[,int.prop.presumtive.tb:=ifelse(tb=='TBD',1,1-spec.symptom)]
  # Prison GP assessment
  D[,int.prop.prison.gp.assessment:=1]
  # D[,soc.prop.prison.gp.assessment:=1]
  # Clinical suspicion of TB disease
  # D[,soc.prop.clinical.tb.suspicion:=ifelse(tb=='TBD',1,1-spec.symptom)]
  # D[,int.prop.clinical.tb.suspicion:=ifelse(tb=='TBD',1,1-spec.any.abn.xray)]
  # Attending NHS referral
  D[,int.prop.attend.nhs.referral:=1]
  # starting & completing ATT
  D[,int.prop.starting.att:=1]
  # D[,int.prop.completing.att:=1]
  # # starting & completing TPT
  D[,int.prop.starting.tpt:=1]
  D[,int.prop.completing.tpt:=1]
} 


if(SA == 'AllattendNHS'){
  # Attending NHS referral
  D[,int.prop.attend.nhs.referral:=1]
}


if(SA == 'ContactTracing'){
  # No contact management costs
  D[,int.prop.attend.nhs.referral:=1]
  D[,cost.contact.management:=0]
}

## checks
D[,sum(value),by=.(isoz,id)] #CHECK
# D[,sum(value),by=.(isoz,id,age)] #CHECK

D[,.(isoz,age,soc.prop.tb.sympt.screen, int.prop.tb.sympt.screen)]
D[,.(ucost.dots)]

## check for leaks
head(SOC.F$checkfun(D)) #SOC arm
head(INT.F$checkfun(D)) #INT arm

names(SOC.F)

## === RUN MODEL
arms <- c('SOC','INT')
D <- runallfuns(D,arm=arms)                      #appends anwers

## restricted trees:
D[['soc_att_check']] <- SOC.att.F$checkfun(D)
D[['soc_att_cost']] <- SOC.att.F$costfun(D)
D[['soc_att_ppd']] <- SOC.att.F$cost.ppdfun(D)
D[['soc_att_nhs']] <- SOC.att.F$cost.nhsfun(D)
D[['soc_tpt_check']] <- SOC.tpt.F$checkfun(D)
D[['soc_tpt_cost']] <- SOC.tpt.F$costfun(D)
D[['soc_tpt_ppd']] <- SOC.tpt.F$cost.ppdfun(D)
D[['soc_tpt_nhs']] <- SOC.tpt.F$cost.nhsfun(D)
D[['soc_notx_check']] <- SOC.notx.F$checkfun(D)
D[['soc_notx_cost']] <- SOC.notx.F$costfun(D)

D[['int_att_check']] <- INT.att.F$checkfun(D)
D[['int_att_cost']] <- INT.att.F$costfun(D)
D[['int_att_ppd']] <- INT.att.F$cost.ppdfun(D)
D[['int_att_nhs']] <- INT.att.F$cost.nhsfun(D)
D[['int_tpt_check']] <- INT.tpt.F$checkfun(D)
D[['int_tpt_cost']] <- INT.tpt.F$costfun(D)
D[['int_tpt_ppd']] <- INT.tpt.F$cost.ppdfun(D)
D[['int_tpt_nhs']] <- INT.att.F$cost.nhsfun(D)
D[['int_notx_check']] <- INT.notx.F$checkfun(D)
D[['int_notx_cost']] <- INT.notx.F$costfun(D)

## cross checks: compute probability of endpoints from whole tree vs pruned trees
head(D[,int_att_check])
head(D[, attend.int])

## NOTE OK
all(D[, attend.int] == D[, int_att_check])
all(D[, attend.soc] == D[, soc_att_check])
all(D[, tptend.int] == D[, int_tpt_check])
all(D[, tptend.soc] == D[, soc_tpt_check])
## NOTE confusingly there is also a notxend variable, which is different
all(D[, notx.int] == D[, int_notx_check])
all(D[, notx.soc] == D[, soc_notx_check])

D[,table(tb)]

## create restricted PSA
DR <- D[,.(id,tb,
           soc_att_check,
           soc_att_cost,
           soc_att_ppdcost=soc_att_ppd/soc_att_cost,
           soc_tpt_check,
           soc_tpt_cost,
           soc_tpt_ppdcost=soc_tpt_ppd/soc_tpt_cost,
           soc_notx_check,
           soc_notx_cost,
           int_att_check,
           int_att_cost,
           int_att_ppdcost=int_att_ppd/int_att_cost,
           int_tpt_check,
           int_tpt_cost,
           int_tpt_ppdcost=int_tpt_ppd/int_tpt_cost,
           int_notx_check,
           int_notx_cost
)]

## condition costs on outcome
DR[,c('soc_att_cost',
      'soc_tpt_cost',
      'soc_notx_cost',
      'int_att_cost',
      'int_tpt_cost',
      'int_notx_cost'):=.(
        soc_att_cost/soc_att_check,
        soc_tpt_cost/soc_tpt_check,
        soc_notx_cost/soc_notx_check,
        int_att_cost/int_att_check,
        int_tpt_cost/int_tpt_check,
        int_notx_cost/int_notx_check
      )]

fn1 <- glue(here('outdata/DR')) + SA + '.Rdata'
save(DR,file=fn1)

## summary
DRS <- DR[,lapply(.SD,mean),.SDcols=names(DR)[-c(1,2)],by=tb]
DRS <- melt(DRS,id='tb')
DRS[,c('arm','outcome','quantity'):=tstrsplit(variable,split='_')]
options(scipen=999)
(DRS <- dcast(data=DRS,formula=arm+quantity+outcome~tb,value.var='value'))
fn1 <- glue(here('outdata/DRS')) + SA + '.csv'
fwrite(DRS,file=fn1)


summary(D[,.(cost.soc, cost.int)])

## cost of getting ATT (from CSV output)
# D[,mean(
#   pDSTB*dstb.visits*(ucost.dstb.opd.visit + ucost.prison.escort) + # DSTB visits
#     (1-pDSTB)*mdrtb.visits*(ucost.mdrtb.opd.visit + ucost.prison.escort) + # MDRTB visits
#     pDSTB*DurDSTB*ucost.dsatt.drugs + (1-pDSTB)*DurMDRTB*ucost.mdratt.drugs + # ATT drugs
#           ucost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) + # DOTS
#           cost.inpatient),
#   by=tb]
# 
# D[,mean(
#   IncompDurDSTB/DurDSTB*(pDSTB*dstb.visits*(ucost.dstb.opd.visit + ucost.prison.escort) +
#                            pDSTB*DurDSTB*(ucost.dsatt.drugs + ucost.dots)) +
#     IncompDurMDRTB/DurMDRTB*((1-pDSTB)*mdrtb.visits*(ucost.mdrtb.opd.visit + ucost.prison.escort)  +
#                                (1-pDSTB)*DurMDRTB*(ucost.mdratt.drugs + ucost.dots)) +
#     cost.inpatient
# ),by=tb]
# 
# D[,mean(durTPT*(ucost.ltbi.drugs + ucost.dots)),by=tb]
# D[,mean(durTPT*(ucost.ltbi.drugs + ucost.dots) + TPT.visits*(ucost.tpt.opd.visit + ucost.prison.escort)),by=tb]
# D[,mean(IncompDurTPT*(durTPT*(ucost.ltbi.drugs + ucost.dots) + TPT.visits*(ucost.tpt.opd.visit + ucost.prison.escort))),by=tb] # BUG fixed
# D[,mean(IncompDurTPT/durTPT*(durTPT*(ucost.ltbi.drugs + ucost.dots) + TPT.visits*(ucost.tpt.opd.visit + ucost.prison.escort))),by=tb]

