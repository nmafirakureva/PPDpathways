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
source(here('R/ppd_pathways_tree.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/ppd_pathways_functions.R'))      #functions for tree parameters

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
  filter(grepl('symptom|any.abn.xray', NAME)) |>
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

D <- makePSA(nreps,P)

# some checks
prob_vrz <- names(P)
prob_vrz <- prob_vrz[!grepl('verbal_screen_time', prob_vrz)]
summary(D[,..prob_vrz])

## some probabilities are > 1, `quick fix` set to 1
# TODO: check for better distribution assumptioms
D[, (prob_vrz) := lapply(.SD, function(x) ifelse(x > 1, 1, x)), .SDcols = prob_vrz]

# Filter variables in prob_vrz that have values > 1
impossible_values <- sapply(D[, ..prob_vrz], function(x) any(x > 1))
filtered_cols <- prob_vrz[impossible_values]

# Summary of filtered columns
summary(D[, ..filtered_cols])

summary(D[,.(ltbi.prev, prog.tb, progInf, progNInf)])

## use these parameters to construct input data by attribute
D <- makeAttributes(D)
D[,sum(value),by=.(isoz,id)] #CHECK
D[,sum(value),by=.(id, tb)] #CHECK

# merge in fixed parameters
D <- cbind(D, PD3)
  
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

## add cost data
D <- merge(D,C,by=c('id'),all.x=TRUE)       # merge into PSA

## compute other parameters (adds by side-effect)
MakeTreeParms(D,P)


names(D)[grepl('cost', names(D))]

## checks
D[,sum(value),by=.(isoz,id)] #CHECK
# D[,sum(value),by=.(isoz,id,age)] #CHECK

D[,.(isoz,age,soc.prop.tb.sympt.screen, int.prop.tb.sympt.screen)]

## check for leaks
head(SOC.F$checkfun(D)) #SOC arm
head(INT.F$checkfun(D)) #INT arm

names(SOC.F)

## === RUN MODEL
arms <- c('SOC','INT')

D <- runallfuns(D,arm=arms)                      #appends anwers

# summary(D)

names(D)[!grepl('cost', names(D))]

## --- run over different countries
cnmz <- names(C) # C=cost PSA data
cnmz <- cnmz[cnmz!=c('id')]
toget <- c('id', 'value',
           'cost.soc','cost.int',
           'cost.screen.soc',	'cost.screen.int',
           'cost.tpt.soc', 'cost.tpt.int',
           'cost.tb.assessment.soc','cost.tb.assessment.int',
           'cost.att.soc','cost.att.int',
           'cost.ppd.soc','cost.ppd.int',	
           'cost.nhs.soc','cost.nhs.int',
           'screen.soc','screen.int',
           'igra.soc','igra.int',
           'ltbi.soc','ltbi.int',
           'noltbi.soc','noltbi.int',
           'noltbinotpt.soc','noltbinotpt.int',
           'tpt.soc','tpt.int',
           'ltbinotpt.soc','ltbinotpt.int',
           'att.soc','att.int',
           'noatt.soc','noatt.int'
           )

toget2 <- c(toget,
            'coprevtb.soc','coprevtb.int')

tosum <- toget
tosum2 <- c(tosum,
            'coprevtb.soc','coprevtb.int')

## heuristic to scale top value for thresholds:
# heur <- c('id','value', 'coprevtb.soc','coprevtb.int')
# out <- D[,..heur]
# out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=c('coprevtb.soc','coprevtb.int'),by=id] #sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.int)]

topl <- 3e4
lz <- seq(from = 0,to=topl,length.out = 1000) #threshold vector for CEACs

## containers & loop
allout <- allpout <- list() #tabular outputs
allout2 <- allpout2 <- list() #tabular outputs
# ceacl <- NMB <- list()             #CEAC outputs etc
psaout <- psapout <- list()
parmsout <- parmsout <- list()

## NOTE I think there was an additional problem here -
## because countries are contained as separate rows we were summing over both in the below computations
## ie we were looping but then operating over data for both countries
## cn <- isoz[1]

for(cn in isoz){
  cn <- 'GBR'
  dc <- D[isoz==cn]
  cat('running model for:',cn, '15-100 years\n')
  ## --- costs
  ## drop previous costs
  # D <- D[age=='0-4',]
  dc[,c(cnmz):=NULL]
  ## add cost data
  C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
  dc <- merge(dc,C,c('id'),all.x=TRUE)        #merge into PSA

  ## --- run model (quietly)
  invisible(capture.output(dc <- runallfuns(dc,arm=arms)))
  ## --- grather outcomes
  out <- dc[,..toget]

  ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
  out <- out[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum,by=id] #sum against popn
  
  ## non-incremental cost per ATT
  out[,costperATT.soc:=cost.soc/att.soc];
  out[,costperATT.int:=cost.int/att.int]; 
  out[,costperTPT.soc:=cost.soc/tpt.soc];
  out[,costperTPT.int:=cost.int/tpt.int]; 
  
  # proportion of cost.ppd & cost.nhs to cost.soc or cost.int
  out[,frac.cost.ppd.soc:=cost.ppd.soc/cost.soc];
  out[,frac.cost.ppd.int:=cost.ppd.int/cost.int];
  out[,frac.cost.nhs.soc:=cost.nhs.soc/cost.soc];
  out[,frac.cost.nhs.int:=cost.nhs.int/cost.int];
  
  ## increments wrt SOC (per person screened)
  out[,Dcost.int:=cost.int-cost.soc];  #inc costs
  out[,Dcost.ppd.int:=cost.ppd.int-cost.ppd.soc];  #inc costs
  out[,Dcost.nhs.int:=cost.nhs.int-cost.nhs.soc];  #inc costs
  out[,Datt.int:=att.int-att.soc];  #inc atts
  out[,Dtpt.int:=tpt.int-tpt.soc];  #inc atts

  ## per whatever
  out[,DcostperATT.int:=Dcost.int/Datt.int];
  out[,DcostperTPT.int:=Dcost.int/Dtpt.int];
  
  ## summarize
  smy <- outsummary(out)
  outs <- smy$outs; pouts <- smy$pouts;
  outs[,iso3:=cn]; pouts[,iso3:=cn]
  
  ## capture tabular
  allout[[cn]] <- outs; allpout[[cn]] <- pouts
  
  ## --- grather outcomes for Table 2
  out2 <- dc[,..toget2]
  
  ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
  out2 <- out2[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum2,by=id] #sum against popn
  
  ## increments
  out2[,Dscreen:=screen.int-screen.soc]
  out2[,Dtpt:=tpt.int-tpt.soc]
  out2[,Datt:=att.int-att.soc]
  out2[,Dcoprevtb:=coprevtb.int-coprevtb.soc]
  
  out2[,Dcost.screen:=cost.screen.int-cost.screen.soc]
  out2[,Dcost.tpt:=cost.tpt.int-cost.tpt.soc]
  out2[,Dcost.att:=cost.att.int-cost.att.soc]
  out2[,Dcost:=cost.int-cost.soc]
  
  ## summarize
  smy2 <- Table2(out2) #NOTE set per 1000 new residents - adjust fac in ppd_pathways_functions.R
  outs2 <- smy2$outs; pouts2 <- smy2$pouts;
  outs2[,iso3:=cn]; pouts2[,iso3:=cn]
  psapout[[cn]] <- out2[,iso3:=cn]
  ## capture tabular
  allout2[[cn]] <- outs2; allpout2[[cn]] <- pouts2

  psaout[[cn]] <- dc[,.(iso3=cn,
                        screen.soc=sum(screen.soc*value),
                        screen.int=sum(screen.int*value),
                        ltbi.soc=sum(ltbi.soc*value),
                        ltbi.int=sum(ltbi.int*value),
                        noltbi.soc=sum(noltbi.soc*value),
                        noltbi.int=sum(noltbi.int*value),
                        noltbinotpt.soc=sum(noltbinotpt.soc*value),
                        noltbinotpt.int=sum(noltbinotpt.int*value),
                        tpt.soc=sum(tpt.soc*value),
                        tpt.int=sum(tpt.int*value),
                        ltbinotpt.soc=sum(ltbinotpt.soc*value),
                        ltbinotpt.int=sum(ltbinotpt.int*value),
                        coprevtb.soc=sum(coprevtb.soc*value),
                        coprevtb.int=sum(coprevtb.int*value),
                        att.soc=sum(att.soc*value),
                        att.int=sum(att.int*value),
                        noatt.soc=sum(noatt.soc*value),
                        noatt.int=sum(noatt.int*value),
                        cost.screen.soc=sum(cost.screen.soc*value),
                        cost.screen.int=sum(cost.screen.int*value),
                        cost.tpt.soc=sum(cost.tpt.soc*value),
                        cost.tpt.int=sum(cost.tpt.int*value),
                        cost.tb.assessment.soc=sum(cost.tb.assessment.soc*value),
                        cost.tb.assessment.int=sum(cost.tb.assessment.int*value),
                        cost.att.soc=sum(cost.att.soc*value),
                        cost.att.int=sum(cost.att.int*value),
                        cost.ppd.soc=sum(cost.ppd.soc*value),
                        cost.ppd.int=sum(cost.ppd.int*value),
                        cost.nhs.soc=sum(cost.nhs.soc*value),
                        cost.nhs.int=sum(cost.nhs.int*value),
                        # frac.cost.ppd.soc=sum(frac.cost.ppd.soc*value),
                        # frac.cost.ppd.int=sum(frac.cost.ppd.int*value),
                        cost.soc=sum(cost.soc*value),
                        cost.int=sum(cost.int*value)),
                     by=.(id, isoz, age, tb)] #PSA summary
  
  psaout2 <- dc[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum[-1],by=c('id', 'tb')] #sum against popn
  
  # 
  # 
  prop_parms <- vrz[grepl('prop', vrz)]
  acc_parms <- names(dc)[grepl('sens.|spec.', names(dc))]
  prev_parms <- names(dc)[grepl('ltbi.prev|prog.', names(dc))]
  ru_parms <- vrz[!grepl('prop|cost', vrz)]
  cost_parms <- vrz[grepl('cost', vrz)]

  # 
  vars <- c(prop_parms, acc_parms, prev_parms, ru_parms, cost_parms)
  parms <- dc[,lapply(.SD,function(x) mean(x, na.rm = T)),.SDcols=vars,by=.(id, isoz, age)]
  
  parms[,iso3:=cn]
  parms[,isoz:=NULL]
  parmsout[[cn]] <- parms
}

summary(D[,..toget])

# 
allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
allout2 <- rbindlist(allout2)
allpout2 <- rbindlist(allpout2)
psaout <- rbindlist(psaout)
# psapout <- rbindlist(psapout)
parmsout <- rbindlist(parmsout)

# quick look at other things
soc <- names(allpout)[grepl('.soc', names(allpout))]
allpout |> 
  select(all_of(soc)) 

int <- names(allpout)[grepl('.int', names(allpout))]
int <- int[!grepl('^D', int)]
allpout |> 
  select(all_of(int))

# #
# allpout2 |> 
#   select(screen.soc, screen.int, att.soc, att.int, Datt, tpt.soc, tpt.int, Dtpt)
# 
# # costs
# allpout2 |> 
#   select(cost.screen.soc, cost.screen.int, cost.att.soc, cost.att.int, Dcost.att, cost.tpt.soc, cost.tpt.int, Dcost.tpt)
# 
# psaout |> 
#   select(ltbi.soc, ltbi.int, tpt.soc, tpt.int, ltbinotpt.soc, ltbinotpt.int, noltbi.soc, noltbi.int,  tb, att.soc, att.int, noatt.soc, noatt.int) |> 
#   group_by(tb) |>
#   summarise(across(everything(), mean))
# 
# psaout2 |> 
#   # filter(id == 1) |>
#   select(tb, ltbi.soc, ltbi.int, tpt.soc, tpt.int, att.soc, att.int) |> 
#   group_by(tb) |>
#   summarise(across(everything(), mean))
# 
# # value not applied
# dc |> 
#   filter(id == 1) |>
#   select(tb, ltbi.soc, ltbi.int, tpt.soc, tpt.int, att.soc, att.int)
# 
# psaout2 |> 
#   filter(id == 1) |>
#   select(tb, ltbi.soc, ltbi.int, tpt.soc, tpt.int, att.soc, att.int)
# 

psaout2 |> 
  # filter(id == 1) |>
  select(tb, screen.soc, screen.int, ltbi.soc, ltbi.int, tpt.soc, tpt.int, att.soc, att.int) |> 
  group_by(tb) |>
  summarise(across(everything(), mean))

