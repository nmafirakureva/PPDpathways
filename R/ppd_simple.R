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
  sacases <- c('', 'AllattendNHS','noltfu', 'FUVisitsCost', 'DOTsCost', 
               'XrayCost', 'PrisonEscort', 'ContactTracing', 'InpatientCost')
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
agelevels <- c('15-100')
isoz <- c('GBR') #relevant countries

## prior parameters
PD <- read.csv(here('indata/ProbParms.csv')) # read in probability parameters
PS <- read.csv(here('indata/ProbParmsFixed.csv')) #read in fixed parameters
RD <- fread(gh('indata/RUParms.csv'))    # read resource use data
CD <- fread(gh('indata/CostParms.csv'))    # read cost data

# pre-process parameters
# names(PD)
# names(PS)
# names(RD)
# names(CD)

# unique(PD$ParameterName)
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

# unique(PD.INT$NAME)
PD0 <- rbind(PD.SOC, PD.INT)

# Fixed parameters to wide format
PS.SOC <- PS |> 
  mutate(ParameterName = paste0('soc.', ParameterName))

PS.INT <- PS |> 
  mutate(ParameterName = paste0('int.', ParameterName))

# unique(PD.INT$NAME)
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

if (SA == "noltfu") {
  # Prison GP assessment
  D[, int.prop.prison.gp.assessment := 1]
  # Attending NHS referral
  D[, int.prop.xray := 1]
  D[, int.prop.attend.nhs.referral := 1]
  # starting & completing ATT
  D[, int.prop.starting.att := 1]
  D[,int.prop.completing.att:=1]
  # # starting & completing TPT
  D[, int.prop.starting.tpt := 1]
  D[, int.prop.completing.tpt := 1]
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

# Extract cascades of care

# just do these for SA=='' or SA=='noltfu'
if (SA %in% c('', 'noltfu')) {
  # names(SOC.att.F)
  
  # ATT
  
  # soc
  D[["soc_att_screen"]] <- SOC.att.F$screenfun(D)
  D[["soc_att_presumtb"]] <- SOC.att.F$presumtbfun(D)
  D[["soc_att_prisongp"]] <- SOC.att.F$prisongpfun(D)
  D[["soc_att_tbsuspicion"]] <- SOC.att.F$tbsuspicionfun(D)
  D[["soc_att_attendnhs"]] <- SOC.att.F$attendnhsfun(D)
  D[["soc_att_tbdx"]] <- SOC.att.F$coprevtbfun(D)
  D[["soc_att_check"]] <- SOC.att.F$checkfun(D)
  D[["soc_att_completed"]] <- SOC.att.F$att_completefun(D)
  
  # int
  D[["int_att_screen"]] <- INT.att.F$screenfun(D)
  D[["int_att_presumtb"]] <- INT.att.F$presumtbfun(D)
  D[["int_att_prisongp"]] <- INT.att.F$prisongpfun(D)
  D[["int_att_tbsuspicion"]] <- INT.att.F$tbsuspicionfun(D)
  D[["int_att_attendnhs"]] <- INT.att.F$attendnhsfun(D)
  D[["int_att_tbdx"]] <- INT.att.F$coprevtbfun(D)
  D[["int_att_check"]] <- INT.att.F$checkfun(D)
  D[["int_att_completed"]] <- INT.att.F$att_completefun(D)
  
  ATTC <- D[, .(
    id, tb,
    soc_att_screen, soc_att_presumtb, soc_att_prisongp, soc_att_tbsuspicion, soc_att_attendnhs, soc_att_tbdx, soc_att_check, soc_att_completed,
    int_att_screen, int_att_presumtb, int_att_prisongp, int_att_tbsuspicion, int_att_attendnhs, int_att_tbdx, int_att_check, int_att_completed
  )]
  
  ## summary
  ATTCS <- ATTC[, lapply(.SD, mean), .SDcols = names(ATTC)[-c(1, 2)], by = tb]
  ATTCS <- melt(ATTCS, id = "tb")
  ATTCS[, c("arm", "outcome", "quantity") := tstrsplit(variable, split = "_")]
  ATTCS[, quantity := factor(quantity, levels = c("screen", "presumtb", "prisongp", "tbsuspicion", "attendnhs", "tbdx", "check", "completed"))]
  options(scipen = 999)
  (ATTCS <- dcast(data = ATTCS, formula = arm + quantity + outcome ~ tb, value.var = "value"))
  fn1 <- glue(here("outdata/ATTCS")) + SA + ".csv"
  fwrite(ATTCS, file = fn1)
  
  # names(SOC.tpt.F)
  
  # TPT
  
  # soc
  D[["soc_tpt_screen"]] <- SOC.tpt.F$screenfun(D)
  D[["soc_tpt_igra"]] <- SOC.tpt.F$igrafun(D)
  D[["soc_tpt_ltbi"]] <- SOC.tpt.F$ltbifun(D)
  D[["soc_tpt_stayingo3m"]] <- SOC.tpt.F$stayingo3mfun(D)
  D[["soc_tpt_attendnhs"]] <- SOC.tpt.F$attendnhsfun(D)
  D[["soc_tpt_check"]] <- SOC.tpt.F$checkfun(D)
  D[["soc_tpt_completed"]] <- SOC.tpt.F$tpt_completefun(D)
  
  # int
  D[["int_tpt_screen"]] <- INT.tpt.F$screenfun(D)
  D[["int_tpt_igra"]] <- INT.tpt.F$igrafun(D)
  D[["int_tpt_ltbi"]] <- INT.tpt.F$ltbifun(D)
  D[["int_tpt_stayingo3m"]] <- INT.tpt.F$stayingo3mfun(D)
  D[["int_tpt_attendnhs"]] <- INT.tpt.F$attendnhsfun(D)
  D[["int_tpt_check"]] <- INT.tpt.F$checkfun(D)
  D[["int_tpt_completed"]] <- INT.tpt.F$tpt_completefun(D)
  
  TPTC <- D[, .(
    id, tb,
    soc_tpt_screen, soc_tpt_igra, soc_tpt_ltbi, soc_tpt_stayingo3m, soc_tpt_attendnhs, soc_tpt_check, soc_tpt_completed,
    int_tpt_screen, int_tpt_igra, int_tpt_ltbi, int_tpt_stayingo3m, int_tpt_attendnhs, int_tpt_check , int_tpt_completed
  )]
  
  ## summary
  TPTCS <- TPTC[, lapply(.SD, mean), .SDcols = names(TPTC)[-c(1, 2)], by = tb]
  TPTCS <- melt(TPTCS, id = "tb")
  TPTCS[, c("arm", "outcome", "quantity") := tstrsplit(variable, split = "_")]
  TPTCS[, quantity := factor(quantity, levels = c("screen", "igra", "ltbi", "stayingo3m", "attendnhs", "check", "completed"))]
  options(scipen = 999)
  (TPTCS <- dcast(data = TPTCS, formula = arm + quantity + outcome ~ tb, value.var = "value"))
  fn1 <- glue(here("outdata/TPTCS")) + SA + ".csv"
  fwrite(TPTCS, file = fn1)
  
  unique(ATTCS$quantity)
  # plot these data
  att_levels <- c('screen','presumtb', 'prisongp', 'tbsuspicion', 'attendnhs', 'tbdx', 'check', 'completed')
  att_labels <- c('TB symptom screening','Presumtive TB', 'Prison GP assessment', 'Clinical TB suspicion',
                  'Attend NHS referral','TB diagnosis','ATT initiation','ATT completion')
  
  tmp <- ATTCS |>
    pivot_longer(cols = -c(arm, outcome, quantity), names_to = "tb", values_to = "value") |>
    mutate(arm = factor(arm, levels = c('soc','int'), labels = c('Standard of care', 'Intervention')),
           quantity = factor(quantity, levels = att_levels, labels = att_labels)) |>
    arrange(arm, quantity) %>%
    group_by(arm, tb) %>%
    mutate(
      value = cummin(value))
  
  setorder(tmp, arm, quantity)
  att_cascade <- tmp %>%
    arrange(quantity) %>%
    group_by(arm, tb) %>%
    mutate(
      # Create transition labels showing the proportion applied at each stage
      Transition_Label = scales::percent(value / lag(value, default = first(value)), accuracy = 1),
      x_start = as.numeric(quantity) - 0.67,  # Start of the arrow (midpoint between stages)
      x_end = as.numeric(quantity) - 0.35,    # End of the arrow (midpoint to next stage)
      Midpoint = (value + lag(value, default = first(value))) / 2  # Position for the label
    )|>
    mutate(Transition_Label = ifelse(grepl('screening', quantity), NA, Transition_Label),
           x_start = ifelse(grepl('screening', quantity), NA, x_start),
           x_end = ifelse(grepl('screening', quantity), NA, x_end))
  
  
  GP <- ggplot(tmp |> filter(tb=='TBD'), aes(x = quantity, y = value, fill = quantity)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)), vjust = -0.5, size = 3) +
    facet_wrap(arm ~ tb) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "ATT cascade ",
         x = "Cascade steps",
         y = "Proportion") +
    # theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "none") +
    scale_x_discrete(labels = ~ str_wrap(as.character(.x), 15))
  
  GP
  fn1 <- glue(here('plots/ATT_cascadeTBD')) + SA + '.png'
  ggsave(GP,file=fn1,w=6,h=5);
  
  # Create the Cascade Plot
  ATT <- ggplot(att_cascade |>  filter(tb=='TBD'), aes(x = quantity, y = value, fill = quantity)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)), vjust = -0.5, size = 2.5) +
    facet_wrap( ~ arm) +
    # Add arrows to indicate transition between stages
    geom_segment(aes(
      x = x_start,
      xend = x_end,
      y = ifelse(arm == "soc", Midpoint*0.65, Midpoint*0.77),
      yend = ifelse(arm == "soc", Midpoint*0.65, Midpoint*0.77)),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      color = "red") +
    # Add descriptive labels showing the proportion applied
    geom_text(aes(
      x = as.numeric(quantity) - 0.5,
      y = Midpoint*0.75,
      label = Transition_Label),
      color = "blue",
      vjust = 1.2,
      size = 2.5) +
    scale_fill_brewer(palette = "Set3") +
    labs(
      title = "ATT",
      x = "",
      y = "Proportion"
    ) +
    # theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    ylim(0, 1) +
    scale_x_discrete(labels = ~ str_wrap(as.character(.x), 8))
  
  ATT
  fn1 <- glue(here('plots/ATT_cascadeTBD_labels')) + SA + '.png'
  ggsave(ATT,file=fn1,w=12,h=5);
  
  unique(TPTCS$quantity)
  
  # plot these data
  tpt_levels <- c('screen','igra', 'ltbi', 'stayingo3m', 'check', 'completed')
  tpt_labels <- c('TB symptom screening','IGRA tested', 'IGRA test positive',
                  'Staying over 3 months', 'TPT initiation', 'TPT completion')
  
  tmp <- TPTCS |>
    pivot_longer(cols = -c(arm, outcome, quantity), names_to = "tb", values_to = "value") |>
    filter(quantity %in% tpt_levels) |>
    mutate(arm = factor(arm, levels = c('soc','int'), labels = c('Standard of care', 'Intervention')),
           quantity = factor(quantity, levels = tpt_levels, labels = tpt_labels)) |>
    arrange(arm, quantity) %>%
    group_by(arm, tb) %>%
    mutate(
      value = cummin(value))
  
  tpt_cascade <- tmp %>%
    arrange(quantity) %>%
    group_by(arm, tb) %>%
    mutate(
      # Create transition labels showing the proportion applied at each stage
      Transition_Label = scales::percent(value / lag(value, default = first(value)), accuracy = 1),
      x_start = as.numeric(quantity) - 0.67,  # Start of the arrow (midpoint between stages)
      x_end = as.numeric(quantity) - 0.35,    # End of the arrow (midpoint to next stage)
      Midpoint = (value + lag(value, default = first(value))) / 2  # Position for the label
    )|>
    mutate(Transition_Label = ifelse(grepl('screening', quantity), NA, Transition_Label),
           x_start = ifelse(grepl('screening', quantity), NA, x_start),
           x_end = ifelse(grepl('screening', quantity), NA, x_end))
  
  GP <- ggplot(tmp |> filter(tb=='TBI'), aes(x = quantity, y = value, fill = quantity)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)), vjust = -0.5, size = 3) +
    facet_wrap(arm ~ tb) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "TPT cascade ",
         x = "Cascade steps",
         y = "Proportion") +
    # theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "none") +
    scale_x_discrete(labels = ~ str_wrap(as.character(.x), 15))
  
  GP
  fn1 <- glue(here('plots/TPT_cascadeTBI')) + SA + '.png'
  ggsave(GP,file=fn1,w=6,h=5);
  
  # Create the Cascade Plot
  TPT <- ggplot(tpt_cascade |> filter(tb=='TBI'), aes(x = quantity, y = value, fill = quantity)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)), vjust = -0.5, size = 2.5) +
    facet_wrap( ~ arm) +
    # Add arrows to indicate transition between stages
    geom_segment(aes(
      x = x_start,
      xend = x_end,
      y = ifelse(arm == "soc", Midpoint*0.3, Midpoint*0.77),
      yend = ifelse(arm == "soc", Midpoint*0.3, Midpoint*0.77)),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      color = "red") +
    # Add descriptive labels showing the proportion applied
    geom_text(aes(
      x = as.numeric(quantity) - 0.5,
      y = ifelse(arm == "soc" & Midpoint>0.1, Midpoint*0.5, Midpoint*0.75),
      label = Transition_Label),
      color = "blue",
      vjust = 1.2,
      size = 2.5) +
    scale_fill_brewer(palette = "Set3") +
    labs(
      title = "TPT",
      x = "",
      y = "Proportion"
    ) +
    # theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    ylim(0, 1) +
    scale_x_discrete(labels = ~ str_wrap(as.character(.x), 8))
  
  TPT
  fn1 <- glue(here('plots/TPT_cascadeTBI_labels')) + SA + '.png'
  ggsave(TPT,file=fn1,w=6,h=5);
  names(SOC.tpt.F)
  
  library(patchwork)
  combined <- ATT / TPT
  combined
  fn1 <- glue(here('plots/ATT_TPT_cascades')) + SA + '.png'
  ggsave(combined,file=fn1,w=11,h=6.5);
  
  fn1 <- glue(here('plots/ATT_TPT_cascades')) + SA + '.eps'
  ggsave(combined,file=fn1,w=11,h=6.5);
  
  # pdf 
  fn1 <- glue(here('plots/ATT_TPT_cascades')) + SA + '.pdf'
  ggsave(combined,file=fn1,w=11,h=6.5);
  
  # Extract costs
  # names(SOC.att.F)[grep("cost.", names(SOC.att.F))]
  
  # ATT cascade
  D[["soc_att_cost.screen"]] <- SOC.att.F$cost.screenfun(D)
  D[["soc_att_cost.prison.gp.assessment"]] <- SOC.att.F$cost.prison.gp.assessmentfun(D)
  D[["soc_att_cost.prison.isolation.assess"]] <- SOC.att.F$cost.prison.isolation.assessfun(D)
  D[["soc_att_cost.chest.xray"]] <- SOC.att.F$cost.chest.xrayfun(D)
  D[["soc_att_cost.nhs.tb.service.assess"]] <- SOC.att.F$cost.nhs.tb.service.assessfun(D)
  D[["soc_att_xray"]] <- SOC.att.F$xrayfun(D)
  D[["soc_att_cost.prison.escort.xray"]] <- SOC.att.F$cost.prison.escort.xrayfun(D)
  D[["soc_att_cost.prison.escort.assess"]] <- SOC.att.F$cost.prison.escort.assessfun(D)
  D[["soc_att_cost.tb.investigations"]] <- SOC.att.F$cost.tb.investigationsfun(D)
  D[["soc_att_cost.igra.test"]] <- SOC.att.F$cost.igra.testfun(D)
  
  D[["soc_att_tbdx"]] <- SOC.att.F$coprevtbfun(D)
  D[["soc_att_isolation"]] <- SOC.att.F$isolatedfun(D)
  D[["soc_att_att"]] <- SOC.att.F$attfun(D)
  D[["soc_att_cost.prison.isolation"]] <- SOC.att.F$cost.prison.isolation.attfun(D)
  D[["soc_att_cost.nhs.tb.service"]] <- SOC.att.F$cost.nhs.tb.servicefun(D)
  D[["soc_att_cost.prison.escort"]] <- SOC.att.F$cost.prison.escortfun(D)
  D[["soc_att_cost.drugs"]] <- SOC.att.F$cost.drugsfun(D)
  D[["soc_att_cost.dots"]] <- SOC.att.F$cost.dotsfun(D)
  D[["soc_att_cost.nhs.ipd"]] <- SOC.att.F$cost.nhs.ipdfun(D)
  D[["soc_att_cost.prison.bedwatch"]] <- SOC.att.F$cost.prison.bedwatchfun(D)
  D[["soc_att_cost.contact.management"]] <- SOC.att.F$cost.contact.managementfun(D)
  D[["soc_att_cost"]] <- SOC.att.F$costfun(D)
  
  # int
  D[["int_att_cost.screen"]] <- INT.att.F$cost.screenfun(D)
  D[["int_att_cost.prison.gp.assessment"]] <- INT.att.F$cost.prison.gp.assessmentfun(D)
  D[["int_att_cost.prison.isolation.assess"]] <- INT.att.F$cost.prison.isolation.assessfun(D)
  D[["int_att_xray"]] <- INT.att.F$xrayfun(D)
  D[["int_att_cost.chest.xray"]] <- INT.att.F$cost.chest.xrayfun(D)
  D[["int_att_cost.nhs.tb.service.assess"]] <- INT.att.F$cost.nhs.tb.service.assessfun(D)
  D[["int_att_cost.prison.escort.xray"]] <- INT.att.F$cost.prison.escort.xrayfun(D)
  D[["int_att_cost.prison.escort.assess"]] <- INT.att.F$cost.prison.escort.assessfun(D)
  D[["int_att_cost.tb.investigations"]] <- INT.att.F$cost.tb.investigationsfun(D)
  D[["int_att_cost.igra.test"]] <- INT.att.F$cost.igra.testfun(D)
  
  D[["int_att_tbdx"]] <- INT.att.F$coprevtbfun(D)
  D[["int_att_isolation"]] <- INT.att.F$isolatedfun(D)
  D[["int_att_att"]] <- INT.att.F$attfun(D)
  D[["int_att_cost.prison.isolation"]] <- INT.att.F$cost.prison.isolation.attfun(D)
  D[["int_att_cost.nhs.tb.service"]] <- INT.att.F$cost.nhs.tb.servicefun(D)
  D[["int_att_cost.prison.escort"]] <- INT.att.F$cost.prison.escortfun(D)
  D[["int_att_cost.drugs"]] <- INT.att.F$cost.drugsfun(D)
  D[["int_att_cost.dots"]] <- INT.att.F$cost.dotsfun(D)
  D[["int_att_cost.nhs.ipd"]] <- INT.att.F$cost.nhs.ipdfun(D)
  D[["int_att_cost.prison.bedwatch"]] <- INT.att.F$cost.prison.bedwatchfun(D)
  D[["int_att_cost.contact.management"]] <- INT.att.F$cost.contact.managementfun(D)
  D[["int_att_cost"]] <- INT.att.F$costfun(D)
  D[["soc_att_igra"]] <- SOC.att.F$igrafun(D)
  D[["int_att_igra"]] <- INT.att.F$igrafun(D)
  
  # unique(ATTCS$quantity)
  # names(SOC.att.F)[!grepl('cost', names(SOC.att.F))]
  
  ## condition costs on outcome
  
  ATTCST <- D[, .(
    id, tb,
    soc_att_cost.screen=soc_att_cost.screen/soc_att_check,
    soc_att_cost.prison.gp.assessment=soc_att_cost.prison.gp.assessment/soc_att_check,
    soc_att_cost.prison.isolation.assessment=soc_att_cost.prison.isolation.assess/soc_att_check,
    soc_att_cost.chest.xray=soc_att_cost.chest.xray/soc_att_check,
    soc_att_cost.prison.escort.xray = soc_att_cost.prison.escort.xray/soc_att_check,
    soc_att_cost.nhs.tb.service.assessment = soc_att_cost.nhs.tb.service.assess/soc_att_check,
    soc_att_cost.prison.escort.assessment = soc_att_cost.prison.escort.assess/soc_att_check,
    soc_att_cost.tb.investigations=soc_att_cost.tb.investigations/soc_att_check,
    soc_att_cost.igra.test=soc_att_cost.igra.test/soc_att_check,
    soc_att_cost.prison.isolation = soc_att_cost.prison.isolation/soc_att_check,
    soc_att_cost.nhs.tb.service = soc_att_cost.nhs.tb.service/soc_att_check,
    soc_att_cost.prison.escort = soc_att_cost.prison.escort/soc_att_check,
    soc_att_cost.drugs = soc_att_cost.drugs/soc_att_check,
    soc_att_cost.dots = soc_att_cost.dots/soc_att_check,
    soc_att_cost.nhs.ipd = soc_att_cost.nhs.ipd/(soc_att_check),
    soc_att_cost.prison.bedwatch = soc_att_cost.prison.bedwatch/(soc_att_check),
    soc_att_cost.contact.management=soc_att_cost.contact.management/soc_att_check,
    soc_att_cost = soc_att_cost/soc_att_check,
    int_att_cost.screen = int_att_cost.screen/int_att_check,
    int_att_cost.prison.gp.assessment = int_att_cost.prison.gp.assessment/int_att_check,
    int_att_cost.prison.isolation.assessment = int_att_cost.prison.isolation.assess/int_att_check,
    int_att_cost.chest.xray = int_att_cost.chest.xray/int_att_check,
    int_att_cost.prison.escort.xray = int_att_cost.prison.escort.xray/int_att_check,
    int_att_cost.nhs.tb.service.assessment = int_att_cost.nhs.tb.service.assess/int_att_check,
    int_att_cost.prison.escort.assessment = int_att_cost.prison.escort.assess/int_att_check,
    int_att_cost.tb.investigations = int_att_cost.tb.investigations/int_att_check,
    int_att_cost.igra.test = int_att_cost.igra.test/int_att_check,
    int_att_cost.prison.isolation = int_att_cost.prison.isolation/int_att_check,
    int_att_cost.nhs.tb.service = int_att_cost.nhs.tb.service/int_att_check,
    int_att_cost.prison.escort = int_att_cost.prison.escort/int_att_check,
    int_att_cost.drugs = int_att_cost.drugs/int_att_check,
    int_att_cost.dots = int_att_cost.dots/int_att_check,
    int_att_cost.nhs.ipd = int_att_cost.nhs.ipd/(int_att_check),
    int_att_cost.prison.bedwatch = int_att_cost.prison.bedwatch/(int_att_check),
    int_att_cost.contact.management = int_att_cost.contact.management/int_att_check,
    int_att_cost = int_att_cost/int_att_check
  )]
  
  ## summary
  ATTCSTS <- ATTCST[, lapply(.SD, mean), .SDcols = names(ATTCST)[-c(1, 2)], by = tb]
  ATTCSTS <- melt(ATTCSTS, id = "tb")
  ATTCSTS[, c("arm", "outcome", "quantity") := tstrsplit(variable, split = "_")]
  
  # unique(ATTCSTS$quantity)
  ATTCSTS[, quantity := factor(quantity,
                               levels = c("cost.screen",
                                          "cost.prison.gp.assessment", "cost.chest.xray","cost.prison.isolation.assessment",
                                          "cost.prison.escort.xray",
                                          "cost.nhs.tb.service.assessment",	"cost.prison.escort.assessment","cost.tb.investigations",
                                          "cost.igra.test",	"cost.prison.isolation","cost.nhs.tb.service",	"cost.prison.escort",
                                          "cost.drugs",	"cost.dots","cost.nhs.ipd",	"cost.prison.bedwatch",
                                          "cost.contact.management","cost"))]
  options(scipen = 999)
  (ATTCSTS <- dcast(data = ATTCSTS, formula = arm + quantity + outcome ~ tb, value.var = "value"))
  fn1 <- glue(here("outdata/ATTCSTS")) + SA + ".csv"
  fwrite(ATTCSTS, file = fn1)
  
  # ATTCSTS |> filter(arm=='int')
  ATTCSTS <- ATTCSTS %>%
    mutate(across(all_of(names(ATTCSTS)[-c(1:3)]), ~ ifelse(is.nan(.), 0, .)))
  
  # # check if sum of costs is equal to total cost
  # ATTCSTS[quantity != "cost", lapply(.SD, sum), .SDcols = names(ATTCSTS)[-c(1:3)], by = arm]
  # 
  # ATTCSTS |>
  #   filter(quantity == "cost") |> select(arm,TBD:noTB)
  # ATTCSTS |>
  #   filter(quantity != "cost") |>
  #   group_by(arm) |>
  #   summarise(across(TBD:noTB, \(x) sum(x, na.rm = TRUE), .names = "sum_{.col}"))
  
  x <- ATTCS |>
    filter(!quantity %in% c("completed", "attendnhs","tbdx")) |>
    pivot_longer(cols = c(TBD:noTB), names_to = "tb", values_to = "value") |>
    group_by(arm,tb) |>
    mutate(last = value[quantity=='check'],
           prop = value/last) |>
    select(-last)
  
  # x |> filter(arm=='int' & tb=='TBD')
  # x |> filter(arm=='int' & tb=='TBI')
  # x |> filter(arm=='int' & tb=='noTB')
  
  tmp <- ATTCSTS |>
    mutate(qty = case_when(
      quantity == "cost.screen" ~ "TB evaluation",
      quantity %in% c('cost.prison.gp.assessment','cost.chest.xray','cost.prison.escort.xray','cost.prison.isolation.assessment',
                      'cost.nhs.tb.service.assessment','cost.prison.escort.assessment','cost.tb.investigations') ~ "TB evaluation",
      quantity %in% c('cost.igra.test','cost.prison.isolation','cost.nhs.tb.service','cost.prison.escort',
                      'cost.drugs','cost.nhs.visits','cost.dots','cost.nhs.ipd','cost.prison.bedwatch',
                      'cost.contact.management') ~ "ATT",
      TRUE ~ 'total'
    ),
    quantity = gsub(".assessment", "", quantity),
    quantity = ifelse(quantity=='cost.prison.gp', 'cost.prison.gp.assessment', 
                      ifelse(quantity=='cost.screen', 'cost.screening', quantity)),
    qty = factor(qty, levels = c('TB evaluation', 'ATT'), ordered = TRUE)) |>
    filter(!quantity %in% c('cost','cost.att')) |>
    mutate(arm = factor(arm, levels = c('soc','int'), labels = c('Standard of care', 'Intervention')),
           quantity = ifelse(quantity=='cost.prison.escort.xray', 'cost.prison.escort', quantity),
           quantity = factor(quantity,
                             levels = c("cost.screening",
                                        "cost.prison.gp.assessment", "cost.chest.xray","cost.prison.isolation.assessment",
                                        "cost.nhs.tb.service.assessment",	"cost.prison.escort.assessment","cost.tb.investigations",
                                        "cost.igra.test",	"cost.prison.isolation","cost.nhs.tb.service",	"cost.prison.escort",
                                        "cost.drugs",	"cost.dots","cost.nhs.ipd",	"cost.prison.bedwatch",
                                        "cost.contact.management"), ordered = TRUE)) |>
    pivot_longer(cols = -c(arm, outcome, quantity, qty), names_to = "tb", values_to = "value") |>
    mutate(tb = factor(tb, levels = c('TBD', 'TBI','noTB'), ordered = TRUE)) |>
    group_by(arm, outcome, quantity, qty,tb) |>
    summarise(value = sum(value, na.rm = TRUE))
  
  lvls <- gsub('cost.', '', unique(tmp$quantity))

  
  
  # library(RColorBrewer)
  # n <- 30
  # colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
  # col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
  # col <- sample(col_vec, n)
  # area <- rep(1,n)
  # pie(area, col = col)
  
  col <- c("#CCECE6","#E7E1EF","#B8E186","#ECE2F0","#CCEBC5","#EF3B2C",
           "#542788","#6BAED6","#7BCCC4","#D4B9DA","#EDF8B1","#33A02C",  "#FDBF6F")
  GP <- tmp |>
    filter(tb == 'TBD') |>
    mutate(quantity = gsub('cost.', '', quantity),
           quantity = factor(quantity,levels = lvls)) |>
    arrange(arm, tb, qty, value) |>
    ggplot(aes(x = qty, y = value, fill = reorder(quantity, value))) +
    geom_bar(stat="identity", position="fill", width = 0.5) +
    facet_grid(.~arm) +
    scale_fill_manual(values=col[1:13],
                      breaks = lvls) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = 'top', 
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    labs(title = "", x = "", y = "Proportion")
  
  GP 
  fn1 <- glue(here('plots/ATT_cascade_costs')) + SA + '.png'
  ggsave(file=fn1,, width = 10, height = 5)
  
  # TPT cascade costs
  
  # TPT cascade
  
  # names(SOC.tpt.F)[grep("cost.", names(SOC.tpt.F))]
  
  D[["soc_tpt_cost.screen"]] <- SOC.tpt.F$cost.screenfun(D)
  D[["soc_tpt_cost.prison.gp.assessment"]] <- SOC.tpt.F$cost.prison.gp.assessmentfun(D)
  D[["soc_tpt_cost.prison.isolation.assess"]] <- SOC.tpt.F$cost.prison.isolation.assessfun(D)
  D[["soc_tpt_cost.chest.xray"]] <- SOC.tpt.F$cost.chest.xrayfun(D)
  D[["soc_tpt_cost.nhs.tb.service.assess"]] <- SOC.tpt.F$cost.nhs.tb.service.assessfun(D)
  D[["soc_tpt_xray"]] <- SOC.tpt.F$xrayfun(D)
  D[["soc_tpt_cost.prison.escort.xray"]] <- SOC.tpt.F$cost.prison.escort.xrayfun(D)
  D[["soc_tpt_cost.prison.escort.assess"]] <- SOC.tpt.F$cost.prison.escort.assessfun(D)
  D[["soc_tpt_cost.tb.investigations"]] <- SOC.tpt.F$cost.tb.investigationsfun(D)
  D[["soc_tpt_cost.igra.test"]] <- SOC.tpt.F$cost.igra.testfun(D)
  
  D[["soc_tpt_tbdx"]] <- SOC.tpt.F$coprevtbfun(D)
  D[["soc_tpt_isolation"]] <- SOC.tpt.F$isolatedfun(D)
  D[["soc_tpt_att"]] <- SOC.tpt.F$attfun(D)
  D[["soc_tpt_cost.prison.isolation"]] <- SOC.tpt.F$cost.prison.isolation.attfun(D)
  D[["soc_tpt_cost.nhs.tb.service"]] <- SOC.tpt.F$cost.nhs.tb.servicefun(D)
  D[["soc_tpt_cost.prison.escort"]] <- SOC.tpt.F$cost.prison.escortfun(D)
  D[["soc_tpt_cost.drugs"]] <- SOC.tpt.F$cost.drugsfun(D)
  D[["soc_tpt_cost.dots"]] <- SOC.tpt.F$cost.dotsfun(D)
  D[["soc_tpt_cost.nhs.ipd"]] <- SOC.tpt.F$cost.nhs.ipdfun(D)
  D[["soc_tpt_cost.prison.bedwatch"]] <- SOC.tpt.F$cost.prison.bedwatchfun(D)
  D[["soc_tpt_cost.contact.management"]] <- SOC.tpt.F$cost.contact.managementfun(D)
  D[["soc_tpt_cost"]] <- SOC.tpt.F$costfun(D)
  
  # int
  D[["int_tpt_cost.screen"]] <- INT.tpt.F$cost.screenfun(D)
  D[["int_tpt_cost.prison.gp.assessment"]] <- INT.tpt.F$cost.prison.gp.assessmentfun(D)
  D[["int_tpt_cost.prison.isolation.assess"]] <- INT.tpt.F$cost.prison.isolation.assessfun(D)
  D[["int_tpt_xray"]] <- INT.tpt.F$xrayfun(D)
  D[["int_tpt_cost.chest.xray"]] <- INT.tpt.F$cost.chest.xrayfun(D)
  D[["int_tpt_cost.nhs.tb.service.assess"]] <- INT.tpt.F$cost.nhs.tb.service.assessfun(D)
  D[["int_tpt_cost.prison.escort.xray"]] <- INT.tpt.F$cost.prison.escort.xrayfun(D)
  D[["int_tpt_cost.prison.escort.assess"]] <- INT.tpt.F$cost.prison.escort.assessfun(D)
  D[["int_tpt_cost.tb.investigations"]] <- INT.tpt.F$cost.tb.investigationsfun(D)
  D[["int_tpt_cost.igra.test"]] <- INT.tpt.F$cost.igra.testfun(D)
  
  D[["int_tpt_tbdx"]] <- INT.tpt.F$coprevtbfun(D)
  D[["int_tpt_isolation"]] <- INT.tpt.F$isolatedfun(D)
  D[["int_tpt_att"]] <- INT.tpt.F$attfun(D)
  D[["int_tpt_cost.prison.isolation"]] <- INT.tpt.F$cost.prison.isolation.attfun(D)
  D[["int_tpt_cost.nhs.tb.service"]] <- INT.tpt.F$cost.nhs.tb.servicefun(D)
  D[["int_tpt_cost.prison.escort"]] <- INT.tpt.F$cost.prison.escortfun(D)
  D[["int_tpt_cost.drugs"]] <- INT.tpt.F$cost.drugsfun(D)
  D[["int_tpt_cost.dots"]] <- INT.tpt.F$cost.dotsfun(D)
  D[["int_tpt_cost.nhs.ipd"]] <- INT.tpt.F$cost.nhs.ipdfun(D)
  D[["int_tpt_cost.prison.bedwatch"]] <- INT.tpt.F$cost.prison.bedwatchfun(D)
  D[["int_tpt_cost.contact.management"]] <- INT.tpt.F$cost.contact.managementfun(D)
  D[["int_tpt_cost"]] <- INT.tpt.F$costfun(D)
  
  D[["soc_tpt_igra"]] <- SOC.tpt.F$igrafun(D)
  D[["int_tpt_igra"]] <- INT.tpt.F$igrafun(D)
  
  # unique(ATTCS$quantity)
  # names(SOC.tpt.F)[!grepl('cost', names(SOC.tpt.F))]
  
  ## condition costs on outcome
  
  TPTCST <- D[, .(
    id, tb,
    soc_tpt_cost.screen=soc_tpt_cost.screen/soc_tpt_check,
    soc_tpt_cost.prison.gp.assessment=soc_tpt_cost.prison.gp.assessment/soc_tpt_check,
    soc_tpt_cost.prison.isolation.assessment=soc_tpt_cost.prison.isolation.assess/soc_tpt_check,
    soc_tpt_cost.chest.xray=soc_tpt_cost.chest.xray/soc_tpt_check,
    soc_tpt_cost.prison.escort.xray = soc_tpt_cost.prison.escort.xray/soc_tpt_check,
    soc_tpt_cost.nhs.tb.service.assessment = soc_tpt_cost.nhs.tb.service.assess/soc_tpt_check,
    soc_tpt_cost.prison.escort.assessment = soc_tpt_cost.prison.escort.assess/soc_tpt_check,
    soc_tpt_cost.tb.investigations=soc_tpt_cost.tb.investigations/soc_tpt_check,
    soc_tpt_cost.igra.test=soc_tpt_cost.igra.test/soc_tpt_check,
    soc_tpt_cost.prison.isolation = soc_tpt_cost.prison.isolation/soc_tpt_check,
    soc_tpt_cost.nhs.tb.service = soc_tpt_cost.nhs.tb.service/soc_tpt_check,
    soc_tpt_cost.prison.escort = soc_tpt_cost.prison.escort/soc_tpt_check,
    soc_tpt_cost.drugs = soc_tpt_cost.drugs/soc_tpt_check,
    soc_tpt_cost.dots = soc_tpt_cost.dots/soc_tpt_check,
    soc_tpt_cost.nhs.ipd = soc_tpt_cost.nhs.ipd/(soc_tpt_check),
    soc_tpt_cost.prison.bedwatch = soc_tpt_cost.prison.bedwatch/(soc_tpt_check),
    soc_tpt_cost.contact.management=soc_tpt_cost.contact.management/soc_tpt_check,
    soc_tpt_cost = soc_tpt_cost/soc_tpt_check,
    int_tpt_cost.screen = int_tpt_cost.screen/int_tpt_check,
    int_tpt_cost.prison.gp.assessment = int_tpt_cost.prison.gp.assessment/int_tpt_check,
    int_tpt_cost.prison.isolation.assessment = int_tpt_cost.prison.isolation.assess/int_tpt_check,
    int_tpt_cost.chest.xray = int_tpt_cost.chest.xray/int_tpt_check,
    int_tpt_cost.prison.escort.xray = int_tpt_cost.prison.escort.xray/int_tpt_check,
    int_tpt_cost.nhs.tb.service.assessment = int_tpt_cost.nhs.tb.service.assess/int_tpt_check,
    int_tpt_cost.prison.escort.assessment = int_tpt_cost.prison.escort.assess/int_tpt_check,
    int_tpt_cost.tb.investigations = int_tpt_cost.tb.investigations/int_tpt_check,
    int_tpt_cost.igra.test = int_tpt_cost.igra.test/int_tpt_check,
    int_tpt_cost.prison.isolation = int_tpt_cost.prison.isolation/int_tpt_check,
    int_tpt_cost.nhs.tb.service = int_tpt_cost.nhs.tb.service/int_tpt_check,
    int_tpt_cost.prison.escort = int_tpt_cost.prison.escort/int_tpt_check,
    int_tpt_cost.drugs = int_tpt_cost.drugs/int_tpt_check,
    int_tpt_cost.dots = int_tpt_cost.dots/int_tpt_check,
    int_tpt_cost.nhs.ipd = int_tpt_cost.nhs.ipd/(int_tpt_check),
    int_tpt_cost.prison.bedwatch = int_tpt_cost.prison.bedwatch/(int_tpt_check),
    int_tpt_cost.contact.management = int_tpt_cost.contact.management/int_tpt_check,
    int_tpt_cost = int_tpt_cost/int_tpt_check
  )]
  
  ## summary
  TPTCSTS <- TPTCST[, lapply(.SD, mean), .SDcols = names(TPTCST)[-c(1, 2)], by = tb]
  TPTCSTS <- melt(TPTCSTS, id = "tb")
  TPTCSTS[, c("arm", "outcome", "quantity") := tstrsplit(variable, split = "_")]
  (TPTCSTS <- dcast(data = TPTCSTS, formula = arm + quantity + outcome ~ tb, value.var = "value"))
  TPTCSTS <- TPTCSTS[TBI!=0, -c("TBD", "noTB")]
  
  # unique(TPTCSTS$quantity)
  TPTCSTS[, quantity := factor(quantity,
                               levels = c("cost.screen",
                                          "cost.prison.gp.assessment", "cost.chest.xray","cost.prison.isolation.assessment",
                                          "cost.prison.escort.xray",
                                          "cost.nhs.tb.service.assessment",	"cost.prison.escort.assessment","cost.tb.investigations",
                                          "cost.igra.test","cost.nhs.tb.service",	"cost.prison.escort",
                                          "cost.drugs",	"cost.dots","cost"))]
  options(scipen = 999)
  fn1 <- glue(here("outdata/TPTCSTS")) + SA + ".csv"
  fwrite(TPTCSTS, file = fn1)
  
  # TPTCSTS |> filter(arm=='int')
  TPTCSTS <- TPTCSTS %>%
    mutate(across(all_of(names(TPTCSTS)[-c(1:3)]), ~ ifelse(is.nan(.), 0, .)))
  
  # # check if sum of costs is equal to total cost
  # TPTCSTS[quantity != "cost", lapply(.SD, sum), .SDcols = names(TPTCSTS)[-c(1:3)], by = arm]
  # 
  # TPTCSTS |>
  #   filter(quantity == "cost") |> select(arm,TBI)
  # TPTCSTS |>
  #   filter(quantity != "cost") |>
  #   group_by(arm) |>
  #   summarise(across(TBI, \(x) sum(x, na.rm = TRUE), .names = "sum_{.col}"))
  
  x <- TPTCS |>
    filter(!quantity %in% c("attendnhs")) |>
    pivot_longer(cols = c(TBD:noTB), names_to = "tb", values_to = "value") |>
    group_by(arm,tb) |>
    mutate(last = value[quantity=='completed'],
           prop = value/last) |>
    select(-last)
  
  # x |> filter(arm=='int' & tb=='TBI')
  # x |> filter(arm!='int' & tb=='TBI')
  
  # unique(TPTCSTS$quantity)
  tmp <- TPTCSTS |>
    mutate(qty = case_when(
      quantity == "cost.screen" ~ "TB evaluation",
      quantity %in% c('cost.prison.gp.assessment','cost.chest.xray','cost.prison.escort.xray','cost.prison.isolation.assessment',
                      'cost.nhs.tb.service.assessment','cost.prison.escort.assessment','cost.tb.investigations') ~ "TB evaluation",
      quantity %in% c('cost.igra.test','cost.nhs.tb.service','cost.prison.escort',
                      'cost.drugs','cost.nhs.visits','cost.dots') ~ "TPT",
      TRUE ~ 'total'
    ),
    quantity = gsub(".assessment", "", quantity),
    quantity = ifelse(quantity=='cost.prison.gp', 'cost.prison.gp.assessment', 
                      ifelse(quantity=='cost.screen', 'cost.screening', quantity)),
    qty = factor(qty, levels = c('screening', 'TB evaluation', 'TPT'), ordered = TRUE)) |>
    filter(!quantity %in% c('cost','cost.att')) |>
    mutate(arm = factor(arm, levels = c('soc','int'), labels = c('Standard of care', 'Intervention')),
           quantity = ifelse(quantity=='cost.prison.escort.xray', 'cost.prison.escort',
                             ifelse(quantity=='cost.prison.isolation.assessment', 'cost.prison.isolation',quantity)),
           quantity = factor(quantity,
                             levels = c("cost.screening",
                                        "cost.prison.gp.assessment", "cost.chest.xray","cost.prison.isolation",
                                        "cost.nhs.tb.service.assessment",	"cost.prison.escort.assessment","cost.tb.investigations",
                                        "cost.igra.test","cost.nhs.tb.service",	"cost.prison.escort",
                                        "cost.drugs",	"cost.dots"), ordered = TRUE)) |>
    pivot_longer(cols = -c(arm, outcome, quantity, qty), names_to = "tb", values_to = "value") |>
    mutate(tb = factor(tb, levels = c('TBD', 'TBI','noTB'), ordered = TRUE)) |>
    group_by(arm, outcome, quantity, qty,tb) |>
    summarise(value = sum(value, na.rm = TRUE))
  
  lvls <- gsub('cost.', '', unique(tmp$quantity))
  
  GP <- tmp |>
    mutate(quantity = gsub('cost.', '', quantity),
           quantity = factor(quantity,levels = lvls)) |>
    ggplot(aes(x = qty, y = value, fill = fct_reorder(quantity, value))) +
    geom_bar(position="fill", stat="identity", width = 0.5) +
    facet_grid(.~arm) +
    scale_fill_manual(values=col[1:10],
                      breaks = lvls) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = 'top', 
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    labs(title = "", x = "", y = "Proportion")
  GP 
  fn1 <- glue(here('plots/TPT_cascade_costs')) + SA + '.png'
  ggsave(file=fn1,, width = 10, height = 5)
  
}
