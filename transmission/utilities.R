## utility/helper functions
library(ecrins) #install from petedodd/ecrins
library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(glue)
library(scico)

source(system.file("extdata", "parameters.R", package="ecrins")) #this loads some fns for PSA
data(parms)                       #some default parameters

## introducing remand TBI prevalence separately
parms$tbip <- 0.1
hyperparms[["tbip"]] <- list(shape1=10,shape2=90)

## prisons parameters
## parms
load(here('churn/outdata/smpsd.Rdata')) #TODO change/update

## tree outputs
load(here('outdata/DR.Rdata'))
DRS <- fread(here('outdata/DRS.csv'))

DR.basecase <- copy(DR) # keep

## pick up the extra sensitivity analyses
saz <- dir(here("outdata"), pattern = ".Rdata")
saz <- saz[saz != "DR.Rdata"]
saznmz <- gsub("DR", "", gsub("\\.Rdata", "", saz))
Dlist <- list()
for(i in 1:length(saz)){
  load(here('outdata',saz[i]))
  Dlist[[saznmz[i]]] <- copy(DR)
}


## helpers
lo <- function(x) quantile(x, probs = 0.025)
hi <- function(x) quantile(x, probs = 1 - 0.025)
## rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
brkt <- function(M, L, H, ndp = 0) {
  paste0(
    round(M, ndp), " (",
    round(L, ndp), " to ",
    round(H, ndp), ")"
  )
}

odds <- function(x) x / (1 - x)
inv.odds <- function(x) x / (1 + x)
odds.ratio <- function(x, y) odds(x) / odds(y)


testfun <- function(times, p) {
  if(p<0) warning('p<0!')
  times
}


## run model but catch warnings
runmodelsafely <- function(times,p) {
   r <-
     tryCatch(
       withCallingHandlers(
         {
           error_val <- 0
           list(y = runmodel(times,p), error_val = error_val)
         },
         warning = function(e) {
           error_val <<- 1
           invokeRestart("muffleWarning")
         }
       ),
       error = function(e) {
         return(list(y = NA, error_val = 2))
       },
       finally = {
       }
     )
   return(r)
}


## functions for reformatting outputs
getpno <- function(x){
  a <- str_extract(x,"\\[(.*?),")
  a <- gsub('\\[','',a)
  a <- gsub(',','',a)
  as.integer(a)
}
gettno <- function(x){
  a <- str_extract(x,",(.*?)\\]")
  a <- gsub('\\]','',a)
  a <- gsub(',','',a)
  as.integer(a)
}
pcts <- c('remand',
          'short stay',
          'long stay',
          'open',
          'previously detained')
tcts <- c('never TPT','on TPT','previous TPT')
getppd <- function(x) pcts[getpno(x)]
gettpt <- function(x) tcts[gettno(x)]
## ## test
## getppd(tail(ydm$variable))
## gettpt(tail(ydm$variable))
## function to map output into dt
output2dt <- function(y){
  yd <- as.data.table(y)
  ydm <- melt(yd,id='t')
  ## transform
  ydm[,population:=getppd(variable)]
  ydm[,tpt:=gettpt(variable)]
  ydm[,variable:=gsub("\\[.+\\]","",variable)]
  ## recat
  ydm$population <- factor(ydm$population)
  ydm$tpt <- factor(ydm$tpt)
  ydm$variable <- factor(ydm$variable)
  ## ## PPD pops
  ## ydm[!is.na(population),unique(variable)] #all TB states
  ydm
}


## CEAC
make.ceac <- function(CEA,lamz){
  crv <- lamz
  for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
  crv
}
cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73",
               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
               "darkorchid1")



revise.flow.parms <- function(parms, # original parameter template
                              smpsd, # sample of flow parameters
                              j      # which row of sample to use
                              ){
  ## safety normalize proportions
  parmfracnmz <- c('parm_frac_SD', 'parm_frac_CD', 'parm_frac_E', 'parm_frac_L',
                   'parm_frac_epTB', 'parm_frac_lpTB', 'parm_frac_U', 'parm_frac_ATT')
  for (nm in parmfracnmz) if (parms[[nm]] < 1e-6) parms[[nm]] <- 1e-6 #safety
  totpop <- with(data = parms, {
    parm_frac_SD + parm_frac_CD+ parm_frac_E + parm_frac_L +
      parm_frac_epTB + parm_frac_lpTB+parm_frac_U+parm_frac_ATT
  })
  parms$parm_frac_SD <- parms$parm_frac_SD / totpop
  parms$parm_frac_CD <- parms$parm_frac_CD / totpop
  parms$parm_frac_E <- parms$parm_frac_E / totpop
  parms$parm_frac_L <- parms$parm_frac_L / totpop
  parms$parm_frac_epTB <- parms$parm_frac_epTB / totpop
  parms$parm_frac_lpTB <- parms$parm_frac_lpTB / totpop
  parms$parm_frac_U <- parms$parm_frac_U / totpop
  parms$parm_frac_ATT <- parms$parm_frac_ATT / totpop
  ## flow parms
  parms$parm_init_PPD <- c(
    smpsd[j]$R,       ## 1 = remand
    smpsd[j]$S,       ## 2 = short stay
    smpsd[j]$L,       ## 3 = long stay
    smpsd[j]$Omega,   ## 4 = open
    0.05              ## 5 = previously detained
  )
  parms$inflow <- smpsd[j]$E                     # inflow rate (less recidivists)
  parms$remand_short <- smpsd[j,ps*(1-pl)/dR]    # remand -> short:   1->2
  parms$remand_long <- smpsd[j,ps*pl/dR]         # remand -> long:    1->3
  parms$remand_release <- smpsd[j,(1-ps)*pl/dR]  # remand -> release: 1->5
  parms$long_short <- smpsd[j,tl]                # long -> short:     3->2
  parms$short_release <- smpsd[j,1/dS]           # short -> release:  2->5
  parms$short_open <- smpsd[j,ts]                # short -> open:     2->4
  parms$long_release <- smpsd[j,1/dL]            # long -> release:   3->5
  parms$open_release <- smpsd[j,1/dom]           # open -> release:   4->5
  parms$previous_remand <- 1e-3                  # previous -> remand 5->1 
  ## return
  parms
}

## TODO revise 
## revise initial states
revise.instates <- function(parms,tscale=20,rescale=TRUE){
  lam <- parms$foi
  ## utilities
  Delt <- 2.3 * lam + 0.1
  ome <- 1 - exp(-lam * tscale)
  ## fractions
  parm_frac_U <- exp(-lam * tscale)
  parm_frac_E <- ome * (2*lam) / Delt
  parm_frac_L <- ome * (0.1) / Delt
  parm_frac_ATT <- ome * (lam/30) / Delt
  parm_frac_SD <- ome * (lam / 30) / Delt
  parm_frac_CD <- ome * (lam / 30) / Delt
  parm_frac_epTB <- ome * (lam / 15) / Delt
  parm_frac_lpTB <- ome * (2 * lam / 15) / Delt
  totparm <- parm_frac_U+parm_frac_E+parm_frac_L+
    parm_frac_ATT+parm_frac_SD+parm_frac_CD+parm_frac_epTB+parm_frac_lpTB
  ## check
  ## print(parm_frac_U+parm_frac_E+parm_frac_L+parm_frac_ATT+
  ##       parm_frac_SD+parm_frac_CD+parm_frac_epTB+parm_frac_lpTB)
  ## implement in parameter object
  ## inflow & initial state
  parms$parm_frac_U <- parms$parm_ifrac_U <- parm_frac_U / totparm
  parms$parm_frac_E <- parms$parm_ifrac_E <- parm_frac_E / totparm
  parms$parm_frac_L <- parms$parm_ifrac_L <- parm_frac_L / totparm
  parms$parm_frac_ATT <- parms$parm_ifrac_ATT <- parm_frac_ATT / totparm
  parms$parm_frac_SD <- parms$parm_ifrac_SD <- parm_frac_SD / totparm
  parms$parm_frac_CD <- parms$parm_ifrac_CD <- parm_frac_CD / totparm
  parms$parm_frac_epTB <- parms$parm_ifrac_epTB <- parm_frac_epTB / totparm
  parms$parm_frac_lpTB <- parms$parm_ifrac_lpTB <- parm_frac_lpTB / totparm
  if(rescale){
    ## scale the inflow separately:
    tbip <- parms$tbip #
    tbip0 <- (1 - parm_frac_U / totparm)
    fac <- tbip / tbip0 #factor to apply to TBI+
    facn <- (1-fac * tbip0) / (1-tbip0) #factor to apply to TBI- st. total=1
    parms$parm_ifrac_U <- parms$parm_ifrac_U * facn
    parms$parm_ifrac_E <- parms$parm_ifrac_E * fac
    parms$parm_ifrac_L <- parms$parm_ifrac_L * fac
    parms$parm_ifrac_ATT <- parms$parm_ifrac_ATT * fac
    parms$parm_ifrac_SD <- parms$parm_ifrac_SD * fac
    parms$parm_ifrac_CD <- parms$parm_ifrac_CD * fac
    parms$parm_ifrac_epTB <- parms$parm_ifrac_epTB * fac
    parms$parm_ifrac_lpTB <- parms$parm_ifrac_lpTB * fac
    if( with(data=parms,{(parm_ifrac_U+
                          parm_ifrac_E+
                          parm_ifrac_L+
                          parm_ifrac_ATT+
                          parm_ifrac_SD+
                          parm_ifrac_CD+
                          parm_ifrac_epTB+
                          parm_ifrac_lpTB -1)^2})>1e-2) warning('Check scaling in revise instates!')
  }
  ## return
  parms
}


## adjust costs and se/sp to represent pre-screening
pre.screen <- function(parms,SE,SP,verbose=FALSE) {
  if (verbose) {
    print('------1: att')
    print(c(
      parms$inflow_toATT_TB1,
      parms$inflow_toATT_L1,
      parms$inflow_toATT_no1
    ))
    print("------1: tpt")
    print(c(
      parms$inflow_toTPT_TB1,
      parms$inflow_toTPT_L1,
      parms$inflow_toTPT_no1
    ))
    print("------1: notx")
    print(c(
      parms$inflow_notx_TB1,
      parms$inflow_notx_L1,
      parms$inflow_notx_no1
    ))
    print("------1: UC")
    print(c(
      parms$uc_entry_notx_TB,
      parms$uc_entry_notx_L,
      parms$uc_entry_notx_no
    ))
  }
  ## save originals for use later in UC adjustments:
  p_notx_TB1 <- parms$inflow_notx_TB1
  p_notx_L1 <- parms$inflow_notx_L1
  p_notx_no1 <- parms$inflow_notx_no1
  ## see notes:
  ## ... TBD
  parms$inflow_toATT_TB1 <- parms$inflow_toATT_TB1 * SE # se
  parms$inflow_toTPT_TB1 <- parms$inflow_toTPT_TB1 * SE # se
  parms$inflow_notx_TB1 <- parms$inflow_notx_TB1 * SE + (1-SE) # (1-sp)
  ## ...TBI
  parms$inflow_toATT_L1 <- parms$inflow_toATT_L1 * SE # se
  parms$inflow_toTPT_L1 <- parms$inflow_toTPT_L1 * SE # se
  parms$inflow_notx_L1 <- parms$inflow_notx_L1 * SE + (1-SE) # (1-sp)
  ## ...NO TB
  parms$inflow_toATT_no1 <- parms$inflow_toATT_no1 * (1 - SP) # 1-sp
  parms$inflow_toTPT_no1 <- parms$inflow_toTPT_no1 * (1 - SP) # 1-sp
  parms$inflow_notx_no1 <- parms$inflow_notx_no1 * (1 - SP) + SP# (1-sp)g + sp
  ## COSTS see notes: only notx differ
  parms$uc_entry_notx_TB <- parms$uc_entry_notx_TB * SE * p_notx_TB1 / parms$inflow_notx_TB1
  parms$uc_entry_notx_L <- parms$uc_entry_notx_L * SE * p_notx_L1 / parms$inflow_notx_L1
  parms$uc_entry_notx_no <- parms$uc_entry_notx_no * (1 - SP) * p_notx_no1 / parms$inflow_notx_no1
  if(verbose){
    print("------2: att")
    print(c(
      parms$inflow_toATT_TB1,
      parms$inflow_toATT_L1,
      parms$inflow_toATT_no1
    ))
    print("------2: tpt")
    print(c(
      parms$inflow_toTPT_TB1,
      parms$inflow_toTPT_L1,
      parms$inflow_toTPT_no1
    ))
    print("------2: notb")
    print(c(
      parms$inflow_notx_TB1,
      parms$inflow_notx_L1,
      parms$inflow_notx_no1
    ))
    print("------2: UC")
    print(c(
      parms$uc_entry_notx_TB,
      parms$uc_entry_notx_L,
      parms$uc_entry_notx_no
    ))
  }
  ## return
  parms
}


## ## check
## revise.instates(parms)

## revise HE stuff
revise.HE.parms <- function(parms, # original parameter template
                            DR,    # sample of outputs from tree
                            j,     # which row of sample to use
                            arm='soc', #soc/int
                            zero.nonscreen.costs=FALSE, #for debugging
                            summarize.HEparms=FALSE     #for recording
                            ){
  ## NOTE notx flow calculated from others in model
  ## NOTE costs only accrue post intervention in model so don't need 2 lots
  ## COLLATE DATA:
  A <- list(
    ## ATT for those found passively within the system
    uc_attppd = DR[id==j & tb=='TBD']$int_att_cost * !zero.nonscreen.costs, 
    uc_attout = 10e3 * !zero.nonscreen.costs, # ATT following release TODO what is the 10e3

    ## ------- BOTH
    ## ATT
    inflow_toATT_TB0 = DR[id==j & tb=='TBD']$soc_att_check, # NOTE fp ATT doesn't affect state
    inflow_toATT_L0 = DR[id==j & tb=='TBI']$soc_att_check, #
    inflow_toATT_no0 = DR[id==j & tb=='noTB']$soc_att_check,
    ## ...TPT
    inflow_toTPT_TB0 = DR[id==j & tb=='TBD']$soc_tpt_check,
    inflow_toTPT_L0 = DR[id==j & tb=='TBI']$soc_tpt_check,
    inflow_toTPT_no0 = DR[id==j & tb=='noTB']$soc_tpt_check,
    ## ...NOTX
    inflow_notx_TB0 = DR[id==j & tb=='TBD']$soc_notx_check,
    inflow_notx_L0 = DR[id==j & tb=='TBI']$soc_notx_check,
    inflow_notx_no0 = DR[id==j & tb=='noTB']$soc_notx_check,

    ## ------- SOC
    ## ... ATT
    inflow_toATT_TB1.soc = DR[id==j & tb=='TBD']$soc_att_check,
    inflow_toATT_L1.soc = DR[id==j & tb=='TBI']$soc_att_check,
    inflow_toATT_no1 = DR[id==j & tb=='noTB']$soc_att_check,
    ## ...TPT
    inflow_toTPT_TB1.soc = DR[id==j & tb=='TBD']$soc_tpt_check,
    inflow_toTPT_L1.soc = DR[id==j & tb=='TBI']$soc_tpt_check,
    inflow_toTPT_no1.soc = DR[id==j & tb=='noTB']$soc_tpt_check,
    ## ...NOTX
    inflow_notx_TB1.soc = DR[id==j & tb=='TBD']$soc_notx_check,
    inflow_notx_L1.soc = DR[id==j & tb=='TBI']$soc_notx_check,
    inflow_notx_no1.soc = DR[id==j & tb=='noTB']$soc_notx_check,
    ## ...unit costs
    ## ......TPT
    uc_entry_tpt_TB.soc = 0,## DR[id==j & tb=='TBD']$soc_tpt_cost
    uc_entry_tpt_L.soc = DR[id==j & tb=='TBI']$soc_tpt_cost,
    uc_entry_tpt_no.soc = 0,## DR[id==j & tb=='noTB']$soc_tpt_cost
    ## ......ATT
    uc_entry_att_TB.soc = DR[id==j & tb=='TBD']$soc_att_cost,
    uc_entry_att_L.soc = DR[id==j & tb=='TBI']$soc_att_cost,
    uc_entry_att_no.soc = DR[id==j & tb=='noTB']$soc_att_cost,
    ## ......NOTX
    uc_entry_notx_TB.soc = DR[id==j & tb=='TBD']$soc_notx_cost,
    uc_entry_notx_L.soc = DR[id==j & tb=='TBI']$soc_notx_cost,
    uc_entry_notx_no.soc = DR[id==j & tb=='noTB']$soc_notx_cost,

    ## ------- INT
    inflow_toATT_TB1.int = DR[id==j & tb=='TBD']$int_att_check,
    inflow_toATT_L1.int = DR[id==j & tb=='TBI']$int_att_check,
    inflow_toATT_no1.int = DR[id==j & tb=='noTB']$int_att_check,
    ## ...TPT
    inflow_toTPT_TB1.int = 0,## DR[id==j & tb=='TBD']$int_tpt_check,
    inflow_toTPT_L1.int = DR[id==j & tb=='TBI']$int_tpt_check,
    inflow_toTPT_no1.int = 0,## DR[id==j & tb=='noTB']$int_tpt_check,
    ## ...NOTX
    inflow_notx_TB1.int = DR[id==j & tb=='TBD']$int_notx_check,
    inflow_notx_L1.int = DR[id==j & tb=='TBI']$int_notx_check,
    inflow_notx_no1.int = DR[id==j & tb=='noTB']$int_notx_check,
    ## ...unit costs
    ## ......TPT
    uc_entry_tpt_TB.int = 0,## DR[id==j & tb=='TBD']$int_tpt_cost
    uc_entry_tpt_L.int = DR[id==j & tb=='TBI']$int_tpt_cost,
    uc_entry_tpt_no.int = 0,## DR[id==j & tb=='noTB']$int_tpt_cost
    ## ......ATT
    uc_entry_att_TB.int = DR[id==j & tb=='TBD']$int_att_cost,
    uc_entry_att_L.int = DR[id==j & tb=='TBI']$int_att_cost,
    uc_entry_att_no.int = DR[id==j & tb=='noTB']$int_att_cost,
    ## ......NOTX
    uc_entry_notx_TB.int = DR[id==j & tb=='TBD']$int_notx_cost,
    uc_entry_notx_L.int = DR[id==j & tb=='TBI']$int_notx_cost,
    uc_entry_notx_no.int = DR[id==j & tb=='noTB']$int_notx_cost
  )

  if(!summarize.HEparms){ #implement changes
    ## BOTH
    ## ... other unit costs
    parms$uc_attppd <- A$uc_attppd
    parms$uc_attout <- A$uc_attout
    ## ... HRQoL
    parms$hrqol <- A$hrqol
    ## ... ATT
    parms$inflow_toATT_TB0 <- A$inflow_toATT_TB0
    parms$inflow_toATT_L0 <- A$inflow_toATT_L0
    parms$inflow_toATT_no0 <- A$inflow_toATT_no0
    ## ...TPT
    parms$inflow_toTPT_TB0 <- A$inflow_toTPT_TB0
    parms$inflow_toTPT_L0 <- A$inflow_toTPT_L0
    parms$inflow_toTPT_no0 <- A$inflow_toTPT_no0
    ## INT
    if(arm=='soc'){
      ## ... ATT
      parms$inflow_toATT_TB1 <- A$inflow_toATT_TB1.soc
      parms$inflow_toATT_L1 <- A$inflow_toATT_L1.soc
      parms$inflow_toATT_no1 <- A$inflow_toATT_no1.soc
      ## ...TPT
      parms$inflow_toTPT_TB1 <- A$inflow_toTPT_TB1.soc
      parms$inflow_toTPT_L1 <- A$inflow_toTPT_L1.soc
      parms$inflow_toTPT_no1 <- A$inflow_toTPT_no1.soc
      ## ...NOTX
      parms$inflow_notx_TB1 <- A$inflow_notx_TB1.soc
      parms$inflow_notx_L1 <- A$inflow_notx_L1.soc
      parms$inflow_notx_no1 <- A$inflow_notx_no1.soc
      ## ...unit costs
      ## ......TPT
      parms$uc_entry_tpt_TB <- A$uc_entry_tpt_TB.soc
      parms$uc_entry_tpt_L <- A$uc_entry_tpt_L.soc
      parms$uc_entry_tpt_no <- A$uc_entry_tpt_no.soc
      ## ......ATT
      parms$uc_entry_att_TB <- A$uc_entry_att_TB.soc
      parms$uc_entry_att_L <- A$uc_entry_att_L.soc
      parms$uc_entry_att_no <- A$uc_entry_att_no.soc
      ## ......NOTX
      parms$uc_entry_notx_TB <- A$uc_entry_notx_TB.soc
      parms$uc_entry_notx_L <- A$uc_entry_notx_L.soc
      parms$uc_entry_notx_no <- A$uc_entry_notx_no.soc
    } else {
      ## ... ATT
      parms$inflow_toATT_TB1 <- A$inflow_toATT_TB1.int
      parms$inflow_toATT_L1 <- A$inflow_toATT_L1.int
      parms$inflow_toATT_no1 <- A$inflow_toATT_no1.int
      ## ...TPT
      parms$inflow_toTPT_TB1 <- A$inflow_toTPT_TB1.int
      parms$inflow_toTPT_L1 <- A$inflow_toTPT_L1.int
      parms$inflow_toTPT_no1 <- A$inflow_toTPT_no1.int
      ## ...NOTX
      parms$inflow_notx_TB1 <- A$inflow_notx_TB1.int
      parms$inflow_notx_L1 <- A$inflow_notx_L1.int
      parms$inflow_notx_no1 <- A$inflow_notx_no1.int
      ## ...unit costs
      ## ......TPT
      parms$uc_entry_tpt_TB <- A$uc_entry_tpt_TB.int
      parms$uc_entry_tpt_L <- A$uc_entry_tpt_L.int
      parms$uc_entry_tpt_no <- A$uc_entry_tpt_no.int
      ## ......ATT
      parms$uc_entry_att_TB <- A$uc_entry_att_TB.int
      parms$uc_entry_att_L <- A$uc_entry_att_L.int
      parms$uc_entry_att_no <- A$uc_entry_att_no.int
      ## ......NOTX
      parms$uc_entry_notx_TB <- A$uc_entry_notx_TB.int
      parms$uc_entry_notx_L <- A$uc_entry_notx_L.int
      parms$uc_entry_notx_no <- A$uc_entry_notx_no.int
    }

    return(parms)
  } else {
    return(A)
  }
}


## Other outputs:
others <- c('dLYL','deaths','notif100k','CC0','CC','qoldec')

## ## ============= diff plot  =============
## make comparison data for SOC vs INT
diffdata <- function(parms,SOCcov=0,INTcov=1){
  parms$inflow_notx_no0 <- parms$inflow_notx_no1 <- 0
  parms$inflow_toTPT_L0 <- parms$inflow_toATT_TB0 <- SOCcov # BL zero
  parms$inflow_toTPT_L1 <- parms$inflow_toATT_TB1 <- SOCcov # OFF
  tt <- seq(from=0, to=120, by=0.1)  #time frame to run over
  ## SOC:
  y <- runmodel(tt,parms)           #run model
  ## INT:
  parms$inflow_toTPT_L1 <- parms$inflow_toATT_TB1 <- INTcov #ON
  yi <- runmodel(tt,parms)           #run model
  ## convert
  ydm <- output2dt(y)
  ydmi <- output2dt(yi)
  ## join
  ydm[,arm:='SOC']
  ydmi[,arm:='INT']
  rbind(ydm,ydmi)
}

## get HE results for SOC/INT
run.HE.socint <- function(parms,DR,j,
                          zero.nonscreen.costs=FALSE, #for debugging/checking
                          ignore.tree.parms=FALSE,    #for debugging/checking
                          end_time=120,int_time=50,
                          static=-1,totpop=87489,prescreen=FALSE){
    tt <- seq(from=0, to=end_time, by=0.1)  #time frame to run over: 70 years after 50 burn
    parms$staticfoi <- static            #dynamic
    parms$int_time <- int_time #fix
    ## SOC: ## sample from tree-derived parms
    if(!ignore.tree.parms)
      parms <- revise.HE.parms(parms, DR, j, arm = "soc", zero.nonscreen.costs = zero.nonscreen.costs)
    test <- parms; test$staticfoi <- NULL
    if(any(unlist(test)<0)) stop(paste0('Parameter<0 for SOC @ run=',j,'\n parm=',
                                        names(test)[which(unlist(test)<0)],'\n'))
    ## y <- runmodel(tt,parms)           #run model
    YO <- runmodelsafely(tt, parms) # run model
    y <- YO$y
    ## INT: ## sample from tree-derived parms
    if (!ignore.tree.parms) {
      parms <- revise.HE.parms(parms, DR, j, arm = "int", zero.nonscreen.costs = zero.nonscreen.costs)
    }
    if (!prescreen) {
      SE1 <- 1
      SP1 <- 0
    } else {
      orv <- runif(n = 1, min = 1, max = 10) #OR
      SP1 <- runif(n = 1, min = 0.5, max = 1) #SP
      SE1 <- runif(n = 1, min = 0, max = 1) # SP
      ## oddssp <- odds(SP1)
      ## SE1 <- inv.odds(orv * oddssp)
      parms <- pre.screen(parms,SE1,SP1)
    }
    test <- parms; test$staticfoi <- NULL
    if(any(unlist(test)<0)) stop(paste0('Parameter<0 for INT @ run=',j,'\n parm=',
                                        names(test)[which(unlist(test)<0)],'\n'))
    YOI <- runmodelsafely(tt, parms) # run model
    yi <- YOI$y
    ## yi <- runmodel(tt,parms)           #run model
    mid <- which(y[,'t']==int_time)    #NOTE could break if don't fit exactly in units of dt
    end <- nrow(y)
    if(j==1)
      save(parms, file = here("transmission/data/first.parms.Rdata"))
    if(nrow(yi)<nrow(y)){
      cat("Fewer rows in yi (", nrow(yi), ") than y (", nrow(y), ")!\n")
      cat("j = ", j, "\n")
      cat("saving parameters!\n")
      save(parms, file = here("transmission/data/last.problem.parms.Rdata"))
      print(c(
        parms$inflow_toATT_TB1,
        parms$inflow_toATT_L1,
        parms$inflow_toATT_no1
      ))
      print(c(
        parms$inflow_toTPT_TB1,
        parms$inflow_toTPT_L1,
        parms$inflow_toTPT_no1
      ))
      print(c(
        parms$uc_entry_notx_TB,
        parms$uc_entry_notx_L,
        parms$uc_entry_notx_no
      ))
    }
    blpop <- y[mid,"ppdpop"]
    ratio <- totpop/blpop
    RES <- data.table(
      blpop=blpop,
      ratio=ratio,
      int.to.end=end_time-int_time,
      mid.notes=y[mid,"notif100k"],
      ## SOC
      soc.CC0=y[end,"CC0"], ## * ratio,
      soc.CC=y[end,"CC"], ## * ratio,
      soc.deaths=y[end,"deaths"],
      soc.qoldec=y[end,"qoldec"],
      soc.dLYL=y[end,"dLYL"],
      soc.ccases=y[end,"cases"],
      soc.ccasesout=y[end,"casesout"],
      soc.cTPT=y[end,'cTPT'],
      soc.cATTtp=y[end,'cATTtp'],
      soc.cATTfp=y[end,'cATTfp'],
      ## INT
      int.CC0=yi[end,"CC0"], ## * ratio,
      int.CC=yi[end,"CC"], ## * ratio,
      int.deaths=yi[end,"deaths"],
      int.qoldec=yi[end,"qoldec"],
      int.dLYL=yi[end,"dLYL"],
      int.ccases=yi[end,"cases"],
      int.ccasesout=yi[end,"casesout"],
      int.cTPT=yi[end,'cTPT'],
      int.cATTtp=yi[end,'cATTtp'],
      int.cATTfp = yi[end, "cATTfp"],
      ## pre-screen
      SE1=SE1,
      SP1=SP1
    )
    RES[,dQ:=(soc.qoldec+soc.dLYL-int.qoldec-int.dLYL)]
    RES[,Q.soc:=soc.qoldec+soc.dLYL]
    RES[,Q.int:=int.qoldec+int.dLYL]
    pbm <- ifelse(YO$error_val + YOI$error_val > 0, 1, 0) #record if problem
    RES[, problem := pbm]
    RES
}

## ## ============= HE workflow =============
PSAloop <- function(Niter=4e3,parms,smpsd,DR,
                    zero.nonscreen.costs=FALSE,verbose=FALSE,
                    static=FALSE,community=TRUE,posttb=TRUE,targeting=TRUE){
  if(!static & community & posttb) cat('Running basecase analysis!\n')
  if (static) cat("Running static model!\n")
  if (!community) cat("Running with community transmission turned off!\n")
  if (!posttb) cat("Running with post-TB effects turned off!\n")
  if (targeting) cat("Running with pre-screen targeting turned on!\n")
  if(Niter>nrow(smpsd)){
    cat('Niter>nrow(smpsd): resampling extra replicates!\n')
    xtra <- smpsd[sample(nrow(smpsd),Niter-nrow(smpsd),replace=TRUE)]
    smpsd <- rbind(smpsd,xtra)
  }
  if(Niter>max(DR$id)){
    cat('Niter>max(DR$id): resampling extra replicates!\n')
    reuse <- sample(max(DR$id),Niter-max(DR$id),replace=TRUE)
    nxt <- max(DR$id)+1
    xtrad <- DR[tb=='TBD'][reuse]; xtrad[,id:=nxt:Niter]
    xtrai <- DR[tb=='TBI'][reuse]; xtrai[,id:=nxt:Niter]
    xtran <- DR[tb=='noTB'][reuse]; xtran[,id:=nxt:Niter]
    DR <- rbindlist(list(DR,xtrad,xtrai,xtran))
  }
  ## loop
  RES <- list()
  for(j in 1:Niter){
    if(verbose) cat('j==',j,'...\n')
    if(!j%%50) print(j)
    ## set parms related to PPD flows
    parms <- revise.flow.parms(parms,smpsd,j)
    ## record for later (redone in run.HE.socint):
    hep <- revise.HE.parms(parms,DR,j,
                           zero.nonscreen.costs=zero.nonscreen.costs,summarize.HEparms = TRUE)
    ## set intervention and HE parms
    ## sample from TB natural hist
    newp <- uv2ps(runif(length(hyperparms)),hyperparms) # natural history
    newp[["m"]] <- 1 + newp[["m"]] # NOTE prior is defined for R-1
    ## safeties:
    newp[["CDR"]] <- min(0.9, newp[["CDR"]])
    newp[["mHR"]] <- max(1, newp[["mHR"]]) #needs truncating to avoid -ve deaths
    newp[["wsn"]] <- max(0.2, newp[["wsn"]])
    ## ## BUG checking
    ## newp[["mHR"]] <- 1
    ## newp[["hrqolptb"]] <- 0
    for(nm in names(newp)) parms[[nm]] <- newp[[nm]] #safety
    ## update initial state:
    parms <- revise.instates(parms)
    ## sensitivity analyses
    parms$staticfoi <- ifelse(static,+1,-1)
    if (!community) parms$m <- 1 #this has 1 added to it
    if (!posttb){
      parms$mHR <- 1
      parms$hrqolptb <- 0
    }
    ## === run model:
    ANS <- run.HE.socint(parms,
                         DR,j,
                         zero.nonscreen.costs=zero.nonscreen.costs,
                         prescreen = targeting)
    ## capture parms also
    V <- unlist(parms)
    nm <- names(V)
    V <- as.list(V)
    ## take out HE and put back in
    hepn <- unique(gsub('\\.int|\\.soc','',names(hep)))
    for(n in names(V)) if(n %in% hepn) V[[n]] <- NULL #remove HE vars
    V <- c(V,hep)                                     #put back in
    nm <- names(V)
    ANS[, c(nm) := V]                                 #record in answers
    ## record
    RES[[j]] <- ANS
  } #end loop
  RES <- rbindlist(RES)
  ## add pre-screen data
  RES[, prevTBI := (parm_ifrac_L + parm_ifrac_E + parm_ifrac_epTB + parm_ifrac_lpTB)]
  ## RES[, SP1.cat := cut(SP1, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), include.lowest = TRUE)]
  ## RES[, SE1.cat := cut(SE1, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), include.lowest = TRUE)]
  RES[, icer := (int.CC - soc.CC) / dQ]
  ## inspect
  cat('--- deaths diff summary ---\n')
  print(RES[,summary(soc.deaths-int.deaths)])
  cat('--- cost diff summary ---\n')
  print(RES[,summary(int.CC0-soc.CC0)])
  if(any(RES$problem>0)) warnings('Some ODE runs had numerical issues! Please check results')
  ## return
  RES
}


## =============== Authors only:
library(googlesheets4) # only for authors to upload

yourl <- "https://docs.google.com/spreadsheets/d/1XYIWB4NAYUzLNESdLBDGqviXcDP_s0RIIDqmQKa5T7I/edit?gid=1037492332#gid=1037492332"

shid <- as.character(as_sheets_id(yourl))

## utility function
upload.to.sheets <- function(basename, filename, sheetid) {
  filename <- gsub("\\.csv$", "", filename) # safety in case csv included at and
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file = fn)
  ## NOTE nonsense to deal with potential Inf error in googlesheets
  for (j in seq_len(ncol(tmp))) {
    if(class(tmp[[j]])=='numeric')
      set(tmp, which(!is.finite(tmp[[j]])), j, -1)
  }
  write_sheet(tmp, sheetid, sheet = filename)
}
