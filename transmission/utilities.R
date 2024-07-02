## utility/helper functions
library(ecrins) #install from petedodd/ecrins
library(data.table)
library(stringr)
library(ggplot2)
library(scales)


source(system.file("extdata", "parameters.R", package="ecrins")) #this loads some fns for PSA
data(parms)                       #some default parameters

## prisons parameters
## parms
load(here('churn/outdata/smpsd.Rdata')) #TODO change/update

## tree outputs
load(here('outdata/DR.Rdata'))
DRS <- fread(here('outdata/DRS.csv'))


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

revise.HE.parms <- function(parms, # original parameter template
                            DR,    # sample of outputs from tree
                            j,     # which row of sample to use
                            arm='soc' #soc/int
                            ){

  ## NOTE notx flow calculated from others in model
  ## NOTE costs only accrue post intervention in model so don't need 2 lots
  ## SOC
  ## ... ATT
  parms$inflow_toATT_TB0 <- DR[id==j & tb=='TBD']$soc_att_check # NOTE fp ATT doesn't affect state
  parms$inflow_toATT_L0 <- DR[id==j & tb=='TBI']$soc_att_check #
  parms$inflow_toATT_no0 <- DR[id==j & tb=='noTB']$soc_att_check
  ## ...TPT
  parms$inflow_toTPT_TB0 <- DR[id==j & tb=='TBD']$soc_tpt_check
  parms$inflow_toTPT_L0 <- DR[id==j & tb=='TBI']$soc_tpt_check
  parms$inflow_toTPT_no0 <- DR[id==j & tb=='noTB']$soc_tpt_check
  ## INT
  if(arm=='soc'){
    ## ... ATT
    parms$inflow_toATT_TB1 <- DR[id==j & tb=='TBD']$soc_att_check
    parms$inflow_toATT_L1 <- DR[id==j & tb=='TBI']$soc_att_check
    parms$inflow_toATT_no1 <- DR[id==j & tb=='noTB']$soc_att_check
    ## ...TPT
    parms$inflow_toTPT_TB1 <- DR[id==j & tb=='TBD']$soc_tpt_check
    parms$inflow_toTPT_L1 <- DR[id==j & tb=='TBI']$soc_tpt_check
    parms$inflow_toTPT_no1 <- DR[id==j & tb=='noTB']$soc_tpt_check
    ## ...unit costs
    ## ......TPT
    parms$uc_entry_tpt_TB <- 0## DR[id==j & tb=='TBD']$soc_tpt_cost
    parms$uc_entry_tpt_L <- DR[id==j & tb=='TBI']$soc_tpt_cost
    parms$uc_entry_tpt_no <- 0## DR[id==j & tb=='noTB']$soc_tpt_cost
    ## ......ATT
    parms$uc_entry_att_TB <- DR[id==j & tb=='TBD']$soc_att_cost
    parms$uc_entry_att_L <- DR[id==j & tb=='TBI']$soc_att_cost
    parms$uc_entry_att_no <- DR[id==j & tb=='noTB']$soc_att_cost
    ## ......NOTX
    parms$uc_entry_notx_TB <- DR[id==j & tb=='TBD']$soc_notx_cost
    parms$uc_entry_notx_L <- DR[id==j & tb=='TBI']$soc_notx_cost
    parms$uc_entry_notx_no <- DR[id==j & tb=='noTB']$soc_notx_cost
  } else {
    ## ... ATT
    parms$inflow_toATT_TB1 <- DR[id==j & tb=='TBD']$int_att_check
    parms$inflow_toATT_L1 <- DR[id==j & tb=='TBI']$int_att_check
    parms$inflow_toATT_no1 <- DR[id==j & tb=='noTB']$int_att_check
    ## ...TPT
    parms$inflow_toTPT_TB1 <- 0## DR[id==j & tb=='TBD']$int_tpt_check
    parms$inflow_toTPT_L1 <- DR[id==j & tb=='TBI']$int_tpt_check
    parms$inflow_toTPT_no1 <- 0## DR[id==j & tb=='noTB']$int_tpt_check
    ## ...unit costs
    ## ......TPT
    parms$uc_entry_tpt_TB <- 0## DR[id==j & tb=='TBD']$int_tpt_cost
    parms$uc_entry_tpt_L <- DR[id==j & tb=='TBI']$int_tpt_cost
    parms$uc_entry_tpt_no <- 0## DR[id==j & tb=='noTB']$int_tpt_cost
    ## ......ATT
    parms$uc_entry_att_TB <- DR[id==j & tb=='TBD']$int_att_cost
    parms$uc_entry_att_L <- DR[id==j & tb=='TBI']$int_att_cost
    parms$uc_entry_att_no <- DR[id==j & tb=='noTB']$int_att_cost
    ## ......NOTX
    parms$uc_entry_notx_TB <- DR[id==j & tb=='TBD']$int_notx_cost
    parms$uc_entry_notx_L <- DR[id==j & tb=='TBI']$int_notx_cost
    parms$uc_entry_notx_no <- DR[id==j & tb=='noTB']$int_notx_cost
  }

  ## other unit costs
  parms$uc_attppd <- 20e3 # ATT for those found passively within the system
  parms$uc_attout <- 15e3 # ATT following release

  ## === HRQoL
  parms$hrqol <- 0.3 # HRQoL decrement while CD
  parms$hrqolptb <- 0.05 # HRQoL decrement while post TB
  parms$m <- 1.0         #multiplier for TB events outside prison

  ## return
  parms
}


## Other outputs:
others <- c('dLYL','deaths','notif100k','CC0','CC','qoldec')

## ## ============= diff plot  =============
## make comparison data for SOC vs INT
diffdata <- function(parms,SOCcov=0,INTcov=1){
  parms$inflow_toTPT_L0 <- parms$inflow_toATT_TB0 <- SOCcov #BL zero
  parms$inflow_toTPT_L1 <- parms$inflow_toATT_TB1 <- SOCcov #OFF
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

## parameter meanings [0/1=BL/INT]:
## inflow_toTPT_L 
## inflow_toATT_TB
## are both coverages of the intervention
## TODO se/sp to come
## TODO check destinations

## ## ============= HE workflow =============
PSAloop <- function(Niter=4e3,parms,smpsd,DR){
  if(Niter>nrow(smpsd)){
    cat('Niter>nrow(smpsd): resampling extra replicates!\n')
    xtra <- smpsd[sample(nrow(smpsd),Niter-nrow(smpsd),replace=TRUE)]
    smpsd <- rbind(smpsd,xtra)
  }
  if(Niter>max(DR$id)){
    ## TODO check this works OK
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
    if(!j%%50) print(j)
    ## set parms related to PPD flows
    parms <- revise.flow.parms(parms,smpsd,j)
    ## set intervention and HE parms
    ## sample from TB natural hist
    newp <- uv2ps(runif(length(hyperparms)),hyperparms) # natural history
    for(nm in names(newp)) parms[[nm]] <- newp[[nm]]
    ## update initial state:
    ## TODO
    ## === run model:
    tt <- seq(from=0, to=120, by=0.1)  #time frame to run over: 70 years after 50 burn
    parms$staticfoi <- -1            #dynamic
    parms$int_time <- 50             #fix
    ## parms$inflow_toTPT_L0 <- parms$inflow_toATT_TB0 <- 0 #BL zero: se
    ## parms$inflow_toTPT_L1 <- parms$inflow_toATT_TB1 <- 0 #OFF: se
    ## SOC:
    parms <- revise.HE.parms(parms,DR,j,arm='soc') ## sample from tree-derived parms
    y <- runmodel(tt,parms)           #run model
    ## INT:
    ## parms$inflow_toTPT_L1 <- parms$inflow_toATT_TB1 <- 1 #ON later : se
    parms <- revise.HE.parms(parms,DR,j,arm='int') ## sample from tree-derived parms
    yi <- runmodel(tt,parms)           #run model
    ## ydm <- output2dt(y)
    mid <- which(y[,'t']==50)
    end <- nrow(y)
    blpop <- y[mid,"ppdpop"]
    totpop <-  87489
    ratio <- totpop/blpop
    RES[[j]] <- data.table(
      blpop=blpop,
      mid.notes=y[mid,"notif100k"],
      ## SOC
      soc.CC0=y[end,"CC0"] * ratio,
      soc.CC=y[end,"CC"] * ratio,
      soc.deaths=y[end,"deaths"],
      soc.qoldec=y[end,"qoldec"],
      soc.dLYL=y[end,"dLYL"],
      ## INT
      int.CC0=yi[end,"CC0"] * ratio,
      int.CC=yi[end,"CC"] * ratio,
      int.deaths=yi[end,"deaths"],
      int.qoldec=yi[end,"qoldec"],
      int.dLYL=yi[end,"dLYL"]
    )
  } #end loop
  RES <- rbindlist(RES)
  RES[,dQ:=(soc.qoldec+soc.dLYL-int.qoldec-int.dLYL)]
  RES[,Q.soc:=soc.qoldec+soc.dLYL]
  RES[,Q.int:=int.qoldec+int.dLYL]
  ## inspect
  cat('--- deaths diff summary ---\n')
  print(RES[,summary(soc.deaths-int.deaths)])
  cat('--- cost diff summary ---\n')
  print(RES[,summary(int.CC0-soc.CC0)])
  ## return
  RES
}


## TODO
## check natural history
## introduce the initial state heuristic as per methods doc
## check FPs!!
## what are our targets
