## TODO questions:
## 
## - stool/sputum
## - flag assumption = groups for SA or pending data
## 

## ========= UTILITIES ===============
logit <- function(x) log(odds(x))
ilogit <- function(x) iodds(exp(x))
AOR <- function(base,OR) ilogit(log(OR) + logit(base))
AOR2 <- function(base,OR) OR*base/(1-base+(base*OR))
odds <- function(x) x/(1-x)
iodds <- function(x) x/(1+x)
lo <- function(x) quantile(x,probs = 0.025, na.rm=TRUE)
hi <- function(x) quantile(x,probs = 1-0.025, na.rm=TRUE)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

brkt <- function(M,L,H,ndp=0) paste0(round(M,ndp),' (',
                                     round(L,ndp),' - ',
                                     round(H,ndp),')')
gm <- function(x) exp(mean(log(x))) #geometric mean
gh <- function(x) glue(here(x))

## ========= OUTCOMES ===============
## TODO - remove excess RNG here
CFRtxY <- function(age,hiv=0,art=0){#NB optimized for clarity not speed
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$ontx.u5$r(length(age))
  tmp[age=='5-14'] <- P$ontx.o5$r(sum(age=='5-14'))  #NB this could be achieved in the tree model
  ## hivartOR
  Z <- P$hivartOR$r(length(age))
  hor <- rep(1,length(age))
  tmp <- logit(tmp)                     #transform
  tmp[hiv>0] <- tmp[hiv>0]+Z[hiv>0,1]
  tmp[art>0] <- tmp[art>0]+Z[art>0,2]
  tmp <- ilogit(tmp)                    #inverse transform
  tmp
}
## CFRtxY(1:10,P)                            #test
## summary(CFRtxY(1:1e3,P))
## summary(CFRtxY(1:1e3,hiv=1))
## summary(CFRtxY(1:1e3,hiv=1,art=1,P))


## == CFR off tx
CFRtxN <- function(age,hiv=0,art=0){
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$notx.u5$r(length(age))          #default a<5 and hiv=art=0
  tmp[age!='5-14' & hiv>0 & art==0] <- P$notxH.u5$r(sum(age!='5-14' & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[age!='5-14' & hiv>0 & art>0] <- P$notxHA.u5$r(sum(age!='5-14' & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[age=='5-14'] <- P$notx.o5$r(sum(age=='5-14'))    #o5, HIV-ve
  tmp[age=='5-14' & hiv>0 & art==0] <- P$notxH.o5$r(sum(age=='5-14' & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[age=='5-14' & hiv>0 & art>0] <- P$notxHA.o5$r(sum(age=='5-14' & hiv>0 & art>0)) #o5,HIV+,ART+
  tmp
}
## CFRtxN(1:10,P)                            #test
## summary(CFRtxN(1:1e3,P))
## summary(CFRtxN(1:1e3,hiv=1,P))
## summary(CFRtxN(1:1e3,hiv=1,art=1,P))

## add CFRs to data by side-effect
AddCFRs <- function(D,P){
  ## d.cfr.notx & d.cfr.tx
  D[,c('cfr.notx','cfr.tx'):=0] #NOTE neglect non-TB mortality
  ## CFR on  ATT
  D[,cfr.tx:=CFRtxY(age,hiv,art)]
  ## CFR w/o ATT
  D[,cfr.notx:=CFRtxN(age,hiv,art)]
}

## ======= COMBINED LABELLER ===========

## additional labels from data (may overwrite some initial version currently)
AddDataDrivenLabels <- function(D){
  
  # Setting some SOC/INT parameters by TB status
  # For ATT: TB symptoms -> Clinical suspicion -> abnormal x-ray -> GP assessment -> tb dx
  # For TPT: IGRA tested -> IGRA positive
  # For setting SOC parameters
  # assuming a baseline some fraction lower than INT for TB screening & prop.xray
  # could also use any TB symptom sensitivity/specificity for TB screening assuming no xray in SOC
  D[,soc.fac:=rbeta(nrow(D),5/1.5,10)]
  # curve(dbeta(x,5/3,10),from=0,to=1,n=200)
  # summary(D$soc.fac)
  # D[,soc.fac:=1]
  
  # TB symptom screening
  D[,soc.prop.tb.sympt.screen:=int.prop.tb.sympt.screen*soc.fac]
  
  # Previous TB diagnosis
  #' `Uncomment these two lines to switched off 'previous TB diagnosis'`
  # D[,soc.prop.prev.tb.dx:=0]
  # D[,int.prop.prev.tb.dx:=0]
  # 
  D[,int.prop.prev.tb.dx:=ifelse(tb=='TBD',int.prop.prev.tb.dx,0)]
  D[,soc.prop.prev.tb.dx:=ifelse(tb=='TBD',soc.prop.prev.tb.dx,0)]
  
  # --------------------------------------------------------------------------------------
  # TB presumption/TB positive symptom screening 
  D[,soc.prop.presumtive.tb:=ifelse(tb=='TBD',sens.symptom,1-spec.symptom)]
  D[,int.prop.presumtive.tb:=ifelse(tb=='TBD',sens.any.abn.xray,1-spec.any.abn.xray)]
  
  # Fraction getting a chest x-ray under SOC
  D[,soc.prop.xray:=int.prop.xray*soc.fac]
  
  # Diagnostic accuracy of Xpert MTB/RIF for pulmonary TB in adults
  D[,soc.prop.prev.tb.tx.symp.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  D[,soc.prop.abn.xray.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  D[,soc.prop.no.xray.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  D[,soc.prop.no.prev.tb.dx.symp.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  
  D[,int.prop.prev.tb.tx.symp.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  D[,int.prop.abn.xray.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  D[,int.prop.no.xray.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  D[,int.prop.no.prev.tb.dx.symp.tb.dx:=ifelse(tb=='TBD',sens.xpert,1-spec.xpert)]
  
  # starting ATT under SOC
  D[,soc.prop.starting.att:=int.prop.starting.att]
  
  # LTBI pathway stuff
  # For TPT: IGRA tested -> IGRA positive -> stying 3+ months -> starting TPT
  
  # IGRA tested 
  #' `set int.prop.igra.tested to 0 switch off LTBI pathway`
  # D[,int.prop.igra.tested:=0] # Temporary assignment to 1
  D[,soc.prop.igra.tested:=int.prop.igra.tested*soc.fac]
  
  # D[,cost.igra.test:=35] # temoporary assignment to 50
  # D[,int.prop.igra.test.positive:=0.15]
  # D[,soc.prop.igra.test.positive:=int.prop.igra.test.positive*soc.fac]
  
  D[,soc.prop.igra.tested:=ifelse(tb=='TBI',soc.prop.igra.tested,0)]
  D[,int.prop.igra.tested:=ifelse(tb=='TBI',int.prop.igra.tested,0)]
  
  # IGRA positive
  # D[,sens.igra:=0.89] # Not being used
  # D[,spec.igra:=0.98]
  # D[,soc.prop.igra.test.positive:=ifelse(tb=='TBI',soc.prop.igra.test.positive,0)]
  # D[,int.prop.igra.test.positive:=ifelse(tb=='TBI',int.prop.igra.test.positive,0)]
  
  D[,soc.prop.igra.test.positive:=ifelse(tb=='TBI',1,0)]
  D[,int.prop.igra.test.positive:=ifelse(tb=='TBI',1,0)]
  
  # (check <- names(D)[grepl('staying.o3.months.tpt|staying.u3.months.uk.tpt', names(D))])
  # summary(D[,..check])
  #
  # The rest of these probably don't need to be changed
  # # staying longer that 3 months
  # D[,int.prop.staying.o3.months:=1] # Temporary assignment to 1
  # D[,soc.prop.staying.o3.months:=1]
  D[,int.prop.staying.o3.months:=ifelse(tb=='TBI',int.prop.staying.o3.months,0)]
  D[,soc.prop.staying.o3.months:=ifelse(tb=='TBI',soc.prop.staying.o3.months,0)]
  
  # #remaining in the UK
  # D[,int.prop.staying.u3.months.uk:=1] # Temporary assignment to 1
  # D[,soc.prop.staying.u3.months.uk:=int.prop.staying.u3.months.uk*soc.fac]
  D[,int.prop.staying.u3.months.uk:=ifelse(tb=='TBI',int.prop.staying.u3.months.uk,0)]
  D[,soc.prop.staying.u3.months.uk:=ifelse(tb=='TBI',soc.prop.staying.u3.months.uk,0)]
  
  # # referral to TPT
  # Nolonger needed
  # D[,tpt.nhs.referral:=1]
  
  # 
  # (check <- names(D)[grepl('int.prop.staying.o3.months.tpt', names(D))])
  # summary(D[,..check])
  # # starting TPT
  # D[,int.prop.staying.o3.months.tpt:=1] # Temporary assignment to 1
  # D[,soc.prop.staying.o3.months.tpt:=int.prop.staying.o3.months.tpt]
  # D[,int.prop.staying.o3.months.tpt:=ifelse(tb=='TBI',int.prop.staying.o3.months.tpt,0)]
  # D[,soc.prop.staying.o3.months.tpt:=ifelse(tb=='TBI',soc.prop.staying.o3.months.tpt,0)]
  # D[,int.prop.staying.u3.months.uk.tpt:=ifelse(tb=='TBI',int.prop.staying.u3.months.uk.tpt,0)]
  # D[,soc.prop.staying.u3.months.uk.tpt:=ifelse(tb=='TBI',soc.prop.staying.u3.months.uk.tpt,0)]
  # 
  # # complete TPT
  # D[,int.prop.staying.o3.months.complete.tpt:=0.99] # Temporary assignment to 1
  # D[,soc.prop.staying.o3.months.complete.tpt:=int.prop.staying.o3.months.complete.tpt]
  # D[,soc.prop.staying.o3.months.complete.tpt:=ifelse(tb=='TBI',soc.prop.staying.o3.months.complete.tpt,0)]
  # D[,int.prop.staying.o3.months.complete.tpt:=ifelse(tb=='TBI',int.prop.staying.o3.months.complete.tpt,0)]
  # 
  # D[,int.prop.staying.u3.months.uk.complete.tpt:=1] # Temporary assignment to 1
  # D[,soc.prop.staying.u3.months.uk.complete.tpt:=int.prop.staying.u3.months.uk.complete.tpt*soc.fac]
  # D[,soc.prop.staying.u3.months.uk.complete.tpt:=ifelse(tb=='TBI',soc.prop.staying.u3.months.uk.complete.tpt,0)]
  # D[,int.prop.staying.u3.months.uk.complete.tpt:=ifelse(tb=='TBI',int.prop.staying.u3.months.uk.complete.tpt,0)]
  # 
  
  # Just a few more changes/naming
  # names(D)[grepl('.dots', names(D))]
  D[,cost.chest.xray:=ucost.chest.xray+ucost.prison.escort]
  D[,cost.prison.gp.assessessment:=pIsolation*ucost.prison.cell.isolation + ucost.prison.gp.assess + int.prop.xray*cost.chest.xray]
  D[,cost.contact.management:=pIsolation*ucost.prison.cell.isolation + nContacts*ucost.contact.tracing]
  D[,cost.tb.evaluation:=ucost.nhs.tb.service + ucost.prison.escort + ucost.tb.investigations]
  D[,cost.attending.nhs.tb.service:=ucost.nhs.tb.service + ucost.prison.escort]
  D[,cost.att.initiation:=ucost.nhs.tb.service + ucost.prison.escort]
  D[,cost.tpt.initiation:=ucost.nhs.tb.service + ucost.prison.escort]
  D[,cost.prison.isolation:=pIsolation*ucost.prison.cell.isolation]
  
  D[,ucost.dots:=ucost.att.dots]
  D[,ucost.att.dots:=NULL]
  D[,ucost.tpt.opd.visit:=ucost.dstb.opd.visit]
  # D[,ucost.mdrtb.opd.visits:=ucost.mdrtb.opd.visit]
  # D[,ucost.opd.visit:=ucost.dstb.opd.visit]
  # D[,ucost.dstb.opd.visit:=NULL]
  D[,durTPT:=DurTPT]
  # D[,cost.inpatient:=(pDSTB*smear.positive*DurDSTBIsolation*(cost.dstb.ipd + cost.prison.bedwatch) + 
  #                       (1-pDSTB)*smear.positive*DurMDRTBIsolation*(cost.mdrtb.ipd.smear.positive + cost.prison.bedwatch))]
  D[,cost.inpatient:=(pDSTB*smear.positive*DurDSTBIsolation*(ucost.dstb.ipd + ucost.prison.bedwatch) + 
                        (1-pDSTB)*smear.positive*DurMDRTBIsolation*(ucost.mdrtb.ipd + ucost.prison.bedwatch))]
  D[,cost.att.complete:=pDSTB*dstb.visits*(ucost.dstb.opd.visit + ucost.prison.escort) + # DSTB outpatient visits
      (1-pDSTB)*mdrtb.visits*(ucost.mdrtb.opd.visit + ucost.prison.escort) + # MDRTB outpatient visits
      pDSTB*DurDSTB*ucost.dsatt.drugs + # DSTB drugs
      (1-pDSTB)*DurMDRTB*ucost.mdratt.drugs + # MDRTB drugs
      ucost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) + cost.inpatient]
  D[,cost.att.incomplete:= IncompDurDSTB/DurDSTB*(
    pDSTB*dstb.visits*(ucost.dstb.opd.visit + ucost.prison.escort) + 
      pDSTB*DurDSTB*(ucost.dsatt.drugs + ucost.dots)) + 
      IncompDurMDRTB/DurMDRTB*(
        (1-pDSTB)*mdrtb.visits*(ucost.mdrtb.opd.visit + ucost.prison.escort)  + 
          (1-pDSTB)*DurMDRTB*(ucost.mdratt.drugs + ucost.dots)) + 
      cost.inpatient]
  D[,cost.tpt.complete:=durTPT*(ucost.ltbi.drugs + ucost.dots) + TPT.visits*(ucost.tpt.opd.visit + ucost.prison.escort)]
  D[,cost.tpt.incomplete:=IncompDurTPT/durTPT*(durTPT*(ucost.ltbi.drugs + ucost.dots) + TPT.visits*(ucost.tpt.opd.visit + ucost.prison.escort))]
  
  # D[,pAttending:=1] # changed to attend.nhs.referral
  D[,cost.tb.sympt.screen:=ucost.tb.sympt.screen+ucost.overheads]
  D[,cost.igra.test:=ucost.igra.test]
  # D[,cost.chest.xray:=ucost.tst.test]
  D[,pIsolation:=0]
  D[,nhs.referral:=1]
  D[,soc.fac:=NULL] # remove temporary variable
  
}


## combined function to add the labels to the tree prior to calculations
MakeTreeParms <- function(D,P){
  ## -- use of other functions
  # AddSampleTests(D) #samples/tests
  # AddCFRs(D,P) #outcomes
  # AddTPTrr(D,P)
  # AddProgProb(D, P)
  # ## new labels from data
  AddDataDrivenLabels(D)
  # AddDetectionLabels(D)
}

## ======= EPIDEMIOLOGY ===========

makeAttributes <- function(D){
    nrep <- nrow(D)
    D[,id:=1:nrep]
    fx <- list(age=agelevels, isoz=c('GBR'), tb=tblevels)
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
  
    # ## --- age
    # D[,value:=ifelse(age=='15-64',0.8,1-0.8)] #NOTE value first set here
    D[,value:=1] # Not using age splits for now
    # ## --- HIV/ART
    # D[,h01:=0]
    # D[age!='5-14',h10:=hivprev.u5*(1-artcov)]
    # D[age=='5-14',h10:=hivprev.o5*(1-artcov)]
    # D[age!='5-14',h00:=1-hivprev.u5]
    # D[age=='5-14',h00:=1-hivprev.o5]
    # D[age!='5-14',h11:=hivprev.u5*artcov]
    # D[age=='5-14',h11:=hivprev.o5*artcov]
    # D[hiv==0 & art==0,value:=value*h00]
    # D[hiv==0 & art==1,value:=value*h01]
    # D[hiv==1 & art==0,value:=value*h10]
    # D[hiv==1 & art==1,value:=value*h11]
    # D[,c('h00','h01','h10','h11'):=NULL]
    ## --- TB
    ## ## Using 3 levels of TB: active TB disease, TB infection, and noTB (no active TB/TB infection)
    # D[,tbi:=tb.prev] # not being used for now
    D[,ltbi.prev:=ltbi.prev]
    D[,prog.tb:=prog.tb]

    D[tb=='TBD',value:=value*ltbi.prev*prog.tb]
    D[tb=='TBI',value:=value*(ltbi.prev-ltbi.prev*prog.tb)]
    D[tb=='noTB',value:=value*(1-ltbi.prev)] # 
    # D[,tbi:=NULL]                            #remove temporary variable
    return(D)
}

unique(vrz.soc[grepl('prop.', vrz.soc)])

# ## Toggle intervention
# # Just changing cascade parameters that are likely to change under intervention
# SetIntervention <- function(D,arm='SOC'){
#   if(arm=='SOC'){
#     D[,prop.tb.sympt.screen:=prop.tb.sympt.screen/10]
#     D[,prop.xray:=prop.xray/10]
#     D[,prop.igra.tested:=prop.igra.tested/100]
#     
#     D[,prop.prev.tb.tx.symp.tb.dx:=ifelse(tb=='TB',sens.symptom,1-spec.symptom)]
#     D[,prop.abn.xray.tb.dx:=ifelse(tb=='TB',sens.symptom,1-spec.symptom)]
#     D[,prop.no.xray.tb.dx:=ifelse(tb=='TB',sens.symptom,1-spec.symptom)]
#     D[,prop.prev.tb.dx:=ifelse(tb=='TB',sens.symptom,1-spec.symptom)]
#     D[,prop.no.prev.tb.dx.symp.tb.dx:=ifelse(tb=='TB',sens.symptom,1-spec.symptom)]
#   } else {
#     D[,prop.tb.sympt.screen:=prop.tb.sympt.screen]
#     D[,prop.xray:=prop.xray]
#     D[,prop.igra.tested:=prop.igra.tested]
#     
#     D[,prop.prev.tb.tx.symp.tb.dx:=ifelse(tb=='TB',sens.any.abn.xray,1-spec.any.abn.xray)]
#     D[,prop.abn.xray.tb.dx:=ifelse(tb=='TB',sens.any.abn.xray,1-spec.any.abn.xray)]
#     D[,prop.no.xray.tb.dx:=ifelse(tb=='TB',sens.any.abn.xray,1-spec.any.abn.xray)]
#     D[,prop.prev.tb.dx:=ifelse(tb=='TB',sens.any.abn.xray,1-spec.any.abn.xray)]
#     D[,prop.no.prev.tb.dx.symp.tb.dx:=ifelse(tb=='TB',sens.any.abn.xray,1-spec.any.abn.xray)]
#   }
# }

## function for generating random sample of costs
MakeCostData <- function(csts,          #base data table of cost data
                         nrep,          #number of replicates being used in PSA
                         anmz=NULL     #attribute names (if any)
                         ){
  if(nrow(csts[cost.sd>0 & cost.m==0])>0) warning(paste0('Some cost input variables have zero mean & SD>0. These will be treated as fixed variables:\n',paste0(csts[cost.sd>0 & cost.m==0,cost],collapse='\n')))
  if(is.null(anmz)& any(csts[,table(cost)]>1)) warning('Some cost names occur >1 times, but no attributes have been specified! This is unlikely to do what you want.')
  csts[cost.m>0,gmsc:=cost.sd^2/cost.m]
  csts[!is.na(gmsc) & gmsc > 0, gmk:=cost.m/gmsc]
  NR <- nrow(csts)
  csts <- csts[rep(1:NR,nrep)]
  csts[,id:=rep(1:nrep,each=NR)]
  csts[,rnd:=!is.na(gmsc) & !is.na(gmk) & gmk>0 & gmsc > 0]
  csts[rnd==TRUE,value:=rgamma(sum(rnd),shape=gmk,scale = gmsc)] #random sample from gamma distribution
  csts[rnd!=TRUE,value:=cost.m]                                  #fixed values
  ## csts[,cnms:=paste0('c_',cost)]
  csts[,cnms:=paste0(cost)]
  F <- 'id '
  if(!is.null(anmz)) F <- paste0(F,'+ ',paste(anmz,collapse='+')) #split out by attributes if included
  F <- paste0(F, ' ~ cnms')
  dcast(csts,as.formula(F),value.var = 'value')      #id ~ cnms
}
## NOTE
## if attributes are included, all costs need to be specified by them even if this means duplicating those without dependence


## making life years
GetLifeYears <- function(isolist,discount.rate,yearfrom){
  ## template:
  LYT <- data.table(age=15:100,
                    age_group=c(rep('15-64',50),rep('65-100',36)),
                    LYS=0.0)
  ## make country/age key
  LYK <- list()
  for(iso in isolist){
    ## iso <- cn
    tmp <- copy(LYT)
    tmp[,iso3:=iso]
    for(ag in tmp$age)
      tmp[age==ag,LYS:=discly::discly(iso3=iso,
                                      age=ag,
                                      yearnow=yearfrom,
                                      sex='Total',
                                      endyear = 2098,
                                      HR=1,
                                      dr=discount.rate,
                                      hiv='both'
      )]
    LYK[[iso]] <- tmp
  }
  LYK <- rbindlist(LYK)
  ## assume unweighted & collapse
  LYK <- LYK[,.(LYS=mean(LYS)),by=.(iso3,age=age_group)]
  setkey(LYK,age)
  LYK
}

## ## scraps for development of below fn
## data <- merge(D,LYK,by='age') #add age
## Kmax <- 1e3
## file.id <- 'test'
## wtp <- 500

## NOTE this is more illustrative for now
## NOTE needs a folder called graphs/ creating (which is currently excluded from the repo)
## some automatic CEA outputs
file.id='';Kmax=5e3;wtp=5e3;
MakeCEAoutputs <- function(data,LY,
                           file.id='',Kmax=5e3,wtp=5e3,
                           arms=c('SOC','INT')){
  data <- merge(data,LY,by='age') #add age
  DS <- D[,.(cost.SOC=sum(cost.soc*value),
                cost.INT=sum(cost.int*value),
                lyl.SOC=sum(deaths.soc*value*LYS),
                lyl.INT=sum(deaths.int*value*LYS)),
             by=id] #PSA summary

  ## prep for BCEA
  LYS <- CST <- matrix(nrow=nreps,ncol=2)
  LYS[,1] <- 1-DS$lyl.SOC #NOTE this is life years lost
  LYS[,2] <- 1-DS$lyl.INT
  CST[,1] <- DS$cost.SOC
  CST[,2] <- DS$cost.INT
  ## BCEA outputs
  M <- bcea(e=LYS,c=CST,ref=1,interventions = arms,Kmax=Kmax)
  print(summary(M))

  fn <- paste0(here('outdata/kstar_'),file.id,'.txt')
  cat(M$kstar,file = fn)
  fn <- paste0(here('outdata/ICER_'),file.id,'.txt')
  cat(M$ICER,file = fn)

  ## NOTE may need more configuration
  ceac.plot(M,graph='ggplot2') +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/CEAC_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  ceplane.plot(M,graph='ggplot2',wtp=wtp)+
    scale_x_continuous(label=comma) +
    theme_classic() +
    theme(legend.position = 'top') + ggpubr::grids()
  fn <- paste0(here('graphs/CE_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  eib.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EIB_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  evi.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EVI_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

}

## --- for reformatting costs
reformatCosts <- function(rcsts){
  iextra <- outer(isoz,c('.lo','.hi','drop'),paste0)
  iextra <- c(t(iextra)); iextra <- rev(rev(iextra)[-1])
  nnmz <- c('drop','DESCRIPTION','NAME',iextra)
  names(rcsts)[1:length(nnmz)] <- nnmz
  drop <- grep('drop',names(rcsts),value=TRUE)
  rcsts[,c(drop):=NULL]
  rcsts[is.na(rcsts)] <- 0.1 #dummy
  rcsts <- melt(rcsts,id=c('NAME','DESCRIPTION'))
  rcsts[,DESCRIPTION:=NULL]
  rcsts[,c('iso3','hilo'):=tstrsplit(variable,split="\\.")]
  rcsts <- dcast(rcsts,iso3 + NAME ~ hilo,value.var = 'value')
  rcsts[,c('cost.m','cost.sd'):=.((lo+hi)/2,(hi-lo)/3.92)]
  rcsts <- rcsts[,.(iso3,cost=NAME,cost.m,cost.sd)]
  rcsts
}

## calculate mid/lo/hi
MLH <- function(dat){
  nnmz <- names(dat)
  lnmz <- paste0(nnmz,'.lo')
  hnmz <- paste0(nnmz,'.hi')
  mnmz <- paste0(nnmz,'.mid')
  L <- dat[,lapply(.SD,lo),.SDcols=nnmz]
  M <- dat[,lapply(.SD,mean),.SDcols=nnmz]
  H <- dat[,lapply(.SD,hi),.SDcols=nnmz]
  setnames(L,nnmz,lnmz); setnames(M,nnmz,mnmz); setnames(H,nnmz,hnmz);
  list(L=L,M=M,H=H)
}
## MLH(out[,.(DcostperATT.int,DcostperATT.soc)]) #test


## =========== output formatters
outsummary <- function(out){

  ## mid/lo/hi
  outa <- MLH(out[,.(costperATT.soc,costperATT.int,
                     DcostperATT.int,
                     costperTPT.soc,costperTPT.int,
                     DcostperTPT.int,
                     # Ddeaths.int,
                     # DLYL.int,
                     # DLYL0.int,
                     # DcostperLYS0.int,
                     # DcostperLYS.int,
                     # Dcostperdeaths.int,
                     Dcost.int)])

  ## more bespoke statistics
  # outi <- out[,.(ICER.int= -mean(Dcost.int) / mean(DLYL.int))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H)) #combine

  ## pretty version
  pouts <- outs[,.(costperATT.soc = brkt(costperATT.soc.mid,costperATT.soc.lo,costperATT.soc.hi),
                   costperATT.int = brkt(costperATT.int.mid,costperATT.int.lo,costperATT.int.hi),
                   DcostperATT.int = brkt(DcostperATT.int.mid,DcostperATT.int.lo,DcostperATT.int.hi), # TODO::quick gap measure check!!
                   costperTPT.soc = brkt(costperTPT.soc.mid,costperTPT.soc.lo,costperTPT.soc.hi),
                   costperTPT.int = brkt(costperTPT.int.mid,costperTPT.int.lo,costperTPT.int.hi),
                   DcostperTPT.int = brkt(DcostperTPT.int.mid,DcostperTPT.int.lo,DcostperTPT.int.hi)
                  )]

  ## return value
  list(outs=outs,pouts=pouts)
}


## ---- utilities for making CEACs
make.ceac <- function(CEA,lamz){
    crv <- lamz
    for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
    crv
}


## additional table 2
## --- cols:
## CMR, UGA x SOC, INT
## --- rows:
## contacts = value //
## TPT courses = tpt //
## ATT courses = att //
## prev TB = coprevtb
## inc TB = inctb
## prev deaths = deaths-incdeaths
## inc tb deaths = incdeaths
## discounted LYL TODO //
## ATT cost TODO
## TPT cost TODO
## total cost TODO //
## ICER TODO



## =========== output formatters
Table2 <- function(dat){
  
  ## mid/lo/hi
  outa <- MLH(dat[,.(
    Dscreen,screen.int,screen.soc,
    Dtpt,tpt.int,tpt.soc,
    Datt,att.int,att.soc,
    # DLYL,LYL.int,LYL.soc,
    Dcoprevtb,coprevtb.int,coprevtb.soc,
    # Dinctb,inctb.int,inctb.soc,
    # Dincdeaths,incdeaths.int,incdeaths.soc,
    # Dprevdeaths,prevdeaths.int,prevdeaths.soc,
    # Ddeaths,deaths.int,deaths.soc,
    Dcost.screen,cost.screen.int,cost.screen.soc,
    Dcost.tpt,cost.tpt.int,cost.tpt.soc,
    Dcost.att,cost.att.int,cost.att.soc,
    # Dcost.inc.att,cost.inc.att.int,cost.inc.att.soc,
    Dcost,cost.int,cost.soc
  )])
  
  ## more bespoke statistics
  # outi <- dat[,.(ICER.int= -mean(Dcost) / mean(DLYL))]
  
  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H)) #combine
  
  ## pretty version
  fac <- 1e3 # per fac new residents
  pouts <- outs[,.(
    Dscreen = brkt(fac*Dscreen.mid,fac*Dscreen.lo,fac*Dscreen.hi),
    screen.int = brkt(fac*screen.int.mid,fac*screen.int.lo,fac*screen.int.hi),
    screen.soc = brkt(fac*screen.soc.mid,fac*screen.soc.lo,fac*screen.soc.hi),
    Dtpt = brkt(fac*Dtpt.mid,fac*Dtpt.lo,fac*Dtpt.hi),
    tpt.int = brkt(fac*tpt.int.mid,fac*tpt.int.lo,fac*tpt.int.hi),
    tpt.soc = brkt(fac*tpt.soc.mid,fac*tpt.soc.lo,fac*tpt.soc.hi),
    Datt = brkt(fac*Datt.mid,fac*Datt.lo,fac*Datt.hi),
    att.int = brkt(fac*att.int.mid,fac*att.int.lo,fac*att.int.hi),
    att.soc = brkt(fac*att.soc.mid,fac*att.soc.lo,fac*att.soc.hi),
    Dcoprevtb = brkt(fac*Dcoprevtb.mid,fac*Dcoprevtb.lo,fac*Dcoprevtb.hi),
    coprevtb.int = brkt(fac*coprevtb.int.mid,fac*coprevtb.int.lo,fac*coprevtb.int.hi),
    coprevtb.soc = brkt(fac*coprevtb.soc.mid,fac*coprevtb.soc.lo,fac*coprevtb.soc.hi),
    Dcost.screen = brkt(fac*Dcost.screen.mid,fac*Dcost.screen.lo,fac*Dcost.screen.hi),
    cost.screen.int = brkt(fac*cost.screen.int.mid,fac*cost.screen.int.lo,fac*cost.screen.int.hi),
    cost.screen.soc = brkt(fac*cost.screen.soc.mid,fac*cost.screen.soc.lo,fac*cost.screen.soc.hi),
    Dcost.tpt = brkt(fac*Dcost.tpt.mid,fac*Dcost.tpt.lo,fac*Dcost.tpt.hi),
    cost.tpt.int = brkt(fac*cost.tpt.int.mid,fac*cost.tpt.int.lo,fac*cost.tpt.int.hi),
    cost.tpt.soc = brkt(fac*cost.tpt.soc.mid,fac*cost.tpt.soc.lo,fac*cost.tpt.soc.hi),
    Dcost.att = brkt(fac*Dcost.att.mid,fac*Dcost.att.lo,fac*Dcost.att.hi),
    cost.att.int = brkt(fac*cost.att.int.mid,fac*cost.att.int.lo,fac*cost.att.int.hi),
    cost.att.soc = brkt(fac*cost.att.soc.mid,fac*cost.att.soc.lo,fac*cost.att.soc.hi),
    Dcost = brkt(fac*Dcost.mid,fac*Dcost.lo,fac*Dcost.hi),
    cost.int = brkt(fac*cost.int.mid,fac*cost.int.lo,fac*cost.int.hi),
    cost.soc = brkt(fac*cost.soc.mid,fac*cost.soc.lo,fac*cost.soc.hi)
    # ICER.int=round(ICER.int,0)
  )]
  
  ## return value
  list(outs=outs,pouts=pouts)
}

