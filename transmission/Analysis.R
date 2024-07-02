library(here) #NOTE this uses petedodd/ecrins: keep versions in sync!
source(here('transmission/utilities.R'))
set.seed(1234)

## example output
tt <- seq(from=0, to=10, by=0.1)  #time frame to run over
y <- runmodel(tt,parms)           #run model
head(y)

## example use of using 1st row of churn samples
parms <- revise.flow.parms(parms,smpsd,1) #change flow parameter to 1st row of samples
y <- runmodel(tt,parms)           #run model
tail(y)

## trying with tree parameters also
parms <- revise.HE.parms(parms,DR,1,arm='soc')
y <- runmodel(tt,parms)           #run model
tail(y)

parms <- revise.HE.parms(parms,DR,1,arm='int')
y <- runmodel(tt,parms)           #run model
tail(y)


## run test NOTE breaks
parms$staticfoi <- -1
tt <- seq(from=0, to=50, by=0.1)  #time frame to run over
y <- runmodel(tt,parms)           #run model
ydm <- output2dt(y)
## ydm[!is.na(population),unique(variable)] #all TB states

## PPD look:
pop <- ydm[!is.na(population),.(value=sum(value)),by=.(t,population)]
popt <- pop[population!='previously detained',.(value=sum(value)),by=t] #total
ggplot(popt,aes(t,value))+geom_line()+ylim(c(0,NA))

ggplot(pop[population!='previously detained'],aes(t,value,col=population))+
  geom_line()+ylim(c(0,NA))

## TB look:
tbpop <- ydm[!is.na(population) & population!='previously detained',
             .(value=sum(value)),
             by=.(t,variable)]
ggplot(tbpop,aes(t,value,col=variable))+
  geom_line()+facet_wrap(~variable,scales='free')

tbpopf <- ydm[!is.na(population) & population!='previously detained',
             .(value=sum(value)),
             by=.(t,variable,population)]
ggplot(tbpopf,aes(t,value,col=population))+
  geom_line()+facet_wrap(~variable,scales='free')

## TPT look:
tptpop <- ydm[!is.na(population) & population!='previously detained',
             .(value=sum(value)),
             by=.(t,tpt,population)]
ggplot(tptpop,aes(t,value,col=tpt))+
  geom_line()+facet_wrap(~population,scales='free')


## Other outputs:
ggplot(ydm[variable %in% others],aes(t,value,col=variable))+
  geom_line()+facet_wrap(~variable,scales='free')

## ============= diff plot  =============

parms$staticfoi <- -1            #dynamic
parms$int_time <- 50             #fix
ydb <- diffdata(parms)

ggplot(ydb[variable %in% others],aes(t,value,col=variable,lty=arm))+
  geom_line()+facet_wrap(~variable,scales='free')

ggsave(here('transmission/plots/p_illustrate_dyn.png'),w=7,h=5)

## ============= HE workflow =============

## NOTE LTBI & foi parameters
parms$parm_frac_L <- 0.2
summary(qfun(runif(1e3),hyperparms$foi))

parms <- revise.flow.parms(parms,smpsd,1) #use some real flow parms

## RUN TEST PSA:
RES <- PSAloop(Niter=2e2,parms,smpsd,DR,zero.nonscreen.costs=TRUE) #TODO update this as we go

## ============= CHECKING ===================
DRS[quantity=='cost' & outcome=='notx'] #pretty similar int/soc
DRS[quantity=='cost' & outcome=='att'] #pretty similar int/soc
DRS[quantity=='cost' & outcome=='tpt'] #pretty similar int/soc


## incoming state?
## inflow split: parms from doing library(ecrins); data(parms)
IFS <- with(data=parms,{
  c(TBD=parm_frac_SD+parm_frac_CD,
    TBI=parm_frac_E+parm_frac_L+parm_frac_epTB+parm_frac_lpTB,
    noTB=parm_frac_U
  )
})

## transition matrices inflow (cols) to outcome (rows)
psoc <- as.matrix(DRS[arm=='soc' & quantity=='check',.(TBD,TBI,noTB)]) #rows att,notx,tpt
pint <- as.matrix(DRS[arm=='int' & quantity=='check',.(TBD,TBI,noTB)]) #rows att,notx,tpt
csoc <- as.matrix(DRS[arm=='soc' & quantity=='cost',.(TBD,TBI,noTB)]) #rows att,notx,tpt
cint <- as.matrix(DRS[arm=='int' & quantity=='cost',.(TBD,TBI,noTB)]) #rows att,notx,tpt
cint[!is.finite(cint)] <- 0 ; csoc[!is.finite(csoc)] <- 0


## average outcomes: almost the same! (transition matrix)
(avosoc <- psoc %*% IFS)
(avoint <- pint %*% IFS)

## average costs by outcome:
(avcsoc <- (psoc*csoc) %*% IFS)
(avcint <- (pint*cint) %*% IFS)

## average cost per person:
sum(avcsoc)
sum(avcint)

## unit cost by outcome
avcsoc/avosoc
avcint/avoint

## TODO check ATT$ for TBI > ATT$ TBD

DRS[quantity=='check' & outcome=='notx',] #pretty similar int/soc

## --- compare ODE HE results
DRtmp <- DR[1:3]
tmp <- melt(DRS,id=c('arm','quantity','outcome'))
tmp <- dcast(tmp,variable ~ arm + outcome + quantity,value.var = 'value')
names(tmp)[1] <- 'tb'
tmp[,id:=1]

## get results for single run using DRS summary parms
res1 <- run.HE.socint(parms,tmp,1,zero.nonscreen.costs=TRUE)
res1[,soc.CC0/(int.to.end*parms$inflow)] #same as below 892
res1[,soc.CC/(int.to.end*parms$inflow)] #discounted ~ 333

res2 <- run.HE.socint(parms,tmp,1,zero.nonscreen.costs=FALSE)
res2[,soc.CC0/(int.to.end*parms$inflow)] # ~20% more than screening only

## ensure using same parms
parms1 <- revise.HE.parms(parms,tmp,1,arm='soc',zero.nonscreen.costs=TRUE)

## from ecrins/tbmod0.R
## cost rate per inflow used in ODEs, SOC
with(parms1,{
  1 * (
    ## TPT
    uc_entry_tpt_TB * (inflow_toTPT_TB0) * (parm_frac_SD+parm_frac_CD) +
    uc_entry_tpt_L * (inflow_toTPT_L0) * (parm_frac_E+parm_frac_L + parm_frac_epTB+parm_frac_lpTB) +
    uc_entry_tpt_no * (inflow_toTPT_no0) * parm_frac_U +
    ## ATT
    uc_entry_att_TB * (inflow_toATT_TB0) * (parm_frac_SD + parm_frac_CD) +
    uc_entry_att_L * (inflow_toATT_L0) * (parm_frac_E+parm_frac_L + parm_frac_epTB+parm_frac_lpTB) + #FP for L
    uc_entry_att_no * (inflow_toATT_no0) * parm_frac_U + #FP for no
    ## no TX
    uc_entry_notx_TB * (1-inflow_toTPT_TB0-inflow_toATT_TB0) * (parm_frac_SD+parm_frac_CD) +
    uc_entry_notx_L * (1-inflow_toTPT_L0-inflow_toATT_L0) * (parm_frac_E+parm_frac_L + parm_frac_epTB+parm_frac_lpTB) +
    uc_entry_notx_no * (1-inflow_toTPT_no0-inflow_toATT_no0) * parm_frac_U
  )
}) #892

psoc1 <- as.matrix(DRS[arm=='soc' & quantity=='check',.(TBD,TBI,noTB)]) #rows att,notx,tpt
csoc1 <- as.matrix(DRS[arm=='soc' & quantity=='cost',.(TBD,TBI,noTB)]) #rows att,notx,tpt
csoc1[!is.finite(csoc1)] <- 0
(avcsoc1 <- (psoc1*csoc1) %*% IFS)
sum(avcsoc1) #892

## NOTE OK perfectly the same
## ---


## TODO impact checks
res1[,soc.CC0/soc.CC] #2.6

res1[,soc.ccasesout/soc.ccases] #over half cases out?
res1[,soc.dLYL/soc.deaths]   #TODO check against discounted LE

## TODO double check
((1-exp(-parms$LifeExp *parms$disc_rate))/parms$disc_rate)*
  (1-exp(-parms$disc_rate*res1$int.to.end))/parms$disc_rate

res1[,soc.deaths/soc.ccases] #BUG CFR too high
res1[,soc.cATTtp/soc.ccases] #86% OK as CDR

## TODO parameters with all on TPT and TPT effect absolute?
## 

## ====================

## RUN REAL PSA:
RES <- PSAloop(Niter=2e3,parms,smpsd,DR,zero.nonscreen.costs=FALSE) #TODO update this as we go


## =========== OUTPUTS & PLOTS ==================

ggplot(RES,aes(Q.soc - Q.int,int.CC-soc.CC,col=mid.notes<100 & mid.notes>30))+
  geom_point(shape=1) +
  scale_y_continuous(label=comma,'Incremental cost',limits=c(0,NA))+
  scale_x_continuous(label=comma,'Health gain',limits=c(0,NA))+
  geom_abline(intercept=0,slope=30e3,col=2)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position='top')

ggsave(file=here('transmission/plots/p_CE1.png'),w=7,h=7)

RES[,mean(int.CC-soc.CC)]
RES[,mean(Q.int - Q.soc)]
RES[,mean(int.CC-soc.CC)/mean(Q.int - Q.soc)]

RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)]
RES[mid.notes<100 & mid.notes>30,mean(Q.int - Q.soc)] #BUG?
RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)/mean(Q.int - Q.soc)]

## make CEAC
ll <- seq(from=0,to=50e3,by=100)
cc <- make.ceac(RES[,.(Q=dQ,P=int.CC-soc.CC)],ll)
ccr <- make.ceac(RES[,## mid.notes<100 & mid.notes>30,
                     .(Q=dQ,P=int.CC-soc.CC)],
                 ll)
D <- data.table(x=ll,all=cc,restricted=ccr)
D <- melt(D,id='x')

## plot
xpad <- 2.5e3
GP <- ggplot(D,aes(x,value,col=variable)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top')+
  scale_y_continuous(label=percent,
                     limits=c(0,1),
                     expand = c(0,0))+
  scale_x_continuous(label=comma,
                     limits=c(0,max(D$x)),
                     expand = c(0,xpad))+
  ## scale_color_colorblind()+
  scale_colour_manual(values=cbPalette)+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  ggpubr::grids()
GP


ggsave(file=here('transmission/plots/p_CEAC1.png'),w=7,h=7)
