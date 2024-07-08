library(here) #NOTE this uses petedodd/ecrins: keep versions in sync!
source(here('transmission/utilities.R'))

sd <- as.integer(1234)
set.seed(sd)

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
parms <- revise.flow.parms(parms,smpsd,1) #use some real flow parms




## testing
## summary(qfun(runif(1e3),hyperparms$foi))
## revise.HE.parms(parms,DR,1)
## revise.HE.parms(parms,DR,1,summarize.HEparms = TRUE)

## ans <- run.HE.socint(parms,DR,1)


## ## RUN TEST PSA:
## RES <- PSAloop(Niter=1e2,parms,smpsd,DR,zero.nonscreen.costs=TRUE)

## ============= CHECKING ===================
## see checks.R

## ====================
## TODO adjust hyperparms for foi and disease to match data
Nruns <- 1e3
parms$staticfoi <- -1            #dynamic=-1
RES <- PSAloop(Niter = Nruns, parms, smpsd, DR, zero.nonscreen.costs = TRUE, verbose = FALSE)
summary(RES$problem)

RES[,problem:=NULL]

## ===== inspect
RES[,mean(int.CC-soc.CC)]
RES[,mean(int.deaths-soc.deaths)] #fewer deaths
RES[, mean(int.ccases - soc.ccases)] # fewer cases
RES[, mean(int.dLYL - soc.dLYL)] # fewer LYL
RES[, mean(Q.int - Q.soc)] # less qol decrement + LYL
RES[,mean(int.CC-soc.CC)/mean(dQ)]/1e3 #ICER ~ 777K; static + w/o multiplier + w/o PTB = 1M

RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)]
RES[mid.notes<100 & mid.notes>30,mean(Q.int - Q.soc)]
RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)/mean(dQ)]/1e3 # 1M




## ======== table
tabkey <- c(
  "Undiscounted Costs" = "CC0",
  "QALYs" = "Q",
  "ATT courses" = "cATT",
  "TPT courses" = "cTPT",
  "Incident TB" = "ccases",
  "Life-years lost to TB (discounted)" = "dLYL",
  "TB deaths" = "deaths",
  "QoL lost to TB (discounted)" = "qoldec"
)
tabkey <- data.table(name=names(tabkey),quantity=tabkey)
tabout <- RES[, .( # entries!
  entries = inflow * 70,
  ## TPT courses
  soc.cTPT, int.cTPT, inc.cTPT = int.cTPT - soc.cTPT,
  ## ATT courses
  soc.cATT = soc.cATTtp + 0, int.cATT = int.cATTtp + 0, inc.cATT = int.cATTtp - soc.cATTtp,
  ## Costs
  soc.CC0, int.CC0, inc.CC0 = int.CC0 - soc.CC0,
  ## Incident TB
  soc.ccases, int.ccases, inc.ccases = int.ccases - soc.ccases,
  ## TB deaths
  soc.deaths, int.deaths, inc.deaths = int.deaths - soc.deaths,
  ## TB LYL
  soc.dLYL, int.dLYL, inc.dLYL = int.dLYL - soc.dLYL,
  ## TB QoL loss
  soc.qoldec, int.qoldec, inc.qoldec = int.qoldec - soc.qoldec,
  ## TB QALYs
  soc.Q = Q.soc, int.Q = Q.int, inc.Q = Q.int - Q.soc
)]
tabout <- tabout[,lapply(.SD,function(x)1e4*x/entries)] #per 10K entries
tabout[, entries := NULL]
tabout[,id:=1:nrow(tabout)]
tabout <- melt(tabout, id = "id")
tabout[, c("arm", "quantity") := tstrsplit(variable, split = "\\.")]
tabout <- tabout[,.(mid=mean(value),lo=lo(value),hi=hi(value)),by=.(arm,quantity)] #
tabout[,txt:=brkt(mid,lo,hi)]
## TODO may need tweaking for costs and for format
tabout <- dcast(tabout, quantity ~ arm, value.var = "txt")
tabout <- merge(tabout, tabkey, by = "quantity")
tabout <- tabout <- tabout[
  c("cTPT", "cATT", "CC0", "ccases", "deaths", "dLYL", "qoldec", "Q"),
  .(name, soc, int, inc)
] # reorder

fwrite(tabout,file=here('transmission/plots/tabout.csv'))


## =========== OUTPUTS & PLOTS ==================

ggplot(RES,aes(Q.soc - Q.int,int.CC-soc.CC,col=mid.notes<100 & mid.notes>30))+
  geom_point(shape=1) +
  scale_y_continuous(label=comma,'Incremental cost',limits=c(0,NA))+
  scale_x_continuous(label=comma,'Health gain',limits=c(0,NA))+
  geom_abline(intercept=0,slope=30e3,col=2)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position='top')

ggsave(file=here('transmission/plots/p_CE1.png'),w=7,h=7)

## make CEAC
ll <- seq(from=0,to=1e5,by=100)
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


## ------------ patterns wrt ICER & TBI
## ICER by prev & TBI
RES[, c("TBD", "TBI") := .(
  parm_frac_SD + parm_frac_CD,
  parm_frac_E + parm_frac_L + parm_frac_epTB + parm_frac_lpTB
  )]
## categories with equal numbers
RES[, TBDC := cut(TBD, quantile(TBD,probs=seq(0,1,l=10)),include.lowest = TRUE)]
RES[, TBIC := cut(TBI, quantile(TBI, probs = seq(0, 1, l = 10)), include.lowest = TRUE)]

## RES[,range(TBD)]
## RES[,unique(TBDC)]
## RES[, table(TBDC)]

smyd <- RES[, .(ICER = mean(int.CC - soc.CC) / mean(dQ)), by = TBDC]
smyi <- RES[, .(ICER = mean(int.CC - soc.CC) / mean(dQ)), by = TBIC]

str(smyd)

ggplot(smyd, aes(TBDC, ICER)) +
  geom_point() +
  expand_limits(y = 0)

ggplot(smyi, aes(TBIC, ICER)) +
  geom_point() +
  expand_limits(y = 0)




## looking at SAVI
save(RES,file='~/Dropbox/Holocron/tmp/RES.Rdata')
load('~/Dropbox/Holocron/tmp/RES.Rdata')


RES

## ============ VOI stuff =========
## effects
mx <- max(RES$Q.int,RES$Q.soc)
EFF <- RES[,.(effect1=mx-Q.soc,effect2=mx-Q.int)]
fwrite(EFF,file='~/Downloads/EFF.csv')

## costs
CST <- RES[,.(costs1=soc.CC,costs2=int.CC)]
fwrite(CST,file='~/Downloads/CST.csv')

## parameters
PMS <- RES[,.(
  parm_frac_U, parm_frac_E,
  parm_frac_L, parm_frac_SD, parm_frac_CD, parm_frac_ATT,
  parm_frac_epTB, parm_frac_lpTB, parm_ifrac_U, parm_ifrac_E,
  parm_ifrac_L, parm_ifrac_SD, parm_ifrac_CD, parm_ifrac_ATT, 
  parm_ifrac_epTB, parm_ifrac_lpTB, parm_ifrac_prevTPT1, parm_ifrac_prevTPT2, 
  parm_ifrac_prevTPT3, inflow, remand_short, remand_long, 
  remand_release, long_short, short_release, short_open, 
  long_release, open_release, previous_remand, parm_init_PPD1, 
  parm_init_PPD2, parm_init_PPD3, parm_init_PPD4, parm_init_PPD5, 
  staticfoi, ptn, foi, stb, 
  prg, eps, rel, CDR, 
  drn, txf, CFR, m, 
  tptHR, tpt_drn, wsn, mHR, 
  att_time, late_post_time, mort, hrqolptb, 
  int_time, disc_rate, LifeExp, uc_attppd, 
  uc_attout, hrqol, inflow_toATT_TB0, inflow_toATT_L0, 
  inflow_toATT_no0, inflow_toTPT_TB0, inflow_toTPT_L0, inflow_toTPT_no0, 
  inflow_toATT_TB1.soc, inflow_toATT_L1.soc, inflow_toATT_no1, inflow_toTPT_TB1.soc, 
  inflow_toTPT_L1.soc, inflow_toTPT_no1.soc, uc_entry_tpt_TB.soc, uc_entry_tpt_L.soc, 
  uc_entry_tpt_no.soc, uc_entry_att_TB.soc, uc_entry_att_L.soc, uc_entry_att_no.soc, 
  uc_entry_notx_TB.soc, uc_entry_notx_L.soc, uc_entry_notx_no.soc, inflow_toATT_TB1.int, 
  inflow_toATT_L1.int, inflow_toATT_no1.int, inflow_toTPT_TB1.int, inflow_toTPT_L1.int, 
  inflow_toTPT_no1.int, uc_entry_tpt_TB.int, uc_entry_tpt_L.int, uc_entry_tpt_no.int, 
  uc_entry_att_TB.int, uc_entry_att_L.int, uc_entry_att_no.int, uc_entry_notx_TB.int, 
  uc_entry_notx_L.int, uc_entry_notx_no.int
)]
fwrite(PMS,file='~/Downloads/PMS.csv')

 



## =========== loop for SAs
## TODO

set.seed(sd)

SA <- list()
## bascase results from above
SA[[1]] <- data.table(
  analysis = "base case",
  ICER = RES[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
## Static model
set.seed(sd)
RES1 <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = TRUE, community = TRUE, posttb = TRUE
)
SA[[2]] <- data.table(
  analysis = "Static model",
  ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES1[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
## Static model, no community transmission
set.seed(sd)
RES1 <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = TRUE, community = FALSE, posttb = TRUE
)
SA[[3]] <- data.table(
  analysis = "Static model, no community transmission",
  ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES1[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
## Static model, no community transmission or post-TB effects
set.seed(sd)
RES1 <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = TRUE, community = FALSE, posttb = FALSE
)
SA[[4]] <- data.table(
  analysis = "Static model, no community transmission or post-TB effects",
  ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES1[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
## join
SA <- rbindlist(SA)
SA[, `ICER, unrestricted TB rates` := paste0(round(ICER))]
SA[, `ICER, restricted TB rates` := paste0(round(ICERr))]
SA[, c("ICER", "ICERr") := NULL]

fwrite(SA,file=here('transmission/plots/SA.csv'))
