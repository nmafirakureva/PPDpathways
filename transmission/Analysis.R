library(here) #NOTE this uses petedodd/ecrins: keep versions in sync!
source(here('transmission/utilities.R'))

sd <- as.integer(1234)
set.seed(sd)

## example output
tt <- seq(from=0, to=10, by=0.1)  #time frame to run over
y <- runmodel(tt,parms)           #run model
head(y)

## ====================
Nruns <- 1e3

## parms$staticfoi <- -1 # dynamic=-1
set.seed(sd)
RES <- PSAloop(
  Niter = 5*Nruns, parms, smpsd, DR,
  zero.nonscreen.costs = FALSE, verbose = FALSE, targeting = FALSE
)
summary(RES$problem)
RES[,problem:=NULL]

## ===== inspect
RES[,mean(int.CC-soc.CC)]
RES[,mean(int.deaths-soc.deaths)] #fewer deaths
RES[, mean(int.ccases - soc.ccases)] # fewer cases
RES[, mean(int.dLYL - soc.dLYL)] # fewer LYL
RES[, mean(Q.int - Q.soc)] # less qol decrement + LYL
RES[,mean(int.CC-soc.CC)/mean(dQ)]/1e3 #54

RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)]
RES[mid.notes<100 & mid.notes>30,mean(Q.int - Q.soc)]
RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)/mean(dQ)]/1e3 # 1M


## ======== table
tabkey <- c(
  "Undiscounted Costs" = "CC0",
  "QALYs lost to TB" = "Q",
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
  soc.cATT = soc.cATTtp + soc.cATTfp, int.cATT = int.cATTtp + int.cATTfp,
  inc.cATT = int.cATTtp + int.cATTfp - soc.cATTtp - soc.cATTfp,
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
tabout[, c("soc.CC0", "int.CC0", "inc.CC0") := .(soc.CC0 / 1e2, int.CC0 / 1e2, inc.CC0 / 1e2)] # NOTE
tabout[, entries := NULL]
tabout[,id:=1:nrow(tabout)]
tabout <- melt(tabout, id = "id")
tabout[, c("arm", "quantity") := tstrsplit(variable, split = "\\.")]
tabout <- tabout[,.(mid=mean(value),lo=lo(value),hi=hi(value)),by=.(arm,quantity)] #
tabout[,txt:=brkt(mid,lo,hi)]
## formatting
tabout <- dcast(tabout, quantity ~ arm, value.var = "txt")
tabout <- merge(tabout, tabkey, by = "quantity")
tabout <- tabout[
  c("cTPT", "cATT", "CC0", "ccases", "deaths", "dLYL", "qoldec", "Q"),
  .(name, soc, int, inc)
] # reorder
print(tabout)

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
ccr <- make.ceac(RES[mid.notes<100 & mid.notes>30,
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
  xlab('Cost-effectiveness threshold (GBP/QALY gained)')+
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

ggplot(smyi, aes(TBIC, ICER,group=1)) +
  geom_point() + geom_line()+
  scale_y_continuous(limits =c(0,NA),label=comma)+
  xlab('TB infection prevalence')+
  ylab('ICER (GBP per QALY gained)')+
  theme_linedraw()

ggsave(file = here("transmission/plots/p_ICERbyTBI.png"), w = 12, h = 7)


## looking at SAVI
fn <- here("transmission/data/RES.Rdata")
save(RES, file = fn)

load(fn)



## =========== loop for SAs
RESA <- SA <- list()
## bascase results from above
SA[[1]] <- data.table(
  analysis = "base case",
  ICER = RES[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
RESA[[1]] <- copy(RES)
RESA[[1]][, analysis := "base case"]
## Static model
set.seed(sd)
RES1 <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = TRUE, community = TRUE, posttb = TRUE, targeting = FALSE
)
SA[[2]] <- data.table(
  analysis = "Static model",
  ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES1[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
RESA[[2]] <- copy(RES1)
RESA[[2]][, analysis := "Static model"]
## Static model, no community transmission
set.seed(sd)
RES1 <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = TRUE, community = FALSE, posttb = TRUE, targeting = FALSE
)
SA[[3]] <- data.table(
  analysis = "Static model, no community transmission",
  ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES1[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
RESA[[3]] <- copy(RES1)
RESA[[3]][, analysis := "Static model, no community transmission"]
## Static model, no community transmission or post-TB effects
set.seed(sd)
RES1 <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = TRUE, community = FALSE, posttb = FALSE, targeting = FALSE
)
SA[[4]] <- data.table(
  analysis = "Static model, no community transmission or post-TB effects",
  ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
  ICERr = RES1[
    mid.notes < 100 & mid.notes > 30,
    mean(int.CC - soc.CC) / mean(dQ)
  ]
)
RESA[[4]] <- copy(RES1)
RESA[[4]][, analysis := "Static model, no community transmission or post-TB effects"]

## join
SAT <- rbindlist(SA)
SAT[, `ICER, unrestricted TB rates` := paste0(round(ICER))]
SAT[, `ICER, restricted TB rates` := paste0(round(ICERr))]
SAT[, c("ICER", "ICERr") := NULL]
print(SAT)

fwrite(SAT, file = here("transmission/plots/SA.csv"))


## extras from DR files
for (snm in saznmz) {
  print(snm)
  set.seed(sd)
  RES1 <- PSAloop(
    Niter = Nruns, parms, smpsd, Dlist[[snm]],
    zero.nonscreen.costs = FALSE, verbose = FALSE, targeting = FALSE
  )
  SA[[snm]] <- data.table(
    analysis = snm,
    ICER = RES1[, mean(int.CC - soc.CC) / mean(dQ)],
    ICERr = RES1[
      mid.notes < 100 & mid.notes > 30,
      mean(int.CC - soc.CC) / mean(dQ)
    ]
  )
  RESA[[snm]] <- copy(RES1)
  RESA[[snm]][, analysis := snm]
}

## join
SATF <- rbindlist(SA)
SATF[, `ICER, unrestricted TB rates` := paste0(round(ICER))]
SATF[, `ICER, restricted TB rates` := paste0(round(ICERr))]
SAT[, c("ICER", "ICERr") := NULL]
print(SATF)

fwrite(SATF, file = here("transmission/plots/SA.full.csv"))


## trim 1st element
dropcolz <- setdiff(names(RESA[[1]]), names(RESA[[2]]))
RESA[[1]][, c(dropcolz) := NULL]
for (i in 2:length(RESA)) RESA[[i]][, problem := NULL]

## join & save
RESA <- rbindlist(RESA)

fn <- here("transmission/data/RESA.Rdata")
save(RESA, file = fn)


## ========== TARGETING
RES <- PSAloop(
  Niter = 10e3, parms, smpsd, DR,
  zero.nonscreen.costs = FALSE, verbose = FALSE, targeting = TRUE
)
summary(RES$problem)
RES[,problem:=NULL]

RES[, OR := odds.ratio(SE1, SP1)]
RES[, frac.screened := prevTBI * SE1 + (1 - prevTBI) * (1 - SP1)]
RES[, OR.cat := cut(OR, breaks = c(0, 1, 2, 4, 10), include.lowest = TRUE)]
RES[, x.cat := cut(frac.screened, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)]
RES[, NB := 30e3 * dQ - (int.CC - soc.CC)]

## look at pre-screen
save(RES, file = here("transmission/data/RES.tgt.Rdata"))

load(file = here("transmission/data/RES.tgt.Rdata"))


summary(RES$NB)
qplot(RES$NB)


## breaks for plotting
nbks <- 5
bks1 <- c(0:5) / 5
bks2 <- bks1/2+0.5
RES[, SP1.cat := cut(SP1, breaks = bks2, include.lowest = TRUE)]
RES[, SE1.cat := cut(SE1, breaks = bks1, include.lowest = TRUE)]
SMY2 <- RES[, .(NB = mean(NB * ratio)), by = .(SE1.cat, SP1.cat)]


ggplot(SMY2, aes(SE1.cat, SP1.cat, fill = NB )) +
  geom_tile() +
  xlab("Pre-screen sensitivity") +
  ylab("Pre-screen specificity") +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_linedraw() +
  labs(fill = "Net benefit")


ggsave(file = here("transmission/plots/p_target_tile_NB.png"), w = 14, h = 7)


ggplot(SMY2, aes(SE1.cat, SP1.cat, fill = NB > 0)) +
  geom_tile() +
  xlab("Pre-screen sensitivity") +
  ylab("Pre-screen specificity") +
  scale_colour_manual(values = c("red", "green")) ## +

ggsave(file = here("transmission/plots/p_target_tile.png"), w = 7, h = 7)


XZ <- (bks1[-1]+bks1[1:5])/2 #SE
YZ <- (bks2[-1] + bks2[1:5]) / 2 # SP
ZZ <- expand.grid(XZ, YZ)
Z <- expand.grid(1:5,1:5)
ZOR <- odds.ratio(ZZ[,1],ZZ[,2])             #OR
ZFR <- (ZZ[, 1] * 0.1 + (1 - ZZ[, 2]) * 0.9) #fraction screened
ZLP <- (ZZ[, 1] * 0.1)/ZFR                   #TBI prev

TXT <- paste0("TBI=", round(1e2 * ZLP), "%\nProportion=", round(1e2 * ZFR), "%\n")
keep <- which(ZLP>0.15)#1:nrow(Z)

GP <- ggplot(SMY2, aes(SE1.cat, SP1.cat, fill = NB)) +
  geom_tile() +
  xlab("Pre-screen sensitivity") +
  ylab("Pre-screen specificity") +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_linedraw() +
  labs(fill = "Net benefit") +
  annotate("text", label = TXT[keep], col = 'yellow', x = Z[keep,1], y = Z[keep,2])
GP

ggsave(GP, file = here("transmission/plots/p_target_tile_NB2.png"), w = 14, h = 7)


## ---- with specific annotations
## se,sp
## 18,95
## 30,70

se1x <- 1 + (18 / 100 - 0.1) / 0.2
sp1x <- 1 + (95 / 100 - 0.55) / 0.1
se2x <- 1 + (30 / 100 - 0.1) / 0.2
sp2x <- 1 + (70 / 100 - 0.55) / 0.1

sz <- 5
cl <- "magenta"
GPA <- GP + annotate(geom = "point", x = se1x, y = sp1x, col = cl, size = sz, shape=1)+
  annotate(geom = "text", x = se1x, y = sp1x, col = cl, size = sz*0.9, label = "1") +
  annotate(geom = "point", x = se2x, y = sp2x, col = cl, size = sz, shape = 1) +
  annotate(geom = "text", x = se2x, y = sp2x, col = cl, size = sz*0.9, label = "2")+
  theme(panel.grid = element_blank(), panel.border = element_blank())
GPA

ggsave(GPA, file = here("transmission/plots/p_target_tile_NB3.png"), w = 14, h = 7)
ggsave(GPA, file = here("transmission/plots/p_target_tile_NB3.pdf"), w = 14, h = 7)
ggsave(GPA, file = here("transmission/plots/p_target_tile_NB3.eps"), w = 14, h = 7)

## Targeting: explicit evaluation of two strategies
TGT <- list()

## group 1
set.seed(sd)
REST <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = FALSE, community = TRUE, posttb = TRUE,
  targeting = TRUE, screen.acc = list(SE1 = 0.18, SP1 = 0.95)
)
gp1 <- copy(REST)

TGT[[1]] <- data.table(
  analysis = "Target group 1",
  ICER = REST[, mean((int.CC - soc.CC) * ratio) / mean(dQ * ratio)],
  ICERr = REST[
    mid.notes < 100 & mid.notes > 30,
    mean((int.CC - soc.CC) * ratio) / mean(dQ * ratio)
  ],
  NB = REST[, 30e3 * mean(dQ * ratio) - mean((int.CC - soc.CC) * ratio)],
  NBr = REST[
    mid.notes < 100 & mid.notes > 30,
    30e3 * mean(dQ * ratio) - mean((int.CC - soc.CC) * ratio)
  ],
  dC = REST[, mean((int.CC - soc.CC) * ratio)],
  dQ = REST[, mean(dQ * ratio)],
  dCr = REST[mid.notes < 100 & mid.notes > 30, mean((int.CC - soc.CC) * ratio)],
  dQr = REST[mid.notes < 100 & mid.notes > 30, mean(dQ * ratio)],
  NB.sd = REST[, 30e3 * sd(dQ * ratio) - sd((int.CC - soc.CC) * ratio)],
  NBr.sd = REST[
    mid.notes < 100 & mid.notes > 30,
    30e3 * sd(dQ * ratio) - sd((int.CC - soc.CC) * ratio)
  ],
  dC.sd = REST[, sd((int.CC - soc.CC) * ratio)],
  dQ.sd = REST[, sd(dQ * ratio)],
  dCr.sd = REST[mid.notes < 100 & mid.notes > 30, sd((int.CC - soc.CC) * ratio)],
  dQr.sd = REST[mid.notes < 100 & mid.notes > 30, sd(dQ * ratio)]
)


## group 2
set.seed(sd)
REST <- PSAloop(
  Niter = Nruns, parms, smpsd, DR,
  static = FALSE, community = TRUE, posttb = TRUE,
  targeting = TRUE, screen.acc = list(SE1 = 0.30, SP1 = 0.70)
)
gp2 <- copy(REST)

TGT[[2]] <- data.table(
  analysis = "Target group 2",
  ICER = REST[, mean((int.CC - soc.CC) * ratio) / mean(dQ * ratio)],
  ICERr = REST[
    mid.notes < 100 & mid.notes > 30,
    mean((int.CC - soc.CC) * ratio) / mean(dQ * ratio)
  ],
  NB = REST[, 30e3 * mean(dQ * ratio) - mean((int.CC - soc.CC) * ratio)],
  NBr = REST[
    mid.notes < 100 & mid.notes > 30,
    30e3 * mean(dQ * ratio) - mean((int.CC - soc.CC) * ratio)
  ],
  dC = REST[, mean((int.CC - soc.CC) * ratio)],
  dQ = REST[, mean(dQ * ratio)],
  dCr = REST[mid.notes < 100 & mid.notes > 30, mean((int.CC - soc.CC) * ratio)],
  dQr = REST[mid.notes < 100 & mid.notes > 30, mean(dQ * ratio)],
  NB.sd = REST[, 30e3 * sd(dQ * ratio) - sd((int.CC - soc.CC) * ratio)],
  NBr.sd = REST[
    mid.notes < 100 & mid.notes > 30,
    30e3 * sd(dQ * ratio) - sd((int.CC - soc.CC) * ratio)
  ],
  dC.sd = REST[, sd((int.CC - soc.CC) * ratio)],
  dQ.sd = REST[, sd(dQ * ratio)],
  dCr.sd = REST[mid.notes < 100 & mid.notes > 30, sd((int.CC - soc.CC) * ratio)],
  dQr.sd = REST[mid.notes < 100 & mid.notes > 30, sd(dQ * ratio)]
)

## combine
TGT <- rbindlist(TGT)

## save out
fwrite(TGT, file=here("transmission/plots/TGT.csv"))
save(gp1, gp2, file = here("transmission/data/gpso.Rdata"))



## === exploring details of group 1:
tmp1 <- gp1[, .(
  group = "group 1",
  blpop, mid.notes,
  ratio, int.to.end,
  soc.CC0, int.CC0,
  soc.CC0inflow, int.CC0inflow,
  soc.CC0inside, int.CC0inside,
  soc.CC0outside, int.CC0outside,
  dCC0 = -soc.CC0 + int.CC0,
  dCC0inflow = -soc.CC0inflow + int.CC0inflow,
  dCC0inside = -soc.CC0inside + int.CC0inside,
  dCC0outside = -soc.CC0outside + int.CC0outside,
  soc.ccases, int.ccases, dccases = -soc.ccases + int.ccases,
  soc.ccasesout, int.ccasesout, dccasesout = -soc.ccasesout + int.ccasesout,
  soc.deaths, int.deaths, ddeaths = -soc.deaths + int.deaths,
  soc.cTPT, int.cTPT, dcTPT = -soc.cTPT + int.cTPT,
  soc.cATT = soc.cATTtp + soc.cATTfp, int.cATT = int.cATTtp + int.cATTfp,
  dcATT = -soc.cATTtp - soc.cATTfp + int.cATTtp + int.cATTfp,
  soc.Cinflow = soc.cInflow, int.Cinflow = int.cInflow
)]


tmp1[, .(87489/blpop,ratio)]

tmp2 <- gp2[, .(
  group = "group 2",
  blpop, mid.notes,
  ratio, int.to.end,
  soc.CC0, int.CC0,
  soc.CC0inflow, int.CC0inflow,
  soc.CC0inside, int.CC0inside,
  soc.CC0outside, int.CC0outside,
  dCC0 = -soc.CC0 + int.CC0,
  dCC0inflow = -soc.CC0inflow + int.CC0inflow,
  dCC0inside = -soc.CC0inside + int.CC0inside,
  dCC0outside = -soc.CC0outside + int.CC0outside,
  soc.ccases, int.ccases, dccases = -soc.ccases + int.ccases,
  soc.ccasesout, int.ccasesout, dccasesout = -soc.ccasesout + int.ccasesout,
  soc.deaths, int.deaths, ddeaths = -soc.deaths + int.deaths,
  soc.cTPT, int.cTPT, dcTPT = -soc.cTPT + int.cTPT,
  soc.cATT = soc.cATTtp + soc.cATTfp, int.cATT = int.cATTtp + int.cATTfp,
  dcATT = -soc.cATTtp - soc.cATTfp + int.cATTtp + int.cATTfp,
  soc.Cinflow = soc.cInflow, int.Cinflow = int.cInflow
)]

save(tmp1,tmp2,file=here('transmission/data/tmporary.Rdata'))



projs <- rbind(
  tmp1[, .(
    group[1],
    dc.m = -mean(ratio * dccases / 70),
    dc.lo = lo(-ratio * dccases / 70), dc.hi = hi(-ratio * dccases / 70),
    dco.m = mean(-ratio * dccasesout / 70),
    dco.lo = lo(-ratio * dccasesout / 70), dco.hi = hi(-ratio * dccasesout / 70),
    dd.m = mean(-ratio * ddeaths / 70),
    dd.lo = lo(-ratio * ddeaths / 70), dd.hi = hi(-ratio * ddeaths / 70),
    da.m = mean(-ratio * dcATT / 70),
    da.lo = lo(-ratio * dcATT / 70), da.hi = hi(-ratio * dcATT / 70),
    dC.m = mean(-ratio * dCC0 / 70) / 1e3,
    dC.lo = lo(-ratio * dCC0 / 70) / 1e3, dC.hi = hi(-ratio * dCC0 / 70) / 1e3,
    dCinflow.m = mean(-ratio * dCC0inflow / 70) / 1e3,
    dCinflow.lo = lo(-ratio * dCC0inflow / 70) / 1e3, dCinflow.hi = hi(-ratio * dCC0inflow / 70) / 1e3,
    dCinside.m = mean(-ratio * dCC0inside / 70) / 1e3,
    dCinside.lo = lo(-ratio * dCC0inside / 70) / 1e3, dCinside.hi = hi(-ratio * dCC0inside / 70) / 1e3,
    dCoutside.m = mean(-ratio * dCC0outside / 70) / 1e3,
    dCoutside.lo = lo(-ratio * dCC0outside / 70) / 1e3,
    dCoutside.hi = hi(-ratio * dCC0outside / 70) / 1e3
  )],
  tmp2[, .(
    group[1],
    dc.m = -mean(ratio * dccases / 70),
    dc.lo = lo(-ratio * dccases / 70), dc.hi = hi(-ratio * dccases / 70),
    dco.m = mean(-ratio * dccasesout / 70),
    dco.lo = lo(-ratio * dccasesout / 70), dco.hi = hi(-ratio * dccasesout / 70),
    dd.m = mean(-ratio * ddeaths / 70),
    dd.lo = lo(-ratio * ddeaths / 70), dd.hi = hi(-ratio * ddeaths / 70),
    da.m = mean(-ratio * dcATT / 70),
    da.lo = lo(-ratio * dcATT / 70), da.hi = hi(-ratio * dcATT / 70),
    dC.m = mean(-ratio * dCC0 / 70) / 1e3,
    dC.lo = lo(-ratio * dCC0 / 70) / 1e3, dC.hi = hi(-ratio * dCC0 / 70) / 1e3,
    dCinflow.m = mean(-ratio * dCC0inflow / 70) / 1e3,
    dCinflow.lo = lo(-ratio * dCC0inflow / 70) / 1e3, dCinflow.hi = hi(-ratio * dCC0inflow / 70) / 1e3,
    dCinside.m = mean(-ratio * dCC0inside / 70) / 1e3,
    dCinside.lo = lo(-ratio * dCC0inside / 70) / 1e3, dCinside.hi = hi(-ratio * dCC0inside / 70) / 1e3,
    dCoutside.m = mean(-ratio * dCC0outside / 70) / 1e3,
    dCoutside.lo = lo(-ratio * dCC0outside / 70) / 1e3,
    dCoutside.hi = hi(-ratio * dCC0outside / 70) / 1e3
  )]
)

projt <- transpose(projs)
projt[, var := names(projs)]
projt[, c("variable", "type") := tstrsplit(var, "\\.")]
projt[, var := NULL]
names(projt)[1:2] <- unlist(projt[1, 1:2])
projt <- projt[!is.na(type)]
projt[, nm := fcase(
  variable == "dc", "TB incidence averted",
  variable == "dco", "TB incidence outside averted",
  variable == "dd", "TB deaths averted",
  variable == "da", "TB treatments averted",
  variable == "dC", "cost averted/K",
  variable == "dCinflow", "inflow cost averted/K",
  variable == "dCinside", "inside cost averted/K",
  variable == "dCoutside", "outside cost averted/K"
)]
projt <- melt(projt[, .(nm, type, `group 1`, `group 2`)], id = c("nm", "type"))
projt <- dcast(projt, nm + variable ~ type, value.var = "value")
setcolorder(projt, c("variable", "nm", "m", "lo", "hi"))
setkey(projt, variable, nm)


fwrite(projt, file = here("transmission/plots/proj.csv"))

projt[variable == "group 1"]



## CE plot?
CE <- rbind(
  RESA[analysis == "base case"][, .(
    scenario = "Base case",
    dC.m = mean(ratio * (-soc.CC + int.CC) ) ,
    dC.lo = lo(ratio * (-soc.CC + int.CC)), dC.hi = hi(ratio * (-soc.CC + int.CC)),
    dQ.m = mean(ratio * dQ),
    dQ.lo = lo(ratio * dQ), dQ.hi = hi(ratio * dQ)
    )],
  ## RESA[analysis == "PrisonEscort"][, .(
  ##   scenario = "No escort costs",
  ##   dC.m = mean(ratio * (-soc.CC + int.CC) ) ,
  ##   dC.lo = lo(ratio * (-soc.CC + int.CC)), dC.hi = hi(ratio * (-soc.CC + int.CC)),
  ##   dQ.m = mean(ratio * dQ),
  ##   dQ.lo = lo(ratio * dQ), dQ.hi = hi(ratio * dQ)
  ## )],
  gp1[, .(
    scenario = "Targeted, group 1",
    dC.m = mean(ratio * (-soc.CC + int.CC) ) ,
    dC.lo = lo(ratio * (-soc.CC + int.CC)), dC.hi = hi(ratio * (-soc.CC + int.CC)),
    dQ.m = mean(ratio * dQ),
    dQ.lo = lo(ratio * dQ), dQ.hi = hi(ratio * dQ)
  )],
  gp2[, .(
    scenario = "Targeted, group 2",
    dC.m = mean(ratio * (-soc.CC + int.CC) ) ,
    dC.lo = lo(ratio * (-soc.CC + int.CC)), dC.hi = hi(ratio * (-soc.CC + int.CC)),
    dQ.m = mean(ratio * dQ),
    dQ.lo = lo(ratio * dQ), dQ.hi = hi(ratio * dQ)
  )]
)
CE

fwrite(CE, file = here("transmission/plots/CEdata.csv"))

## ======== CE simple plot

library(ggplot2)
library(scales)
library(geomtextpath)


CE <- read.csv(here("transmission/plots/CEdata.csv"))
nnmz <- c(
  "No targeting",
  "High/medium TB incidence country-of-birth",
  "High/medium TB incidence country-of-birth,\nor a history of homelessness, injecting drug-use,\nor problematic alcohol use"
)
CE$scenario <- nnmz
CE$scenario <- factor(CE$scenario,levels=nnmz)


colz <- c("#000000", "#E69F00", "#56B4E9","#009E73")


ggplot(CE, aes(
  x = dQ.m, y = dC.m / 1e6,
  xmin = dQ.lo, xmax = dQ.hi,
  ymin = dC.lo / 1e6, ymax = dC.hi / 1e6,
  col = scenario
)) +
  geom_textabline(
    intercept = 0, slope = 30e-3, col = 2, lty = 2,
    label = "30K GBP/QALY", vjust = 1.1, fontface = "italic"
  ) +
  geom_abline(intercept = 0, slope = 0, col = 2, lty = 2) +
  geom_vline(xintercept = 0, col = 2, lty = 2) +
  geom_errorbar(width = 1) +
  geom_errorbarh(height = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = colz) +
  theme_classic() +
  ggpubr::grids() +
  scale_x_continuous(label = comma) +
  scale_y_continuous(label = comma) +
  xlab("Discounted quality-adjusted life-years (QALY) gained over 70 years") +
  ylab("Discounted incremental cost over 70 years (million GBP)") +
  annotate(x = 1e4, y = 10, col = 2, geom = "text", label = "Costs more", fontface = "italic") +
  annotate(x = 1e4, y = -10, col = 2, geom = "text", label = "Costs less", fontface = "italic") +
  geom_segment(aes(x = 1e4, y = 20, xend = 1e4, yend = 60),
    col = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_segment(aes(x = 1e4, y = -20, xend = 1e4, yend = -60),
    col = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  )  +
  theme(legend.position = "top",legend.title = element_blank())



ggsave(file = here("transmission/plots/CEsimple.png"), w = 9, h = 6)
ggsave(file = here("transmission/plots/CEsimple.pdf"), w = 9, h = 6)



## ========== Authors only:
upload.to.sheets(here("outdata/"), "DRS.csv", shid)

flz <- c("SA.csv", "SA.full.csv", "tabout.csv", "TGT.csv", "SDRtab.csv", "proj.csv")
for(fn in flz)
  upload.to.sheets(here('transmission/plots/'),fn,shid)
