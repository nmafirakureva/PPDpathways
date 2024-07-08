## this if for various checks and any debugging work
DRS[quantity=='cost' & outcome=='notx'] #
DRS[quantity=='cost' & outcome=='att'] #
DRS[quantity=='cost' & outcome=='tpt'] #


## incoming state?
## inflow split: parms from doing library(ecrins); data(parms)
parms$foi <- qfun(runif(1), hyperparms[["foi"]]) # median
parms <- revise.instates(parms)                  #uses FOI
IV <- with(data = parms, {
  c(
    parm_frac_SD, parm_frac_CD,
    parm_frac_E, parm_frac_L, parm_frac_epTB, parm_frac_lpTB,
    parm_frac_U,parm_frac_ATT
  )
})
sum(IV)

IFS <- with(data = parms, {
  c(
    TBD = parm_frac_SD + parm_frac_CD,
    TBI = parm_frac_E + parm_frac_L + parm_frac_epTB + parm_frac_lpTB,
    noTB = parm_frac_U
  )
})
IFS

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

## --------- back of envelop CE ----------
## want ~1M with effects off

(dLE <- (1-exp(-40*3.5/1e2))/(3.5/100)) #discounted life-exp
cdr <- 0.9                              #CDR
drn <- 3                                #untreated TB durn
tef <- 0.17                             #efficacy TPT in IGRA+
pgn <- 0.1                              #progression risk
txf <- 0.05                             #CFR on tx
cfr <- 0.4 * (1-cdr)+0*cdr*txf            #CFR
yil <- drn * (1-cdr)                    #duration ill
dec <- 1/3                              #qol decrement
tptcov <- 0.3                          #TPT cov in TBI+ INT (assumed 0 SOC)
attcov <- 0.4                           #ATT cov in TBD+ INT (assumed 0 SOC)
(UCatt <- (avcint/avoint)[1])             #unit cost ATT


dq <- IFS[1] * attcov * (cfr) * dLE + # mortality = TBD x fraction immediately found x cfr=(CFRx(1-CDR)) x dLE
  IFS[2] * (tptcov) * (pgn*tef) * (cfr) * dLE + #mortality = TBI x fraction TPT x TPT eff x progn x cfr x dLE
  IFS[1] * dec * yil +  #qol=TBD x dec x (durn as T x 1-CDR)
  IFS[2] * dec * (pgn*tef) * (yil) #prevented cases qol
(dc <- sum(avcint-avcsoc))         #difference per av. person: 'on the door costs'
dc <- dc - IFS[2] * tptcov * (pgn*tef) * cdr * UCatt #UC saved of treating detected now prevented
dc <- dc - IFS[1] * attcov * cdr * UCatt #UC saved of treating earlier (rather than later)
(dc/dq)                              #90 M



1e4*IFS[1] * attcov * (cfr) * dLE # mortality = TBD x fraction immediately found x cfr=(CFRx(1-CDR)) x dLE
1e4*IFS[2] * (tptcov) * (pgn * tef) * (cfr) * dLE  # mortality = TBI x fraction TPT x TPT eff x progn x cfr x dLE
1e4*IFS[1] * dec * yil  # qol=TBD x dec x (durn as T x 1-CDR)
1e4*IFS[2] * dec * (pgn * tef) * (yil)

## TODO check ATT$ for TBI > ATT$ TBD

DRS[quantity=='check' & outcome=='notx',] #

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


## impact/outcome checks
## inputs
parms3 <- parms
parms3$mHR <- 1 # for meaningful CFR checks
res3 <- run.HE.socint(parms3, tmp, 1, zero.nonscreen.costs = TRUE)

## checks
res3[,soc.CC0/soc.CC] #2.6 for info
res3[,soc.ccasesout/soc.ccases] #over half cases out?

## analytical dLYL/deaths if mortality constant:
((1-exp(-parms$LifeExp *parms$disc_rate))/parms$disc_rate)*
  (1-exp(-parms$disc_rate*res3$int.to.end))/(parms$disc_rate*res3$int.to.end)
res3[,soc.dLYL/soc.deaths]   #dLYL/deaths -- pretty close

res3[,soc.deaths/soc.ccases] #CFR OK; too high NOTE this was because of the post-TB mHR
res3[,soc.cATTtp/soc.ccases] #86% OK as CDR

## PERFECT vs PERFECTLY AWFUL
parms3$tptHR <- 0; parms3$tpt_drn <- 50#potent, durable tpt
parms3$inflow_toATT_TB0 <- parms3$inflow_toATT_L0 <- parms3$inflow_toATT_no0 <- 0 #no ATT
parms3$inflow_toTPT_TB0 <- parms3$inflow_toTPT_L0 <- parms3$inflow_toTPT_no0 <- 0 # no TPT
parms3$inflow_toATT_TB1 <- parms3$inflow_toATT_L1 <- parms3$inflow_toATT_no1 <- 0 # no ATT
parms3$inflow_toTPT_TB1 <- parms3$inflow_toTPT_L1 <- parms3$inflow_toTPT_no1 <- 0 # no TPT
res3a <- run.HE.socint(parms3, tmp, 1, zero.nonscreen.costs = TRUE, ignore.tree.parms = TRUE) #run
parms3$inflow_toATT_TB1 <- 1 # PERFECT!
parms3$inflow_toTPT_L1 <- 1  # PERFECT!
res3b <- run.HE.socint(parms3, tmp, 1, zero.nonscreen.costs = TRUE, ignore.tree.parms = TRUE) #run
## comparison
CF <- data.table(
  res3a[,.(soc.deaths,soc.qoldec,soc.dLYL,soc.ccases,soc.ccasesout,soc.cTPT,soc.cATTtp)],
  res3b[, .(int.deaths, int.qoldec, int.dLYL, int.ccases, int.ccasesout, int.cTPT, int.cATTtp)]
)
CF <- CF[, .(soc.deaths, int.deaths, soc.qoldec, int.qoldec, soc.dLYL, int.dLYL,
             soc.ccases, int.ccases, soc.ccasesout,int.ccasesout,
             soc.cTPT, int.cTPT, soc.cATTtp, int.cATTtp)]
CF

## TODO not prev TB?
## NOTE still all the latent in prison

## test problem catching
parms3$CDR <- 1
res3b <- run.HE.socint(parms3, tmp, 1, zero.nonscreen.costs = TRUE, ignore.tree.parms = TRUE) # run

res3b$problem
1

## ====================== DEBUGGING =========

## RUN REAL PSA:
RES <- PSAloop(Niter=2e3,parms,smpsd,DR,zero.nonscreen.costs=FALSE,verbose=FALSE) #TODO update this as we go

summary(RES$problem)

## TODO seeing negative deaths
## capture inputs to analyse

summary(RES)

RES[soc.deaths<0]
which(RES$soc.deaths<0)


## SO <- runmodelsafely(tt,parms)
## SO$error_val

## SO <- runmodelsafely(tt, -11)
## SO
## SO$error_val



## J <- 1222
## parms <- revise.flow.parms(parms, smpsd, J) # TODO BUG
## newp <- uv2ps(runif(length(hyperparms)), hyperparms) # natural history
## for (nm in names(newp)) parms[[nm]] <- newp[[nm]] # safety
## ans <- run.HE.socint(parms, DR, J, zero.nonscreen.costs = FALSE)

## ## get RES use in parms
## J <- which(RES$soc.deaths < 0)[3]
## PMZ <- as.list(RES[J])
## PMZ[["parm_ifrac_prevTPT"]] <- c(PMZ$parm_ifrac_prevTPT1, PMZ$parm_ifrac_prevTPT3, PMZ$parm_ifrac_prevTPT3)
## PMZ[["parm_ifrac_prevTPT1"]] <- PMZ[["parm_ifrac_prevTPT2"]] <- PMZ[["parm_ifrac_prevTPT3"]] <- NULL
## PMZ[["parm_init_PPD"]] <- c(PMZ$parm_init_PPD1, PMZ$parm_init_PPD3, PMZ$parm_init_PPD3,
##                             PMZ$parm_init_PPD4, PMZ$parm_init_PPD5)
## PMZ[["parm_init_PPD1"]] <- PMZ[["parm_init_PPD2"]] <- PMZ[["parm_init_PPD3"]] <- PMZ[["parm_init_PPD4"]] <- PMZ[["parm_init_PPD5"]]  <- NULL
## for(i in 1:26) PMZ[[1]] <- NULL
## PMZ <- revise.flow.parms(PMZ, smpsd, J) # TODO BUG
## ans <- run.HE.socint(PMZ, DR, J, zero.nonscreen.costs = FALSE)



save(RES,file='~/Downloads/RES.Rdata')

load(file = "~/Downloads/RES.Rdata")

RES[, bad.soc := ifelse(soc.deaths < 0, 2, 1)]

nm <- names(parms)
nm <- nm[!grepl('uc',nm)] #non-flow parms excluding unit costs
nfs <- nm[grepl("inflow_", nm)] # inflow parms
nm <- setdiff(nm, nfs)
its <- nm[grepl("parm_", nm)]
its <- setdiff(its, c("parm_ifrac_prevTPT"))
nm <- setdiff(nm, its)
nm <- setdiff(nm, c(
  "int_time", "disc_rate", "hrqol", "LifeExp", "mort",
  "att_time", "tpt_drn", "late_post_time","staticfoi"
))
prns <- nm[1:10]
nm <- setdiff(nm, prns)
colz <- RES$bad
shps <- colz; shps[shps!=1] <- 16
az <- ifelse(colz == 1, 0.3, 1)

## parms
## png('~/Dropbox/Holocron/tmp/pairs.png',h=1000,w=1000)
pairs(as.matrix(RES[, ..nm]), col = alpha(colz, az), cex = colz, pch = shps)
## dev.off()
## NOTE hrqol & mHR seem worst

## prns
pairs(as.matrix(RES[, ..prns]), col = alpha(colz, az), cex = colz, pch = shps)
## NOTE no variation

## its
pairs(as.matrix(RES[, ..its]), col = alpha(colz, az), cex = colz, pch = shps)
## NOTE no variation

## nfs
pairs(as.matrix(RES[, ..nfs]), col = alpha(colz, az), cex = colz, pch = shps)
## NOTE no variation


x <- function(i) {
  if (i < 10) warning("A warning")
  i
}
xl <- function(i) {
  list(x(i),0)
}



tryCatch(x(5), warning = function(w) {
  return(list(x(5), w))
})

tryCatch(xl(10), warning = function(w) {
  list(x(5), ifelse(is.null(w),0,1))
})

tryCatch(x(5), warning = function(w) {
  list(x(5), ifelse(!is.null(w), 1, 0))
})




laus <- function(x) {
  r <-
    tryCatch(
      withCallingHandlers(
        {
          error_text <- "No error."
          list(value = hurz(x), error_text = error_text)
        },
        warning = function(e) {
          error_text <<- trimws(paste0("WARNING: ", e))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        return(list(value = NA, error_text = trimws(paste0("ERROR: ", e))))
      },
      finally = {
      }
    )

  return(r)
}

## [[1]]
## [1] 5
##
## [[2]]
## <simpleWarning in x(5): A warning>
