## TODO questions:
##
## - stool/sputum
## - flag assumption = groups for SA or pending data
##

## ========= UTILITIES ===============
logit <- function(x) log(odds(x))
ilogit <- function(x) iodds(exp(x))
AOR <- function(base, OR) ilogit(log(OR) + logit(base))
AOR2 <- function(base, OR) OR * base / (1 - base + (base * OR))
odds <- function(x) x / (1 - x)
iodds <- function(x) x / (1 + x)
lo <- function(x) quantile(x, probs = 0.025, na.rm = TRUE)
hi <- function(x) quantile(x, probs = 1 - 0.025, na.rm = TRUE)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

brkt <- function(M, L, H, ndp = 0) {
  paste0(
    round(M, ndp), " (",
    round(L, ndp), " - ",
    round(H, ndp), ")"
  )
}
gm <- function(x) exp(mean(log(x))) # geometric mean
gh <- function(x) glue(here(x))

## ======= COMBINED LABELLER ===========

## additional labels from data (may overwrite some initial version currently)
AddDataDrivenLabels <- function(D) {
  
  # Setting some SOC/INT parameters by TB status
  # For ATT: TB symptoms -> Clinical suspicion -> abnormal x-ray -> GP assessment -> tb dx
  # For TPT: IGRA tested -> IGRA positive
  # For setting SOC parameters
  # assuming a baseline some fraction lower than INT for TB screening & prop.xray
  # could also use any TB symptom sensitivity/specificity for TB screening assuming no xray in SOC
  D[, soc.fac := rbeta(nrow(D), 5 / 1.5, 10)]

  #' `TB symptom screening`
  D[, soc.prop.tb.sympt.screen := int.prop.tb.sympt.screen * soc.fac]

  # --------------------------------------------------------------------------------------
  # TB presumption/TB positive symptom screening
  D[, soc.prop.presumtive.tb := ifelse(tb == "TBD", sens.symptom, 1 - spec.symptom)]
  D[, int.prop.presumtive.tb := ifelse(tb == "TBD", sens.any.abn.xray, 1 - spec.any.abn.xray)]

  # Fraction getting a chest x-ray under SOC
  D[, soc.prop.xray := int.prop.xray * soc.fac]

  #' `Diagnostic accuracy of Xpert MTB/RIF for pulmonary TB in adults`
  D[, soc.prop.prev.tb.tx.symp.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]
  D[, soc.prop.abn.xray.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]
  D[, soc.prop.no.xray.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]
  D[, soc.prop.no.prev.tb.dx.symp.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]

  D[, int.prop.prev.tb.tx.symp.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]
  D[, int.prop.abn.xray.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]
  D[, int.prop.no.xray.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]
  D[, int.prop.no.prev.tb.dx.symp.tb.dx := ifelse(tb == "TBD", sens.xpert, 1 - spec.xpert)]

  #' `LTBI pathway stuff`
  # For TPT: IGRA tested -> IGRA positive -> stying 3+ months -> starting TPT

  # IGRA tested
  D[, soc.prop.igra.tested := ifelse(tb == "TBI", soc.prop.igra.tested, 0)]
  D[, int.prop.igra.tested := ifelse(tb == "TBI", int.prop.igra.tested, 0)]

  # IGRA positive
  D[, soc.prop.igra.test.positive := ifelse(tb == "TBI", 1, 0)]
  D[, int.prop.igra.test.positive := ifelse(tb == "TBI", 1, 0)]

  # # staying longer that 3 months
  D[, int.prop.staying.o3.months := ifelse(tb == "TBI", int.prop.staying.o3.months, 0)]
  D[, soc.prop.staying.o3.months := ifelse(tb == "TBI", soc.prop.staying.o3.months, 0)]

  # #remaining in the UK: only for TBI
  D[, int.prop.staying.u3.months.uk := ifelse(tb == "TBI", int.prop.staying.u3.months.uk, 0)]
  D[, soc.prop.staying.u3.months.uk := ifelse(tb == "TBI", soc.prop.staying.u3.months.uk, 0)]

  # Costs
  D[, cost.tb.sympt.screen := ucost.tb.sympt.screen + ucost.overheads]
  D[, pxray := 1]
  D[, cost.chest.xray := pxray * (ucost.chest.xray + ucost.prison.escort)]
  D[, cost.prison.isolation := pIsolation * ucost.prison.cell.isolation]
  D[, cost.prison.gp.assessment := ucost.prison.gp.assess]
  D[, cost.contact.management := smear.positive * (cost.prison.isolation + nContacts * ucost.contact.tracing)]
  D[, cost.nhs.tb.service := ucost.nhs.tb.service]
  D[, cost.prison.escort := ucost.prison.escort]
  D[, cost.attending.nhs.tb.service := ucost.nhs.tb.service + ucost.prison.escort]
  D[, cost.tb.evaluation := cost.attending.nhs.tb.service + ucost.tb.investigations]
  D[, cost.tb.investigations := ucost.tb.investigations]
  D[, cost.att.initiation := cost.attending.nhs.tb.service]
  D[, cost.tpt.initiation := cost.attending.nhs.tb.service]
  D[, ucost.dots := ucost.att.dots]
  D[, ucost.att.dots := NULL]
  D[, ucost.tpt.opd.visit := ucost.dstb.opd.visit]
  D[, durTPT := DurTPT]
  D[, DurMDRTB := DurDSTB / DurMDRTBfactor]
  D[, IncompDurDSTB := DurDSTB * IncompTxfactor]
  D[, IncompDurMDRTB := DurMDRTB * IncompTxfactor]
  D[, IncompDurTPT := durTPT * IncompTxfactor]
  D[, cost.inpatient := (pDSTB * smear.positive * DurDSTBIsolation * (ucost.dstb.ipd + ucost.prison.bedwatch) +
    (1 - pDSTB) * smear.positive * DurMDRTB2Isolation * (ucost.mdrtb.ipd + ucost.prison.bedwatch) +
    (1 - pDSTB) * (1 - smear.positive) * DurMDRTBIsolation * (ucost.mdrtb.ipd + ucost.prison.bedwatch))]
  D[, cost.att.complete := pDSTB * dstb.visits * (ucost.dstb.opd.visit + ucost.prison.escort) + # DSTB outpatient visits
    (1 - pDSTB) * mdrtb.visits * (ucost.mdrtb.opd.visit + ucost.prison.escort) + # MDRTB outpatient visits
    pDSTB * DurDSTB * ucost.dsatt.drugs + # DSTB drugs
    (1 - pDSTB) * DurMDRTB * ucost.mdratt.drugs + # MDRTB drugs
    ucost.dots * (pDSTB * DurDSTB + (1 - pDSTB) * DurMDRTB) + # DOTS
      cost.inpatient]
  D[, cost.att.incomplete := IncompDurDSTB / DurDSTB * pDSTB * (
    dstb.visits * (ucost.dstb.opd.visit + ucost.prison.escort) +
      DurDSTB * (ucost.dsatt.drugs + ucost.dots)) +
    IncompDurMDRTB / DurMDRTB * (1 - pDSTB) * (
      mdrtb.visits * (ucost.mdrtb.opd.visit + ucost.prison.escort) +
        DurMDRTB * (ucost.mdratt.drugs + ucost.dots)) +
    cost.inpatient]
  D[, cost.tpt.complete := durTPT * (ucost.ltbi.drugs + ucost.dots) + TPT.visits * (ucost.tpt.opd.visit + ucost.prison.escort)]
  D[, cost.tpt.incomplete := IncompDurTPT / durTPT * (durTPT * (ucost.ltbi.drugs + ucost.dots) + TPT.visits * (ucost.tpt.opd.visit + ucost.prison.escort))]
  D[, cost.igra.test := ucost.igra.test]
  D[, nhs.referral := 1]
  D[, soc.fac := NULL] # remove temporary variable
}

## combined function to add the labels to the tree prior to calculations
MakeTreeParms <- function(D, P) {
  ## -- use of other functions
  AddDataDrivenLabels(D)
}

## ======= EPIDEMIOLOGY ===========

makeAttributes <- function(D) {
  nrep <- nrow(D)
  D[, id := 1:nrep]
  fx <- list(age = agelevels, isoz = c("GBR"), tb = tblevels)
  cofx <- expand.grid(fx)
  cat("Attribute combinations used:\n")
  print(cofx)
  D <- D[rep(1:nrow(D), each = nrow(cofx))] # expand out data
  D[, names(cofx) := cofx[rep(1:nrow(cofx), nrep), ]]

  # ## --- age
  # D[,value:=ifelse(age=='15-64',0.8,1-0.8)] #NOTE value first set here
  D[, value := 1] # Not using age splits for now
  ## --- TB
  ## ## Using 3 levels of TB: active TB disease, TB infection, and noTB (no active TB/TB infection)
  D[, ltbi.prev := ltbi.prev]
  D[, prog.tb := prog.tb]

  D[tb == "TBD", value := value * ltbi.prev * prog.tb]
  D[tb == "TBI", value := value * (ltbi.prev - ltbi.prev * prog.tb)]
  D[tb == "noTB", value := value * (1 - ltbi.prev)] 
  return(D)
}


## function for generating random sample of costs
MakeCostData <- function(csts, # base data table of cost data
                         nrep, # number of replicates being used in PSA
                         anmz = NULL # attribute names (if any)
) {
  if (nrow(csts[cost.sd > 0 & cost.m == 0]) > 0) warning(paste0("Some cost input variables have zero mean & SD>0. These will be treated as fixed variables:\n", paste0(csts[cost.sd > 0 & cost.m == 0, cost], collapse = "\n")))
  if (is.null(anmz) & any(csts[, table(cost)] > 1)) warning("Some cost names occur >1 times, but no attributes have been specified! This is unlikely to do what you want.")
  csts[cost.m > 0, gmsc := cost.sd^2 / cost.m]
  csts[!is.na(gmsc) & gmsc > 0, gmk := cost.m / gmsc]
  NR <- nrow(csts)
  csts <- csts[rep(1:NR, nrep)]
  csts[, id := rep(1:nrep, each = NR)]
  csts[, rnd := !is.na(gmsc) & !is.na(gmk) & gmk > 0 & gmsc > 0]
  csts[rnd == TRUE, value := rgamma(sum(rnd), shape = gmk, scale = gmsc)] # random sample from gamma distribution
  csts[rnd != TRUE, value := cost.m] # fixed values
  ## csts[,cnms:=paste0('c_',cost)]
  csts[, cnms := paste0(cost)]
  F <- "id "
  if (!is.null(anmz)) F <- paste0(F, "+ ", paste(anmz, collapse = "+")) # split out by attributes if included
  F <- paste0(F, " ~ cnms")
  dcast(csts, as.formula(F), value.var = "value") # id ~ cnms
}
## NOTE
## if attributes are included, all costs need to be specified by them even if this means duplicating those without dependence


