#!/bin/bash
# NOTE (after changing shell flag in ppd_simple.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi dscr, cdr (higher cdr for incidence), txd (completion of ATT/TPT included)

R --slave --vanilla --args <ppd_simple.R screenAll & R --slave --vanilla --args <ppd_simple.R noltfu &
# R --slave --vanilla --args <ppd_simple.R presumptiveTB & R --slave --vanilla --args <ppd_simple.R prisonGP &
# R --slave --vanilla --args <ppd_simple.R clinicalSuspicion & R --slave --vanilla --args <ppd_simple.R attendNHS &
# R --slave --vanilla --args <ppd_simple.R startATT & 
R --slave --vanilla --args <ppd_simple.R DOTsCost & R --slave --vanilla --args <ppd_simple.R XrayCost &
R --slave --vanilla --args <ppd_simple.R FUVisitsCost & R --slave --vanilla --args <ppd_simple.R PrisonEscort & 
R --slave --vanilla --args <ppd_simple.R none


