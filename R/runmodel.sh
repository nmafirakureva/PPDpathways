#!/bin/bash
# NOTE (after changing shell flag in ppd_simple.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi dscr, cdr (higher cdr for incidence), txd (completion of ATT/TPT included)

R --slave --vanilla --args <ppd_simple.R AllattendNHS & R --slave --vanilla --args <ppd_simple.R noltfu &
R --slave --vanilla --args <ppd_simple.R DOTsCost & R --slave --vanilla --args <ppd_simple.R XrayCost &
R --slave --vanilla --args <ppd_simple.R FUVisitsCost & R --slave --vanilla --args <ppd_simple.R InpatientCost &
R --slave --vanilla --args <ppd_simple.R none


