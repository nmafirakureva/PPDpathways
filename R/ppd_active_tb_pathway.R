# rm(list=ls())
library(here)
library(HEdtree)
library(discly)
library(data.tree)
library(data.table)
library(glue)
## NOTE these packages are only needed if wanting to output graphs etc
library(BCEA)
library(ggplot2)
library(scales)

# set_here('/Users/nyashamarfirakureva/My Drive/PPD modelling/Modelling')

## === outcomes subtree ===
tb <- txt2tree(here('indata/4_TB_outcomes.txt')) # tb dx
notbdxo <- txt2tree(here('indata/tbnotx.txt')) # no tx

# notb <- txt2tree(here('indata/noTB.outcomes.txt')) # no tb
tpt <- txt2tree(here('indata/3_TPT_outcomes.txt')) # tpt

## default prob/cost:
tb$Set(p=1)
notbdxo$Set(p=1)
tpt$Set(p=1)
tb$Set(cost=0)
notbdxo$Set(cost=0)
# notb$Set(cost=0)
tpt$Set(cost=0)


# TB outcomes
tbtxo <- Node$new('TB')
tbtxo$AddChildNode(tb)

## -- no tx:
notbdxo$`No TB treatment`$Dies$p <- 'p.cfr.notx'
notbdxo$`No TB treatment`$Survives$p <- '1-p.cfr.notx'

# no tb outcomes
notb <- Node$new('Not TB')
ltbi <- notb$AddChild('LTBI pathway')

## restrict to no TB tx (rather than dx)
notbtxo <- top(notbdxo) #remove top

## ====== function to add outcomes & counters
AddOutcomes <- function(D){
  ## === cost and probs (defaults)
  D$Set(p=1)
  D$Set(cost=0)

  ## === merge 'New people in PPDs' with 'LTBI screening pathway' to create final tree ===
  # MergeByName(D,notb,'Not TB', leavesonly = TRUE) # first add a branch to 'No TB' for easy merging in the next step
  MergeByName(D,No_urgent_GP_referral,'No urgent prison GP referral',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'Refer to TB services',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'did not Attend',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'No clinical suspicion of TB',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'No prison GP assessment',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'Not assessed at next GP appointment',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'Advice',leavesonly = TRUE)
  MergeByName(D,ltbi_pathway,'LTBI pathway',leavesonly = TRUE)

  ## final outcomes
  MergeByName(D,tpt,'TPT outcomes',leavesonly = TRUE)
  MergeByName(D,tbtxo,'TB',leavesonly = TRUE)

  ## ===========  other counters
  ## check
  D$Set(check=1)
  D$Set(check=0,filterFun=function(x) length(x$children)>0)

  ## TB screening
  D$Set(screened=0)
  D$Set(screened=1,filterFun=function(x) x$name=='TB screening')

  ## TB dx
  D$Set(prevtb=0)
  D$Set(prevtb=1,filterFun=function(x) x$name=='TB')

  ## ATT courses
  D$Set(att=0)
  D$Set(att=1,
        filterFun=function(x)x$name=='ATT')

  ## tested on IGRA
  D$Set(igra=0)
  D$Set(igra=1,filterFun=function(x) x$name=='IGRA test')

  ## TPT courses
  D$Set(tpt=0)
  D$Set(tpt=1,
        filterFun=function(x)x$name=='gets TPT')

  return(D)
}

# main tree structures
# new_people_in_PPDs <- txt2tree(here('indata/1.2_new_people_in_PPDs.txt'))
Refer_to_TB_services <- txt2tree(here('indata/1.2.1_Refer_to_TB_services.txt'))
No_urgent_GP_referral <- txt2tree(here('indata/1.2.2_No_urgent_GP_referral.txt'))
ltbi_pathway <- txt2tree(here('indata/2_LTBI_pathway.txt'))
active_tb_pathway <- txt2tree(here('indata/5b_People_who_develop_symptoms_in_PPDs.txt'))

## merge in extras, make model of care branches, write out
tempTree <- AddOutcomes(active_tb_pathway)

## === SOC
SOC_ATB <- Node$new('Standard of care pathway')
SOC_ATB$AddChildNode(tempTree)

SOC_ATB$name <- 'Standard of care pathway'


## # defining tree quantities 
## # p = probability/proportion
## # cost = cost
## # screen = screened for TB
## # igra = received IGRAS test
## # tpt = tpt initiation 
## # prevtb = coprevalent TB diagnosed
## # att = anti-TB treatments

tree2file(SOC_ATB,filename = here('indata/CSV/SOC_ATB.csv'),
          'p','cost', 'prevtb','att', 'igra', 'tpt','check')

labdat <- c('p','cost','cost.att','cost.tpt', 'att','igra', 'tpt','prevtb','check')


## create version with probs/costs
fn <- here('indata/CSV/SOC_ATB1.csv')
if(file.exists(fn)){
  ## read
  labz <- fread(fn)
  labz$p <- gsub('<3','u3.',labz$p)
  labz$p <- gsub('3\\+', 'o3.', labz$p)
  # labz$cost <- gsub('<','u',labz$cost)
  LabelFromData(SOC,labz[,..labdat]) #add label data
  ## NOTE checks need redoing
  SOC$Set(check=1)
  SOC$Set(check=0,filterFun=function(x) length(x$children)>0)
  ## save out
  tree2file(SOC,filename = here('indata/CSV/SOC_ATB2.csv'),
            'p','cost','cost.screen',	'cost.tpt','cost.att','screen', 'igra', 'tpt','prevtb','att','check')
}


## === INT
INT_ATB <- Clone(SOC_ATB)
INT_ATB$name <- 'Intervention model'

tree2file(INT,filename = here('indata/CSV/INT_ATB.csv'),
          'p','cost', 'screen', 'igra', 'tpt', 'prevtb','att','check')


## create version with probs/costs
fn <- here('indata/CSV/INT_ATB1.csv')
if(file.exists(fn)){
  ## read
  labz <- fread(fn)
  labz$p <- gsub('<3','u3.',labz$p)
  labz$p <- gsub('3\\+', 'o3.', labz$p)
  # labz[,cost:=q]; labz[,q:=NULL]     #rename
  LabelFromData(INT_ATB,labz[,..labdat]) #add label data
  ## NOTE checks need redoing
  INT_ATB$Set(check=1)
  INT_ATB$Set(check=0,filterFun=function(x) length(x$children)>0)
  ## save out
  tree2file(INT_ATB,filename = here('indata/CSV/INT_ATB2.csv'),
            'p','cost','cost.screen',	'cost.tpt','cost.att','screen', 'igra', 'tpt','prevtb','att','check')
}


## make functions
fnmz <- labdat[-1]

SOC_ATB.F <- makeTfuns(SOC_ATB,fnmz)
INT_ATB.F <- makeTfuns(INT_ATB,fnmz)

## running all function
runallfuns <- function(D,arm='all'){
        done <- FALSE
        if('SOC_ATB' %in% arm | arm[1]=='all'){
                cat('Running functions for SOC_ATB:\n')
                for(nm in names(SOC_ATB.F)){
                        snm <- gsub('fun','',nm)
                        snma <- paste0(snm,'.soc')
                        D[[snma]] <- SOC_ATB.F[[nm]](D)
                        cat('...',snm,' run...\n')
                        done <- TRUE
                }
        }
        if('INT_ATB' %in% arm | arm[1]=='all'){
                cat('Running functions for INT_ATB:\n')
                for(nm in names(INT_ATB.F)){
                        snm <- gsub('fun','',nm)
                        snma <- paste0(snm,'.int')
                        D[[snma]] <- INT_ATB.F[[nm]](D)
                        cat('...',snm,' run...\n')
                        done <- TRUE
                }
        }


        if(!done)stop('Functions not run! Likely unrecognised arm supplied.')
        return(D)
}

## checking
vrz.soc <- showAllParmz(SOC_ATB)
vrz.int <- showAllParmz(INT_ATB)

vrz <- c(vrz.soc,
         vrz.int
)

vrz <- unique(vrz) #NOTE

# write out as csv
write.csv(data.frame(vrz),here('indata/CSV/parmz_atb.csv'),row.names=FALSE)

# save all as R.data

# Save an object to a file
save.image(file = here("outdata/temp.RData"))

# Restore the object
load(here("outdata/temp.RData"))

A <- makeTestData(5e3,vrz)

# Switch different parms on and off to test
A$prop.prev.tb.dx <- 1
A$prop.prev.tb.dx.started.att <- 1
A$prop.prev.tb.dx.still.on.att <- 1
A$prop.prev.tb.dx.on.attend.referral <- 1
# A$prop.prev.tb.tx.symp <- 1
# A$prop.prev.xray <- 1
# A$prop.abnormal.xray <- 1
# A$prop.abn.xray.nhs.referral <- 1
# A$prop.no.prev.tb.tx.symp <- 0
# A$prop.no.prev.tb.tx.symp.gp.assess <- 0
# A$prop.prev.xray <- 0

## checks
any(SOC_ATB.F$checkfun(A)!=1) #NOTE OK
any(round(SOC_ATB.F$checkfun(A))!=1) # Looks like there are some rounding errors

# INT_ATB.F$checkfun(A) #NOTE OK
all(abs(SOC_ATB.F$checkfun(A)-1)<1e-10) #NOTE OK
all(abs(INT_ATB.F$checkfun(A)-1)<1e-10) #NOTE OK
all(SOC_ATB.F$attfun(A)>0)
all(INT_ATB.F$attfun(A)>0)

## full graph out
## plotter(SOC_ATB)
## plotter(INT_ATB)
## full graph out
# DiagrammeR::export_graph(ToDiagrammeRGraph(SOC_ATB),
#              file_name=here('plots/SOC_ATB.pdf'))
# 
# DiagrammeR::export_graph(ToDiagrammeRGraph(ltbi_pathway),
#                          file_name=here('plots/ltbi_pathway.pdf'))
# 
# DiagrammeR::export_graph(ToDiagrammeRGraph(new_people_in_PPDs),
#                          file_name=here('plots/new_people_in_PPDs.pdf'))
