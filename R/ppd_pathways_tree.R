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

# notb <- txt2tree(here('indata/noTB.outcomes.txt')) # no tb
tpt <- txt2tree(here('indata/3_TPT_outcomes.txt')) # tpt

## default prob/cost:
tb$Set(p=1)
# notb$Set(p=1)
tpt$Set(p=1)
tb$Set(cost=0)
# notb$Set(cost=0)
tpt$Set(cost=0)

# TB outcomes
tbtxo <- Node$new('TB')
tbtxo$AddChildNode(tb)

# tbtx <- Node$new('TB')
# tbtxo <- tbtx$AddChild('TB outcomes')

# no tb outcomes
notb <- Node$new('Not TB')
ltbi <- notb$AddChild('LTBI pathway')

## ====== function to add outcomes & counters
AddOutcomes <- function(D){
  ## === cost and probs (defaults)
  D$Set(p=1)
  D$Set(cost=0)
  
  ## === merge 'New people in PPDs' with 'LTBI screening pathway' to create final tree ===
  # MergeByName(D,notb,'Not TB', leavesonly = TRUE) # first add a branch to 'No TB' for easy merging in the next step
  MergeByName(D,No_urgent_GP_referral,'No urgent prison GP referral',leavesonly = TRUE)
  MergeByName(D,Refer_to_TB_services,'Refer to TB services',leavesonly = TRUE)
  MergeByName(D,ltbi_pathway,'LTBI pathway',leavesonly = TRUE)
  
  ## final outcomes
  MergeByName(D,tpt,'TPT outcomes',leavesonly = TRUE)
  MergeByName(D,tbtxo,'TB',leavesonly = TRUE)
  
  ## ===========  other counters
  ## check
  D$Set(check=1)
  D$Set(check=0,filterFun=function(x) length(x$children)>0)
  
  ## TB screening
  D$Set(screen=0)
  D$Set(screen=1,filterFun=function(x) x$name=='TB screening')
  
  ## TB dx
  D$Set(prevtb=0)
  D$Set(prevtb=1,filterFun=function(x) x$name=='TB')
  
  ## ATT courses
  D$Set(att=0)
  D$Set(att=1,
        filterFun=function(x)x$name=='ATT')
  
  ## prevtb no ATT 
  D$Set(noatt=0)
  D$Set(noatt=1,
        filterFun=function(x)x$name=='no ATT')
  
  ## tested on IGRA
  D$Set(igra=0)
  D$Set(igra=1,filterFun=function(x) x$name=='IGRA test')
  
  ## positive IGRA
  D$Set(ltbi=0)
  D$Set(ltbi=1,filterFun=function(x) x$name=='IGRA test +')
  
  ## negative IGRA
  D$Set(noltbi=0)
  D$Set(noltbi=1,filterFun=function(x) x$name=='IGRA test -')
  
  ## TPT courses
  D$Set(tpt=0)
  D$Set(tpt=1,
        filterFun=function(x)x$name=='gets TPT')
  
  ## positive IGRA no TPT
  D$Set(ltbinotpt=0)
  D$Set(ltbinotpt=1,filterFun=function(x) x$name=='does not get TPT')
  
  ## negative IGRA no TPT
  D$Set(noltbinotpt=0)
  D$Set(noltbinotpt=1,filterFun=function(x) x$name=='IGRA test -')
  
  ## attend
  D$Set(attend=0)
  D$Set(attend=1,filterFun=function(x) x$name=='completed ATT')
  D$Set(attend=1,filterFun=function(x) x$name=='did not complete ATT')
  
  ## tptend
  D$Set(tptend=0)
  D$Set(tptend=1,filterFun=function(x) x$name=='completed TPT')
  D$Set(tptend=1,filterFun=function(x) x$name=='did not complete TPT')
  
  ## no TPT & no ATT
  D$Set(notxend=0)
  D$Set(notxend=1,filterFun=function(x) x$name=='no ATT')
  D$Set(notxend=1,filterFun=function(x) x$name=='does not get TPT')
  
  return(D)
}

# main tree structures
new_people_in_PPDs <- txt2tree(here('indata/1.2_new_people_in_PPDs.txt'))
Refer_to_TB_services <- txt2tree(here('indata/1.2.1_Refer_to_TB_services.txt'))
No_urgent_GP_referral <- txt2tree(here('indata/1.2.2_No_urgent_GP_referral.txt'))
ltbi_pathway <- txt2tree(here('indata/2_LTBI_pathway.txt'))

# DiagrammeR::export_graph(ToDiagrammeRGraph(new_people_in_PPDs),
#                          file_name=here('plots/new_people_in_PPDs.pdf'))
# DiagrammeR::export_graph(ToDiagrammeRGraph(ltbi_pathway),
#                          file_name=here('plots/ltbi_pathway.pdf'))

## merge in extras, make model of care branches, write out
tempTree <- AddOutcomes(new_people_in_PPDs)

## === SOC
SOC <- Node$new('Standard of care pathway')
SOC$AddChildNode(tempTree)

SOC$name <- 'Standard of care pathway'


## # defining tree quantities 
## # p = probability/proportion
## # cost = cost
## # cost.screen = screening cost
## # cost.tpt = TPT cost
## # cost.tb.assessment = TB assessment cost
## # cost.att = ATT cost
## # cost.ppd = PPD cost
## # cost.nhs = NHS cost
## # screen = screened for TB
## # prevtbdx = previous TB diagnosed
## # igra = IGRAS tested
## # ltbi = IGRAS test positive
## # noltbi = IGRAS test negative
## # tpt = TPT initiation 
## # ltbinotpt = IGRAS test positive no TPT
## # noltbinotpt = IGRAS test negative no TPT
## # coprevtb = coprevalent TB diagnosed
## # att = anti-TB treatments
## # noatt = TB no ATT
## # check = check

tree2file(SOC,filename = here('indata/CSV/SOC.csv'),
          'p','cost','cost.screen',	'cost.tpt','cost.tb.assessment','cost.att','cost.ppd', 'cost.nhs',	
          'screen', 'igra',	'ltbi', 'noltbi', 'tpt', 'ltbinotpt', 'noltbinotpt', 'prevtbdx', 'coprevtb', 'att', 'noatt',  'attend', 'tptend', 'notxend','check')

labdat <- c('p','cost','cost.screen',	'cost.tpt','cost.tb.assessment','cost.att','cost.ppd', 'cost.nhs',	
            'screen', 'igra',	'ltbi', 'noltbi', 'tpt', 'ltbinotpt', 'noltbinotpt', 'prevtbdx', 'coprevtb','att', 'noatt', 'attend', 'tptend', 'notxend', 'check')


## create version with probs/costs
fn <- here('indata/CSV/SOC1.csv')
if(file.exists(fn)){
  ## read
  labz <- fread(fn)
  labz$p <- gsub('<3','u3.',labz$p)
  labz$p <- gsub('3\\+', 'o3.', labz$p)
  # labz$p <- gsub('soc.', '', labz$p)
  # labz$cost <- gsub('<','u',labz$cost)
  LabelFromData(SOC,labz[,..labdat]) #add label data
  ## NOTE checks need redoing
  SOC$Set(check=1)
  SOC$Set(check=0,filterFun=function(x) length(x$children)>0)
  ## save out
  tree2file(SOC,filename = here('indata/CSV/SOC2.csv'),
            'p','cost.screen',	'cost.tpt','cost.tb.assessment',
            'cost.att','cost.ppd', 'cost.nhs','cost',	
            'screen', 'igra',	'ltbi', 'noltbi', 'tpt',
            'ltbinotpt', 'noltbinotpt', 'prevtbdx', 'coprevtb','att', 'noatt',  'attend', 'tptend', 'notxend', 'check')
}


## NOTE this would ideally be moved up into the workflow above
## add a notx variable = no ATT *and* no TPT
leaves <- as.integer(SOC$Get('check')) #indicator for being a leaf
sum(leaves) == SOC$leafCount
notx <- as.integer((!SOC$Get('attend')) * (!SOC$Get('tptend'))) * leaves #only 1 on leaves

SOC$Set(notx = 0)
SOC$Set(notx = notx)

## this gives us 3 outcome functions: tpt,att, notx, which are exhaustive
labz[,sum(tptend==1)] + labz[,sum(attend==1)] + sum(notx) == sum(leaves) #only 1 on leaves
# labz[,sum(tptend==1)] + labz[,sum(attend==1)] + labz[,sum(notxend==1)] == nrow(labz) #check
##  & exclusive:
labz[tptend>1 & attend > 1]
labz[,table(tptend,attend)]
labz[,table(tptend,notx)]
labz[,table(attend,notx)]

## === INT
INT <- Clone(SOC)
INT$name <- 'Intervention model'

tree2file(INT,filename = here('indata/CSV/INT.csv'),
          'p','cost','cost.screen',	'cost.tpt','cost.tb.assessment','cost.att','cost.ppd', 'cost.nhs',	
          'screen', 'igra',	'ltbi', 'noltbi', 'tpt', 'ltbinotpt', 'noltbinotpt', 'prevtbdx', 'coprevtb','att', 'noatt',  'attend', 'tptend', 'notxend', 'check')


## create version with probs/costs
fn <- here('indata/CSV/INT1.csv')
if(file.exists(fn)){
  ## read
  labz <- fread(fn)
  labz$p <- gsub('<3','u3.',labz$p)
  labz$p <- gsub('3\\+', 'o3.', labz$p)
  # labz$p <- gsub('int.', '', labz$p)
  # labz$cost <- gsub('int.', '', labz$cost)
  # labz[,cost:=q]; labz[,q:=NULL]     #rename
  LabelFromData(INT,labz[,..labdat]) #add label data
  ## NOTE checks need redoing
  INT$Set(check=1)
  INT$Set(check=0,filterFun=function(x) length(x$children)>0)
  ## save out
  tree2file(INT,filename = here('indata/CSV/INT2.csv'),
            'p','cost.screen',	'cost.tpt','cost.tb.assessment','cost.att','cost.ppd', 'cost.nhs','cost',	
            'screen', 'igra',	'ltbi', 'noltbi', 'tpt', 'ltbinotpt', 'noltbinotpt', 'prevtbdx', 'coprevtb','att',  'attend', 'tptend', 'notxend', 'noatt', 'check')
}


## NOTE this would ideally be moved up into the workflow above
## add a notx variable = no ATT *and* no TPT
leaves <- as.integer(INT$Get("check")) # indicator for being a leaf
sum(leaves) == INT$leafCount

notx <- as.integer((!INT$Get('attend')) * (!INT$Get('tptend'))) * leaves
INT$Set(notx = 0)
INT$Set(notx = notx)
## this gives us 3 outcome functions: tpt,att, notx, which are exhaustive
labz[,sum(tptend==1)] + labz[,sum(attend==1)] + sum(notx) == sum(leaves) #check
##  & exclusive:
labz[tptend>1 & attend > 1]
labz[,table(tptend,attend)]
labz[,table(tptend,notx)]
labz[,table(attend,notx)]


## make functions
fnmz <- labdat[-1]
fnmz <- c(fnmz,'notx')

## full tree
SOC.F <- makeTfuns(SOC,fnmz)
INT.F <- makeTfuns(INT,fnmz)


## NOTE making pruned trees conditioned on outcomes (subtrees ending variable > 0)
SOC.att <- PruneByOutcome(SOC,'attend')
SOC.tpt <- PruneByOutcome(SOC,'tptend')
SOC.notx <- PruneByOutcome(SOC, "notx")
INT.att <- PruneByOutcome(INT, "attend")
INT.tpt <- PruneByOutcome(INT, "tptend")
INT.notx <- PruneByOutcome(INT, "notx")

## checking...
leaves <- as.integer(SOC.notx$Get("check")) # indicator for being a leaf
sum(leaves) == SOC.notx$leafCount

notx <- as.integer(SOC.notx$Get("notx"))
attend <- as.integer(SOC.notx$Get("attend"))
tptend <- as.integer(SOC.notx$Get("tptend"))

sum(notx+attend+tptend)==sum(leaves)  #each leaf has outcome
sum((notx + attend + tptend)*!leaves) #only on leaves
which(leaves == 1 & (notx + attend + tptend) == 0) #but I see them in the CSV??

tree2file(SOC.att,
  filename = here("indata/CSV/SOC.att.csv"),
  "notx", "attend", "tptend", "check"
)
tree2file(SOC.tpt,
  filename = here("indata/CSV/SOC.tpt.csv"),
  "notx", "attend", "tptend",  "check"
)
tree2file(SOC.notx,
  filename = here("indata/CSV/SOC.notx.csv"),
  "notx", "attend", "tptend", "check"
)


## retricted trees:
SOC.att.F <- makeTfuns(SOC.att,fnmz)
SOC.tpt.F <- makeTfuns(SOC.tpt,fnmz)
SOC.notx.F <- makeTfuns(SOC.notx,fnmz)
INT.att.F <- makeTfuns(INT.att,fnmz)
INT.tpt.F <- makeTfuns(INT.tpt,fnmz)
INT.notx.F <- makeTfuns(INT.notx,fnmz)



## running all function
runallfuns <- function(D,arm='all'){
  done <- FALSE
  if('SOC' %in% arm | arm[1]=='all'){
    cat('Running functions for SOC:\n')
    for(nm in names(SOC.F)){
      snm <- gsub('fun','',nm)
      snma <- paste0(snm,'.soc')
      D[[snma]] <- SOC.F[[nm]](D)
      cat('...',snm,' run...\n')
      done <- TRUE
    }
  }
  if('INT' %in% arm | arm[1]=='all'){
    cat('Running functions for INT:\n')
    for(nm in names(INT.F)){
      snm <- gsub('fun','',nm)
      snma <- paste0(snm,'.int')
      D[[snma]] <- INT.F[[nm]](D)
      cat('...',snm,' run...\n')
      done <- TRUE
    }
  }
  
  
  if(!done)stop('Functions not run! Likely unrecognised arm supplied.')
  return(D)
}

## checking
vrz.soc <- showAllParmz(SOC)
vrz.int <- showAllParmz(INT)

vrz <- c(vrz.soc,
         vrz.int
)

vrz <- unique(vrz) #NOTE
vrz.soc[grepl('prop',vrz.soc)]
vrz.soc[!grepl('cost|prop',vrz.soc)]
vrz.soc[grepl('cost',vrz.soc)]

# # write out as csv
# write.csv(data.frame(vrz),here('indata/CSV/parmz.csv'),row.names=FALSE)
# 
# # save all as R.data
# 
# # Save an object to a file
# save.image(file = here("outdata/temp.RData"))
# 
# # Restore the object
# load(here("outdata/temp.RData"))

A <- makeTestData(5e3,vrz)


## checks
any(SOC.F$checkfun(A)!=1) #NOTE OK
any(round(SOC.F$checkfun(A))!=1) # Looks like there are some rounding errors

# INT.F$checkfun(A) #NOTE OK
all(abs(SOC.F$checkfun(A)-1)<1e-10) #NOTE OK
all(abs(INT.F$checkfun(A)-1)<1e-10) #NOTE OK

summary(SOC.F$checkfun(A)) # NOTE OK
summary(INT.F$checkfun(A)) # NOTE OK

## BUG? up to 1% off

all(SOC.F$attfun(A)>0)
all(INT.F$attfun(A)>0)

# plotter(new_people_in_PPDs)
## full graph out
## plotter(SOC)
## plotter(INT)
## full graph out
# DiagrammeR::export_graph(ToDiagrammeRGraph(SOC),
#              file_name=here('plots/SOC.pdf'))

# DiagrammeR::export_graph(ToDiagrammeRGraph(SOC.att),
#                          file_name=here('plots/SOC_att.pdf'))

# 
# DiagrammeR::export_graph(ToDiagrammeRGraph(ltbi_pathway),
#                          file_name=here('plots/ltbi_pathway.pdf'))
# 
# DiagrammeR::export_graph(ToDiagrammeRGraph(new_people_in_PPDs),
#                          file_name=here('plots/new_people_in_PPDs.pdf'))


