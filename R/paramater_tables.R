## prior parameters
PD <- read.csv(here('indata/ProbParms2.csv')) #read in probability parameters
AD <- read.csv(here('indata/DiagnosticAccuracy.csv')) #read in accuracy parameters
RD <- fread(gh('indata/RUParms.csv'))    #read resource use data
CD <- fread(gh('indata/CostParms.csv'))    #read cost data

names(PD)
names(RD)
names(CD)

names <- c("ParameterName", "Mean", "Range","Description","Source","SourceFull")
names(PD) <- names(RD) <- names(CD) <- names

# PD1 <- PD |> 
#   filter(ParameterName != 'tb.prev') |>
#   mutate(ParameterName = paste0('soc.', ParameterName))

PD2 <- PD |>
  filter(ParameterName != 'tb.prev') 

PD3 <- PD |> 
  filter(ParameterName %in% c('tb.prev', 'ltbi.prev', 'prog.tb', 'progInf', 'progNInf')) 

PD0 <- rbind(PD2, PD3, RD)

PD0$NAME <- PD0$ParameterName

PD0 <- PD0 |> 
  mutate(Mean = as.numeric(Mean), Range = as.character(Range),
         Median = ifelse(Range!='', paste0(Mean, ' (', Range, ')'), Range),
         Median = ifelse(!is.na(as.numeric(Range)), '', Median)) 

PD1 <- PD0 |> 
  filter(Median != '') # parameters to be sampled

PD2 <- PD0 |> 
  filter(Median == '') # parameters to be fixed (no uncertainty estimates)

tmp <- PD1 %>%
  extract(Median, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") %>%
  separate(lo,c("lo","hi"),"-") %>%
  mutate_at(c("mid", "lo","hi"), as.numeric)

tmp <- setDT(tmp)
tmp1 <- getLNparms(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2,med=FALSE)
tmp[,DISTRIBUTION:=paste0("LN(",round(tmp1$mu,3),",",round(tmp1$sig,3),")")] #LN distributions

# tmp1 <- getAB(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2)
# tmp[,DISTRIBUTION:=paste0("B(",tmp1$a,",",tmp1$b,")")] # Beta distributions
# tmp[,DISTRIBUTION:=paste0("G(",tmp[,mid]^2/((tmp[,hi]-tmp[,lo])^2/3.92^2),",",((tmp[,hi]-tmp[,lo])^2/3.92^2)/tmp[,mid],")")] # Gamma distributions
PD1 <- PD1 |> 
  mutate(NAME=ParameterName,
         DISTRIBUTION=tmp$DISTRIBUTION,
         # DISTRIBUTION=ifelse(DISTRIBUTION=='LN(NA,NA)', NA, DISTRIBUTION))
         DISTRIBUTION=ifelse(DISTRIBUTION=='B(NA,NA)', NA, DISTRIBUTION))

## 
PD1 <- PD1 |> 
  filter(DISTRIBUTION!="")

# Diagnostic accuracy 
tmp <- AD %>%
  extract(mqrng, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") %>%
  separate(lo,c("lo","hi"),"to") %>%
  mutate_at(c("mid", "lo","hi"), as.numeric)

tmp <- setDT(tmp)
tmp1 <- getAB(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2)
tmp[,DISTRIBUTION:=paste0("B(",round(tmp1$a),",",round(tmp1$b),")")] # Beta distributions

ParmsTab1 <- PD1 |> 
  select(NAME, Description, DISTRIBUTION, Median, Source) |> 
  mutate(NAME = gsub('prop.', '', NAME),
         NAME = gsub('prev.tb.dx.|prev.tb.dx.', 'prev.tb.', NAME),
         NAME = case_when(
           grepl('symp$', NAME) ~ 'tb.symptoms',
           grepl('gp.assess', NAME) ~ 'gp.assessment',
           grepl('tb.suspicion', NAME) ~ 'clinical.tb.suspicion',
           grepl('attend.referral|nhs.referral', NAME) ~ 'attend.nhs.referral',
           grepl('starting.att', NAME) ~ 'att.initiation',
           .default = NAME
         ))

unique(ParmsTab1$NAME)

keep <- c('prev.tb.started.att', 'prev.tb.still.on.att', 'prev.tb.continue.att', 
          'tb.symptoms', 'gp.assessment', 'clinical.tb.suspicion', 'xray',
          'attend.nhs.referral', 'att.initiation',  
          'igra.tested', 'igra.test.positive', 'ltbi.prev', 
          'verbal_screen_time', 'nContacts')
ParmsTab1 <- ParmsTab1 |> 
  filter(NAME %in% keep) |>
  mutate(NAME = factor(NAME, levels = keep),
         Description = case_when(
           # NAME == 'prev.tb.started.att' ~ 'prev.tb.started.att',
           # NAME == 'prev.tb.still.on.att' ~ 'prev.tb.still.on.att',
           NAME == 'prev.tb.continue.att' ~ 'Past TB diagnosis, continue anti-TB treatment',
           NAME == 'tb.symptoms' ~ 'TB symptoms',
           NAME == 'gp.assessment' ~ 'Prison GP assessment',
           NAME == 'clinical.tb.suspicion' ~ 'Clinical TB suspicion',
           # NAME == 'xray' ~ 'xray',
           NAME == 'attend.nhs.referral' ~ 'Attend NHS referral',
           NAME == 'att.initiation' ~ 'Starting anti-TB treatment (ATT)',
           .default = Description
         )) |>
  distinct(NAME, .keep_all = TRUE) |> 
  arrange(NAME) |> 
  rename_with(~ c('Name', 'Description', 'Distribution', 'Estimate (Range)', 'Source'), 
              c('NAME', 'Description', 'DISTRIBUTION','Median', 'Source'))

unique(ParmsTab1$Name)

ParmsTab1 <- ParmsTab1 |> 
  mutate(Distribution = case_when(
    Name == 'attend.nhs.referral' ~ 'B(7.6673,2.3715)',
    Name == 'nContacts' ~ 'LN(1.84,0.1)',
    Name == 'att.initiation' ~ 'B(8,2)',
    Name == 'gp.assessment' ~ 'B(30,2)',
    Name == 'prev.tb.continue.att' ~ 'B(10,2)',
    Name == 'prev.tb.started.att' ~ 'B(10,2)',
    Name == 'prev.tb.still.on.att' ~ 'B(10,2)',
    Name == 'xray' ~ 'B(19.7,5.87)',
    .default = Distribution))

ParmsTab2 <- tmp |>  
  mutate(Median = paste0(mid, ' (', lo, '-', hi)) |>
  select(NAME, DESCRIPTION, DISTRIBUTION, Median, SOURCE) |> 
  rename_with(~ c('Name', 'Description', 'Distribution', 'Estimate (Range)', 'Source'), 
              c('NAME', 'DESCRIPTION', 'DISTRIBUTION','Median', 'SOURCE'))

# Fixed parameters
names(PD2)

ParmsTab3 <- PD2 |> 
  select(NAME, Description, Mean, Source) |> 
  mutate(NAME = gsub('prop.', '', NAME),
         NAME = gsub('prev.tb.dx.|prev.tb.dx.', 'prev.tb.', NAME),
         NAME = case_when(
           grepl('attend.referral$', NAME) ~ 'attend.nhs.referral',
           grepl('.tp$', NAME) ~ 'tpt.initiation',
           grepl('ompleting.att', NAME) ~ 'att.completion',
           grepl('complete.tpt', NAME) ~ 'tpt.completion',
           .default = NAME
         ))

unique(ParmsTab3$NAME)

keep <- c('prev.tb.dx', 'prev.xray', 'abnormal.xray', 'pIsolation', 
          'tpt.initiation', 'tpt.completion', 'att.completion',
          'pDSTB', 'dstb.visits', 'mdrtb.visits', 'TPT.visits', 
          'DurDSTB', 'DurMDRTB', 'DurTPT', 'IncompDurDSTB', 'IncompDurMDRTB', 'IncompDurTPT',
          'smear.positive', 'DurDSTBIsolation', 'DurMDRT2BIsolation')
ParmsTab3 <- ParmsTab3 |> 
  filter(NAME %in% keep) |>
  mutate(NAME = factor(NAME, levels = keep),
         Distribution = 'Fixed',
         Description = case_when(
           grepl('attend.nhs.referral', NAME) ~ 'Attending NHS referral',
           grepl('tpt.initiation', NAME) ~ 'Initiating TPT',
           grepl('tpt.completion', NAME) ~ 'Completing TPT',
           grepl('att.completion', NAME) ~ 'Completing ATT',
           grepl('DurDSTB', NAME) ~ 'Duration of DSTB treatment',
           grepl('DurMDRTB', NAME) ~ 'Duration of MDR-TB treatment',
           grepl('DurTPT', NAME) ~ 'Duration of TPT',
           grepl('IncompDurDSTB', NAME) ~ 'Duration of incomplete DSTB treatment',
           grepl('IncompDurMDRTB', NAME) ~ 'Duration of incomplete MDR-TB treatment',
           grepl('IncompDurTPT', NAME) ~ 'Duration of incomplete TPT',
           grepl('smear.positive', NAME) ~ 'Proportion of smear positive patients',
           grepl('DurDSTBIsolation', NAME) ~ 'Duration of DSTB isolation',
           grepl('DurMDRT2BIsolation', NAME) ~ 'Duration of MDR-TB isolation',
           .default = Description
         )) |> 
  distinct(NAME, .keep_all = TRUE) |> 
  arrange(NAME) |> 
  rename_with(~ c('Name', 'Description', 'Distribution', 'Estimate (Range)', 'Source'), 
              c('NAME', 'Description', 'Distribution','Mean', 'Source')) 

ru <- c('verbal_screen_time', 'nContacts',
          'pDSTB', 'dstb.visits', 'mdrtb.visits', 'TPT.visits', 
          'DurDSTB', 'DurMDRTB', 'DurTPT', 'IncompDurDSTB', 'IncompDurMDRTB', 'IncompDurTPT',
          'smear.positive', 'DurDSTBIsolation', 'DurMDRT2BIsolation')

ProbParmsTab <- rbind(ParmsTab1, ParmsTab2, ParmsTab3) |>
  filter(!Name %in% ru) 

RUParmsTab <- rbind(ParmsTab1, ParmsTab2, ParmsTab3) |>
  filter(Name %in% ru) 

# costs 
names(rcsts)
CostsParmsTab <- rcsts |> 
  select(ParameterName, Description, cost.m, cost.sd, Source) |>
  rename_with(~ c('Name', 'Description', 'Mean', 'SD', 'Source'), 
              c('ParameterName', 'Description', 'cost.m', 'cost.sd', 'Source')) 

check <- ProbParmsTab |> 
  filter(Distribution != 'Fixed') |>
  rename(NAME=Name, DISTRIBUTION=Distribution) |> 
  select(NAME, DISTRIBUTION) |>
  parse.parmtable(outfile=here('outdata/out.csv'),
                  testdir = here('plots/test'))
                  
# save for later use in creating parameter table for Appendix
write.csv(ProbParmsTab, file = here("outdata/ProbParmsTab.csv"), row.names = FALSE)
write.csv(RUParmsTab, file = here("outdata/RUParmsTab.csv"), row.names = FALSE)
write.csv(CostsParmsTab, file = here("outdata/CostsParmsTab.csv"), row.names = FALSE)
