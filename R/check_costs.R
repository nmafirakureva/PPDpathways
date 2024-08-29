# # shorteds ATT pathway
# TB screening
summary(D[,.(int.prop.tb.sympt.screen)])
(p_scr <- D |> group_by(tb) |> summarise(p_scr = round(mean(int.prop.tb.sympt.screen),2)))
(uc_scr <- D |> group_by(tb) |> summarise(uc_scr = round(mean((cost.tb.sympt.screen + cost.overheads)))))
(c_scr <- D |> group_by(tb) |> 
    summarise(c_scr = round(mean(int.prop.tb.sympt.screen*
                             (cost.tb.sympt.screen + cost.overheads)))))
# No previous TB diagnosis
summary(D[,.(int.prop.prev.tb.dx)]) # 0, just checking it's off
# TB symptoms
summary(D[,.(int.prop.prev.tb.tx.symp)]) # doesn't matter since int.prop.prev.tb.dx is off
(p_sym <- D |> group_by(tb) |> summarise(p_sym = round(mean(int.prop.no.prev.tb.dx.symp),2)))
(uc_sym <- D |> group_by(tb) |> summarise(uc_sym = round(0*mean((cost.tb.sympt.screen + cost.overheads)))))
(c_sym <- D |> group_by(tb) |> 
    summarise(c_sym = 0*mean(int.prop.tb.sympt.screen*
                               int.prop.no.prev.tb.dx.symp*
                               (cost.tb.sympt.screen + cost.overheads))))
# Prison GP assessment
(p_gp.assess <- D |> group_by(tb) |> summarise(p_gp.assess = mean(int.prop.no.prev.tb.dx.symp.gp.assess)))

summary(D[,.(pIsolation*cost.prison.cell.isolation)]) # 0, just checking it's off
summary(D[,.(cost.prison.gp.assess)]) 
summary(D[,.(int.prop.xray*cost.chest.xray)]) 
summary(D[,.(cost.prison.gp.assess+int.prop.xray*cost.chest.xray)]) 
(uc_gp.assess <- D |> group_by(tb) |> 
    summarise(uc_gp.assess = mean((pIsolation*cost.prison.cell.isolation + 
                                    cost.prison.gp.assess + 
                                    int.prop.xray*cost.chest.xray))))  # mainly gp assessment
(c_gp.assess <- D |> group_by(tb) |> 
    summarise(c_gp.assess = mean(int.prop.tb.sympt.screen*
                                   int.prop.no.prev.tb.dx.symp*
                                   int.prop.no.prev.tb.dx.symp.gp.assess*
                                   (pIsolation*cost.prison.cell.isolation +
                                      cost.prison.gp.assess +
                                      int.prop.xray*cost.chest.xray))))  # mainly gp assessment
# Clinical suspicion of TB
(p_clin.susp <- D |> group_by(tb) |> summarise(p_clin.susp = mean(int.prop.no.prev.tb.dx.symp.tb.suspicion)))
(uc_clin.susp <- D |> group_by(tb) |> summarise(uc_clin.susp = 0*mean((cost.prison.gp.assess + 
                                                                   int.prop.xray*cost.chest.xray))))
(c_clin.susp <- D |> group_by(tb) |> 
    summarise(c_clin.susp = 0*mean(int.prop.tb.sympt.screen*
                                     int.prop.no.prev.tb.dx.symp*
                                     int.prop.no.prev.tb.dx.symp.gp.assess*
                                     int.prop.no.prev.tb.dx.symp.tb.suspicion*
                                     (cost.prison.gp.assess + int.prop.xray*cost.chest.xray))))
                                                                  
# Urgent referral to local NHS TB service - this is set to 1
# Attended
(p_attend <- D |> group_by(tb) |> summarise(p_attend = mean(attend.nhs.referral)))
summary(D[,.(cost.nhs.tb.service)]) 
summary(D[,.(cost.prison.escort)])
summary(D[,.(cost.tb.investigations)])
(uc_attend <- D |> group_by(tb) |> summarise(uc_attend = mean((cost.nhs.tb.service + cost.prison.escort + cost.tb.investigations))))
(c_attend <- D |> group_by(tb) |> 
    summarise(c_attend = mean(int.prop.tb.sympt.screen*
                                int.prop.no.prev.tb.dx.symp*
                                int.prop.no.prev.tb.dx.symp.gp.assess*
                                int.prop.no.prev.tb.dx.symp.tb.suspicion*
                                attend.nhs.referral*
                                (cost.nhs.tb.service + cost.prison.escort + cost.tb.investigations))))
# TB 
(p_tbdx <- D |> group_by(tb) |> summarise(p_tbdx = mean(int.prop.no.prev.tb.dx.symp.tb.dx)))

# actually cost.contact.tracing is per person so underestimating the cost here
(uc_tbdx <- D |> group_by(tb) |> summarise(uc_tbdx = mean(1*
                                                          (pIsolation*cost.prison.cell.isolation + nContacts*cost.contact.tracing))))
(c_tbdx <- D |> group_by(tb) |> 
    summarise(c_tbdx = mean(int.prop.tb.sympt.screen*
                              int.prop.no.prev.tb.dx.symp*
                              int.prop.no.prev.tb.dx.symp.gp.assess*
                              int.prop.no.prev.tb.dx.symp.tb.suspicion*
                              attend.nhs.referral*
                              int.prop.no.prev.tb.dx.symp.tb.dx*
                              (pIsolation*cost.prison.cell.isolation + nContacts*cost.contact.tracing))))
# ATT initiation
(p_att_init <- D |> group_by(tb) |> summarise(p_att_init = mean(int.prop.starting.att)))
(uc_att_ini <- D |> group_by(tb) |> summarise(uc_att_init = mean(1 * 
                                                                 (cost.nhs.tb.service + cost.prison.escort))))
(c_att_ini <- D |> group_by(tb) |> 
    summarise(c_att_init = mean(int.prop.tb.sympt.screen*
                                  int.prop.no.prev.tb.dx.symp*
                                  int.prop.no.prev.tb.dx.symp.gp.assess*
                                  int.prop.no.prev.tb.dx.symp.tb.suspicion*
                                  attend.nhs.referral*
                                  int.prop.no.prev.tb.dx.symp.tb.dx*
                                  int.prop.starting.att * 
                                  (cost.nhs.tb.service + cost.prison.escort))))
# completed ATT
(p_att_comp <- D |> group_by(tb) |> summarise(p_att_comp = mean(int.prop.completing.att)))
(uc_att_comp <- D |> group_by(tb) |> 
    summarise(uc_att_comp = mean(int.prop.completing.att * 
                                  (
                                    (pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) +            # DSTB OPD visits
                                      (1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visit + cost.prison.escort)) +    # MDRTB OPD visits
                                     (pDSTB*DurDSTB*cost.dsatt.drugs + (1-pDSTB)*DurMDRTB*cost.mdratt.drugs) + # Drugs
                                     cost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) +                          # DOTS
                                      cost.inpatient                                                           # Inpatient
                                    ))))
(c_att_comp <- D |> group_by(tb) |> 
    summarise(c_att_comp = mean(int.prop.tb.sympt.screen*
                                  int.prop.no.prev.tb.dx.symp*
                                  int.prop.no.prev.tb.dx.symp.gp.assess*
                                  int.prop.no.prev.tb.dx.symp.tb.suspicion*
                                  attend.nhs.referral*
                                  int.prop.no.prev.tb.dx.symp.tb.dx*
                                  int.prop.starting.att * 
                                  int.prop.completing.att * 
                                  (
                                    (pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) +            # DSTB OPD visits
                                       (1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visit + cost.prison.escort)) +    # MDRTB OPD visits
                                      (pDSTB*DurDSTB*cost.dsatt.drugs + (1-pDSTB)*DurMDRTB*cost.mdratt.drugs) + # Drugs
                                      cost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) +                          # DOTS
                                      cost.inpatient                                                           # Inpatient
                                  ))))
# did not complete ATT
(p_att_incomp <- D |> group_by(tb) |> summarise(p_att_incomp = mean(1-int.prop.completing.att)))
(uc_att_incomp <- D |> group_by(tb) |> 
    summarise(uc_att_incomp = mean((1-int.prop.completing.att) * 
                                  (
                                    IncompDurDSTB/DurDSTB*(pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + 
                                                             pDSTB*DurDSTB*(cost.dsatt.drugs + cost.dots)) + 
                                      IncompDurMDRTB/DurMDRTB*((1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visits + cost.prison.escort)  + 
                                                                 (1-pDSTB)*DurMDRTB*(cost.mdratt.drugs + cost.dots)) + 
                                      cost.inpatient                                                          # Inpatient
                                  ))))

(c_att_incomp <- D |> group_by(tb) |> 
    summarise(c_att_incomp = mean(int.prop.tb.sympt.screen*
                                    int.prop.no.prev.tb.dx.symp*
                                    int.prop.no.prev.tb.dx.symp.gp.assess*
                                    int.prop.no.prev.tb.dx.symp.tb.suspicion*
                                    attend.nhs.referral*
                                    int.prop.no.prev.tb.dx.symp.tb.dx*
                                    int.prop.starting.att *
                                    (1-int.prop.completing.att) * 
                                    (
                                      IncompDurDSTB/DurDSTB*(pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + 
                                                               pDSTB*DurDSTB*(cost.dsatt.drugs + cost.dots)) + 
                                        IncompDurMDRTB/DurMDRTB*((1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visits + cost.prison.escort)  + 
                                                                   (1-pDSTB)*DurMDRTB*(cost.mdratt.drugs + cost.dots)) + 
                                        cost.inpatient                                                          # Inpatient
                                    )
                                  )))
# check
ps <- reduce(list(p_scr, p_sym, p_gp.assess, p_clin.susp, p_attend, p_tbdx, p_att_init, p_att_comp, p_att_incomp), full_join, by = "tb")
ucosts <- reduce(list(uc_scr, uc_sym , uc_gp.assess, uc_clin.susp, uc_attend, uc_tbdx, uc_att_ini, uc_att_comp, uc_att_incomp), full_join, by = "tb")
costs <- reduce(list(c_scr, c_sym, c_gp.assess, c_clin.susp, c_attend, c_tbdx, c_att_ini, c_att_comp, c_att_incomp), full_join, by = "tb")

ps
ucosts
costs1 <- ucosts[, -1] * ps[, -1]
(costs1 <- cbind(tb = ps$tb, costs1))
costs1 <- costs1 |> rowwise() |> mutate(payoff = sum(c_across(uc_scr:uc_att_comp)))
ucosts <- ucosts |> rowwise() |> mutate(payoff = sum(c_across(uc_scr:uc_att_incomp)))
ps <- ps |> rowwise() |> mutate(payoff = prod(c_across(p_scr:p_att_init)))
costs <- costs |> rowwise() |> mutate(payoff = sum(c_across(c_scr:c_att_incomp)))


ps
ucosts
costs1
costs
costs$payoff/ps$payoff
DRS[outcome=='att'] 


ucosts <- ucosts |> mutate(metric = 'unit_cost')
ps <- ps |> rowwise() |> mutate(metric = 'probability')
costs <- costs |> rowwise() |> mutate(metric = 'cost')
payoff_result <- round(costs$payoff / ps$payoff)
costs <- ps |> select(tb, payoff) |> rename(po=payoff)|> left_join(costs) |> mutate(payoff = payoff/po) |> select(-po)

names(ps) <- names(ucosts) <- names(costs) <- gsub('p_', '', names(ps))
xr <- rbind(ps, ucosts, costs)
xr <- data.frame(xr)
xr <- xr |> 
  select(all_of(c('tb', 'metric', 'scr', 'sym', 'gp.assess', 'clin.susp', 'attend', 'tbdx', 'att_init', 'att_comp', 'att_incomp', 'payoff')))

# Function to round values based on the metric
round_by_metric <- function(metric, value) {
  if (metric == "probability") {
    return(round(value, 2))
  } else if (metric == "unit_cost" || metric == "cost") {
    return(round(value, 0))
  } else {
    return(value)
  }
}

xr$payoff <- as.numeric(xr$payoff)
# Apply the rounding function
xr <- xr %>%
  rowwise() %>%
  mutate(across(
    where(is.numeric), 
    ~ round_by_metric(metric, .), 
    .names = "rounded_{col}"
  )) %>%
  ungroup() %>%
  select(-contains("rounded_"))

# Print the rounded data frame
print(xr)

xr |> 
  filter(tb=='TBD')
DRS[outcome=='att'] 

# # shorteds ATT pathway
# TB screening
summary(D[,.(soc.prop.tb.sympt.screen)])
(p_scr <- D |> group_by(tb) |> summarise(p_scr = round(mean(soc.prop.tb.sympt.screen),2)))
(uc_scr <- D |> group_by(tb) |> summarise(uc_scr = round(mean((cost.tb.sympt.screen + cost.overheads)))))
(c_scr <- D |> group_by(tb) |> 
    summarise(c_scr = round(mean(soc.prop.tb.sympt.screen*
                                   (cost.tb.sympt.screen + cost.overheads)))))
# No previous TB diagnosis
summary(D[,.(soc.prop.prev.tb.dx)]) # 0, just checking it's off
# TB symptoms
summary(D[,.(soc.prop.prev.tb.tx.symp)]) # doesn't matter since soc.prop.prev.tb.dx is off
(p_sym <- D |> group_by(tb) |> summarise(p_sym = round(mean(soc.prop.no.prev.tb.dx.symp),2)))
(uc_sym <- D |> group_by(tb) |> summarise(uc_sym = round(0*mean((cost.tb.sympt.screen + cost.overheads)))))
(c_sym <- D |> group_by(tb) |> 
    summarise(c_sym = 0*mean(soc.prop.tb.sympt.screen*
                               soc.prop.no.prev.tb.dx.symp*
                               (cost.tb.sympt.screen + cost.overheads))))
# Prison GP assessment
(p_gp.assess <- D |> group_by(tb) |> summarise(p_gp.assess = mean(soc.prop.no.prev.tb.dx.symp.gp.assess)))

summary(D[,.(pIsolation*cost.prison.cell.isolation)]) # 0, just checking it's off
summary(D[,.(cost.prison.gp.assess)]) 
summary(D[,.(soc.prop.xray*cost.chest.xray)]) 
summary(D[,.(cost.prison.gp.assess+soc.prop.xray*cost.chest.xray)]) 
(uc_gp.assess <- D |> group_by(tb) |> 
    summarise(uc_gp.assess = mean((pIsolation*cost.prison.cell.isolation + 
                                     cost.prison.gp.assess + 
                                     soc.prop.xray*cost.chest.xray))))  # mainly gp assessment
(c_gp.assess <- D |> group_by(tb) |> 
    summarise(c_gp.assess = mean(soc.prop.tb.sympt.screen*
                                   soc.prop.no.prev.tb.dx.symp*
                                   soc.prop.no.prev.tb.dx.symp.gp.assess*
                                   (pIsolation*cost.prison.cell.isolation +
                                      cost.prison.gp.assess +
                                      soc.prop.xray*cost.chest.xray))))  # mainly gp assessment
# Clinical suspicion of TB
(p_clin.susp <- D |> group_by(tb) |> summarise(p_clin.susp = mean(soc.prop.no.prev.tb.dx.symp.tb.suspicion)))
(uc_clin.susp <- D |> group_by(tb) |> summarise(uc_clin.susp = 0*mean((cost.prison.gp.assess + 
                                                                         soc.prop.xray*cost.chest.xray))))
(c_clin.susp <- D |> group_by(tb) |> 
    summarise(c_clin.susp = 0*mean(soc.prop.tb.sympt.screen*
                                     soc.prop.no.prev.tb.dx.symp*
                                     soc.prop.no.prev.tb.dx.symp.gp.assess*
                                     soc.prop.no.prev.tb.dx.symp.tb.suspicion*
                                     (cost.prison.gp.assess + soc.prop.xray*cost.chest.xray))))

# Urgent referral to local NHS TB service - this is set to 1
# Attended
(p_attend <- D |> group_by(tb) |> summarise(p_attend = mean(attend.nhs.referral)))
summary(D[,.(cost.nhs.tb.service)]) 
summary(D[,.(cost.prison.escort)])
summary(D[,.(cost.tb.investigations)])
(uc_attend <- D |> group_by(tb) |> summarise(uc_attend = mean((cost.nhs.tb.service + cost.prison.escort + cost.tb.investigations))))
(c_attend <- D |> group_by(tb) |> 
    summarise(c_attend = mean(soc.prop.tb.sympt.screen*
                                soc.prop.no.prev.tb.dx.symp*
                                soc.prop.no.prev.tb.dx.symp.gp.assess*
                                soc.prop.no.prev.tb.dx.symp.tb.suspicion*
                                attend.nhs.referral*
                                (cost.nhs.tb.service + cost.prison.escort + cost.tb.investigations))))
# TB 
(p_tbdx <- D |> group_by(tb) |> summarise(p_tbdx = mean(soc.prop.no.prev.tb.dx.symp.tb.dx)))

# actually cost.contact.tracing is per person so underestimating the cost here
(uc_tbdx <- D |> group_by(tb) |> summarise(uc_tbdx = mean(1*
                                                            (pIsolation*cost.prison.cell.isolation + nContacts*cost.contact.tracing))))
(c_tbdx <- D |> group_by(tb) |> 
    summarise(c_tbdx = mean(soc.prop.tb.sympt.screen*
                              soc.prop.no.prev.tb.dx.symp*
                              soc.prop.no.prev.tb.dx.symp.gp.assess*
                              soc.prop.no.prev.tb.dx.symp.tb.suspicion*
                              attend.nhs.referral*
                              soc.prop.no.prev.tb.dx.symp.tb.dx*
                              (pIsolation*cost.prison.cell.isolation + nContacts*cost.contact.tracing))))
# ATT initiation
(p_att_init <- D |> group_by(tb) |> summarise(p_att_init = mean(soc.prop.starting.att)))
(uc_att_ini <- D |> group_by(tb) |> summarise(uc_att_init = mean(1 * 
                                                                   (cost.nhs.tb.service + cost.prison.escort))))
(c_att_ini <- D |> group_by(tb) |> 
    summarise(c_att_init = mean(soc.prop.tb.sympt.screen*
                                  soc.prop.no.prev.tb.dx.symp*
                                  soc.prop.no.prev.tb.dx.symp.gp.assess*
                                  soc.prop.no.prev.tb.dx.symp.tb.suspicion*
                                  attend.nhs.referral*
                                  soc.prop.no.prev.tb.dx.symp.tb.dx*
                                  soc.prop.starting.att * 
                                  (cost.nhs.tb.service + cost.prison.escort))))
# completed ATT
(p_att_comp <- D |> group_by(tb) |> summarise(p_att_comp = mean(soc.prop.completing.att)))
(uc_att_comp <- D |> group_by(tb) |> 
    summarise(uc_att_comp = mean(soc.prop.completing.att * 
                                   (
                                     (pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) +            # DSTB OPD visits
                                        (1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visit + cost.prison.escort)) +    # MDRTB OPD visits
                                       (pDSTB*DurDSTB*cost.dsatt.drugs + (1-pDSTB)*DurMDRTB*cost.mdratt.drugs) + # Drugs
                                       cost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) +                          # DOTS
                                       cost.inpatient                                                           # Inpatient
                                   ))))
(c_att_comp <- D |> group_by(tb) |> 
    summarise(c_att_comp = mean(soc.prop.tb.sympt.screen*
                                  soc.prop.no.prev.tb.dx.symp*
                                  soc.prop.no.prev.tb.dx.symp.gp.assess*
                                  soc.prop.no.prev.tb.dx.symp.tb.suspicion*
                                  attend.nhs.referral*
                                  soc.prop.no.prev.tb.dx.symp.tb.dx*
                                  soc.prop.starting.att * 
                                  soc.prop.completing.att * 
                                  (
                                    (pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) +            # DSTB OPD visits
                                       (1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visit + cost.prison.escort)) +    # MDRTB OPD visits
                                      (pDSTB*DurDSTB*cost.dsatt.drugs + (1-pDSTB)*DurMDRTB*cost.mdratt.drugs) + # Drugs
                                      cost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) +                          # DOTS
                                      cost.inpatient                                                           # Inpatient
                                  ))))
# did not complete ATT
(p_att_incomp <- D |> group_by(tb) |> summarise(p_att_incomp = mean(1-soc.prop.completing.att)))
(uc_att_incomp <- D |> group_by(tb) |> 
    summarise(uc_att_incomp = mean((1-soc.prop.completing.att) * 
                                     (
                                       IncompDurDSTB/DurDSTB*(pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + 
                                                                pDSTB*DurDSTB*(cost.dsatt.drugs + cost.dots)) + 
                                         IncompDurMDRTB/DurMDRTB*((1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visits + cost.prison.escort)  + 
                                                                    (1-pDSTB)*DurMDRTB*(cost.mdratt.drugs + cost.dots)) + 
                                         cost.inpatient                                                          # Inpatient
                                     ))))

(c_att_incomp <- D |> group_by(tb) |> 
    summarise(c_att_incomp = mean(soc.prop.tb.sympt.screen*
                                    soc.prop.no.prev.tb.dx.symp*
                                    soc.prop.no.prev.tb.dx.symp.gp.assess*
                                    soc.prop.no.prev.tb.dx.symp.tb.suspicion*
                                    attend.nhs.referral*
                                    soc.prop.no.prev.tb.dx.symp.tb.dx*
                                    soc.prop.starting.att *
                                    (1-soc.prop.completing.att) * 
                                    (
                                      IncompDurDSTB/DurDSTB*(pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + 
                                                               pDSTB*DurDSTB*(cost.dsatt.drugs + cost.dots)) + 
                                        IncompDurMDRTB/DurMDRTB*((1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visits + cost.prison.escort)  + 
                                                                   (1-pDSTB)*DurMDRTB*(cost.mdratt.drugs + cost.dots)) + 
                                        cost.inpatient                                                          # Inpatient
                                    )
    )))
# check
ps_soc <- reduce(list(p_scr, p_sym, p_gp.assess, p_clin.susp, p_attend, p_tbdx, p_att_init, p_att_comp, p_att_incomp), full_join, by = "tb")
ucosts <- reduce(list(uc_scr, uc_sym , uc_gp.assess, uc_clin.susp, uc_attend, uc_tbdx, uc_att_ini, uc_att_comp, uc_att_incomp), full_join, by = "tb")
costs <- reduce(list(c_scr, c_sym, c_gp.assess, c_clin.susp, c_attend, c_tbdx, c_att_ini, c_att_comp, c_att_incomp), full_join, by = "tb")

ps_soc
ucosts
costs1 <- ucosts[, -1] * ps_soc[, -1]
(costs1 <- cbind(tb = ps_soc$tb, costs1))
costs1 <- costs1 |> rowwise() |> mutate(payoff = sum(c_across(uc_scr:uc_att_comp)))
ucosts <- ucosts |> rowwise() |> mutate(payoff = sum(c_across(uc_scr:uc_att_incomp)))
ps_soc <- ps_soc |> rowwise() |> mutate(payoff = prod(c_across(p_scr:p_att_init)))
costs <- costs |> rowwise() |> mutate(payoff = sum(c_across(c_scr:c_att_incomp)))

# # shorteds ATT pathway
# No previous TB diagnosis
# No TB symptoms
# Chest x-ray < 6 months
# IGRA test
# Previous chest x-ray
# Normal chest x-ray
# In prison/IRC for >3 months
# Refer for treatment
# gets TPT
# completed TPT
# did not complete TPT
# does not get TPT

# # shortened LTBI pathway

# TB screening
summary(D[,.(int.prop.tb.sympt.screen)])
(p_scr <- D |> group_by(tb) |> summarise(p_scr = mean(int.prop.tb.sympt.screen)))
(uc_scr <- D |> group_by(tb) |> summarise(uc_scr = mean((cost.tb.sympt.screen + cost.overheads))))
(c_scr <- D |> group_by(tb) |> summarise(c_scr = mean(int.prop.tb.sympt.screen*(cost.tb.sympt.screen + cost.overheads))))

# No previous TB diagnosis
(p_prev_tbdx <- D |> group_by(tb) |> summarise(p_sym = mean(int.prop.prev.tb.dx)))

# No TB symptoms
(p_sym <- D |> group_by(tb) |> summarise(p_sym = mean(int.prop.no.prev.tb.dx.symp)))

# Chest x-ray < 6 months
(p_cxr6m <- D |> group_by(tb) |> summarise(p_cxr = mean(int.prop.prev.xray)))

# Chest x-ray 
(p_cxr <- D |> group_by(tb) |> summarise(p_cxr = mean(int.prop.xray)))

# IGRA test
(p_igra <- D |> group_by(tb) |> summarise(p_igra = mean(int.prop.igra.tested)))
(uc_igra <- D |> group_by(tb) |> summarise(uc_igra = mean((cost.igra.test))))
(c_igra <- D |> group_by(tb) |> 
    summarise(c_igra = mean(int.prop.tb.sympt.screen*
                              int.prop.igra.tested*
                              (cost.igra.test))))

# IGRA positive
(p_igra_pos <- D |> group_by(tb) |> summarise(p_igra_pos = mean(int.prop.igra.test.positive)))
(uc_igra_pos <- D |> group_by(tb) |> summarise(uc_igra_pos = 0*mean((cost.igra.test))))
(c_igra_pos <- D |> group_by(tb) |> summarise(c_igra_pos = 0*mean((cost.igra.test))))

# Attended
(p_attend <- D |> group_by(tb) |> summarise(p_attend = mean(attend.nhs.referral)))
(uc_attend <- D |> group_by(tb) |> summarise(uc_attend = mean((cost.nhs.tb.service + cost.prison.escort + cost.tb.investigations))))
(c_attend <- D |> group_by(tb) |> 
    summarise(c_attend = mean(int.prop.tb.sympt.screen*
                                int.prop.igra.tested*
                                int.prop.igra.test.positive*
                                attend.nhs.referral*
                                (cost.nhs.tb.service + cost.prison.escort + cost.tb.investigations))))

# No previous chest x-ray -> No chest x-ray -> Attended -> TB 
(p_tbdx <- D |> group_by(tb) |> summarise(p_tbdx = mean(int.prop.no.prev.tb.dx.symp.tb.dx)))

# In prison/IRC for >3 months
(p_long_stay <- D |> group_by(tb) |> summarise(p_long_stay = mean(int.prop.staying.o3.months)))
(uc_long_stay <- D |> group_by(tb) |> summarise(uc_long_stay = 0*mean((cost.igra.test))))
(c_long_stay  <- D |> group_by(tb) |> summarise(c_long_stay = 0*mean((cost.igra.test))))

# Refer for treatment : this is set to 1 & potentially problematic since the cost depends on attend.nhs.referral
# TODO: remove attend.nhs.referral for now and assign cost to 'gets TPT
(p_tpt <- D |> group_by(tb) |> summarise(p_tpt = mean(int.prop.staying.o3.months.tpt)))
(uc_tpt_init <- D |> group_by(tb) |> 
    summarise(c_tpt_attend = mean(1*(cost.nhs.tb.service + cost.prison.escort))))
(c_tpt_init <- D |> group_by(tb) |> 
    summarise(c_tpt_attend = mean(int.prop.tb.sympt.screen*
                                    int.prop.igra.tested*
                                    int.prop.igra.test.positive*
                                    attend.nhs.referral*
                                    (1-int.prop.no.prev.tb.dx.symp.tb.dx)*
                                    int.prop.staying.o3.months.tpt*
                                    (cost.nhs.tb.service + cost.prison.escort))))

# completed TPT
(p_tpt_comp <- D |> group_by(tb) |> summarise(p_tpt_comp = mean(int.prop.staying.o3.months.complete.tpt)))
p_tpt_comp <- p_tpt_comp |> 
  mutate(p_tpt_comp = ifelse(tb=='TBI', 0.939, p_tpt_comp))
(uc_tpt_comp <- D |> group_by(tb) |> 
    summarise(c_tpt_comp = mean(1*
                                  (durTPT*(cost.ltbi.drugs + cost.dots) + TPT.visits*(cost.tpt.opd.visit + cost.prison.escort)))))
(c_tpt_comp <- D |> group_by(tb) |> 
    summarise(c_tpt_comp = mean(int.prop.tb.sympt.screen*
                                  int.prop.igra.tested*
                                  int.prop.igra.test.positive*
                                  attend.nhs.referral*
                                  (1-int.prop.no.prev.tb.dx.symp.tb.dx)*
                                  int.prop.staying.o3.months.tpt*
                                  0.939*
                                  (durTPT*(cost.ltbi.drugs + cost.dots) + TPT.visits*(cost.tpt.opd.visit + cost.prison.escort)))))
# did not complete TPT
(p_tpt_incomp <- D |> group_by(tb) |> summarise(p_tpt_incomp = mean(1-0.939)))
(c_tpt_incomp <- D |> group_by(tb) |> 
    summarise(c_tpt_incomp = mean(int.prop.tb.sympt.screen*
                                    int.prop.igra.tested*
                                    int.prop.igra.test.positive*
                                    attend.nhs.referral*
                                    (1-int.prop.no.prev.tb.dx.symp.tb.dx)*
                                    int.prop.staying.o3.months.tpt*
                                    (1-0.939)*
                                    (IncompDurTPT*(durTPT*(cost.ltbi.drugs + cost.dots) + 
                                                     TPT.visits*(cost.tpt.opd.visit + cost.prison.escort))))))

# there is actually a BUG here: costs should be multiplied by fraction of treatment completed not the incomplete duration
# IncompDurTPT/durTPT and not IncompDurTPT
(uc_tpt_incomp <- D |> group_by(tb) |> 
    summarise(c_tpt_incomp = mean(1*
                                    (IncompDurTPT/durTPT*(durTPT*(cost.ltbi.drugs + cost.dots) + 
                                                            TPT.visits*(cost.tpt.opd.visit + cost.prison.escort))))))
(c_tpt_incomp <- D |> group_by(tb) |> 
    summarise(c_tpt_incomp = mean(int.prop.tb.sympt.screen*
                                    int.prop.igra.tested*
                                    int.prop.igra.test.positive*
                                    attend.nhs.referral*
                                    (1-int.prop.no.prev.tb.dx.symp.tb.dx)*
                                    int.prop.staying.o3.months.tpt*
                                    (1-int.prop.staying.o3.months.complete.tpt)*
                                    (IncompDurTPT/durTPT*(durTPT*(cost.ltbi.drugs + cost.dots) + 
                                                            TPT.visits*(cost.tpt.opd.visit + cost.prison.escort))))))
tpt_ps <- reduce(list(p_scr, p_igra, p_igra_pos, p_attend, p_long_stay, p_tpt, p_tpt_comp, p_tpt_incomp), full_join, by = "tb")
tpt_uc <- reduce(list(uc_scr, uc_igra, uc_igra_pos, uc_attend, uc_long_stay, uc_tpt_init, uc_tpt_comp, uc_tpt_incomp), full_join, by = "tb")
tpt_cs <- reduce(list(c_scr, c_igra, c_igra_pos, c_attend, c_long_stay, c_tpt_init, c_tpt_comp, c_tpt_incomp), full_join, by = "tb")

tpt_ps <- tpt_ps |> rowwise() |> mutate(payoff = prod(c_across(p_scr:p_tpt_comp)))
tpt_cs <- tpt_cs |> rowwise() |> mutate(payoff = sum(c_across(c_scr:c_tpt_incomp)))
tpt_uc <- tpt_uc |> rowwise() |> mutate(payoff = '-')

tpt_uc <- tpt_uc |> mutate(metric = 'unit_cost')
tpt_ps <- tpt_ps |> mutate(metric = 'probability')
tpt_cs <- tpt_cs |> mutate(metric = 'cost')

tpt_cs <- tpt_ps |> select(tb, payoff) |> rename(po=payoff)|> left_join(tpt_cs) |> mutate(payoff = payoff/po) |> select(-po)

names(tpt_ps) <- names(tpt_cs) <- names(tpt_uc) <- gsub('p_', '', names(tpt_ps))
tpxr <- rbind(tpt_ps, tpt_uc, tpt_cs)
tpxr <- data.frame(tpxr)
tpxr <- tpxr |> 
  select(all_of(c('tb', 'metric', 'scr', 'igra', 'igra_pos', 'long_stay', 'tpt', 'tpt_comp', 'tpt_incomp', 'payoff')))

# Function to round values based on the metric
round_by_metric <- function(metric, value) {
  if (metric == "probability") {
    return(round(value, 2))
  } else if (metric == "unit_cost" || metric == "cost") {
    return(round(value, 0))
  } else {
    return(value)
  }
}

tpxr$payoff <- as.numeric(tpxr$payoff)
# Apply the rounding function
tpxr <- tpxr %>%
  rowwise() %>%
  mutate(across(
    where(is.numeric), 
    ~ round_by_metric(metric, .), 
    .names = "rounded_{col}"
  )) %>%
  ungroup() %>%
  select(-contains("rounded_"))

# Print the rounded data frame
print(tpxr)

tpxr |> 
  filter(tb=='TBI')
DRS[outcome=='tpt'] 

# Previous chest x-ray
# Normal chest x-ray
# In prison/IRC for >3 months
# Refer for treatment
# gets TPT
# completed TPT
# did not complete TPT
# does not get TPT