############################################################
# R-project                                                #
# Program      : a2-EarlyEfficacyStop.r                    #
# Protocol     :                                           #
# Date         : <2020/06/01 (Mon) 08:45>                  #
# Last         :                                           #
# Programmer   : yoshifumi ukyo                            #
#                                                          #
############################################################
# [Ver.0000]                                               #
# Memorandom   :                                           #
#                                                          #
############################################################


#----- clean-up working directory
rm(list = (ls(all = TRUE)))
#----- library assignment 
base_dir <- ""
setwd(base_dir)
set.seed(4989)



library(gsDesign)
gsDesign(k = 2, test.type = 1, alpha = 0.025, beta = 0.2,  
         n.fix = 694, timing = c(208/694, 1), sfu = sfHSD, sfupar = -4)


library(dplyr)

N_total <- 694 
delta <- seq(from = 1.6, to = 2.0, by = 0.1)  
sigma <- 7.5 

#----- wo: without overrun, wt: with overrun 
sim_summary <- data.frame(N_total = N_total, delta = delta, sigma = sigma, 
                          power = 0, ASN_wo = 0, ASN_wt = 0, IA_stop_prob = 0)

for (k in 1:nrow(sim_summary)){ 
  
  nIter <- 10000 
  sim_result <- data.frame(nIter = 1:nIter, 
                           IA_N_wo = 0, IA_N_wt = 0, IA_stop = 0, FA_sig = 0, 
                           ASN_wo = 0, ASN_wt = 0)
  
  for (i in 1:nrow(sim_result)){ 
    
    N_total <- sim_summary$N_total[k]
    delta   <- sim_summary$delta[k]
    sigma   <- sim_summary$sigma[k]
    
    enroll_time <- c(seq(from = 0, to = 86, length = 688), seq(from = 86, to = 87, length = 6))
    treat_group <- rbinom(n = N_total, size = 1, prob = 0.5)
    
    y <- rep(0, N_total)
    for (j in 1:N_total){ 
      if (treat_group[j] == 0) {
        y[j] <- rnorm(n = 1, mean = 0, sd = sigma)
      } else if (treat_group[j] == 1) {
        y[j] <- rnorm(n = 1, mean = delta, sd = sigma)
      }
    }
    
    d_sim <- data.frame(id = 1:N_total, 
                        enroll_time = enroll_time, 
                        treat_group = treat_group, 
                        y = y)
    
    IA_enroll_time <- d_sim[208, ]$enroll_time
    IA_cutoff_time <- IA_enroll_time + 26 
    
    d_sim <- d_sim %>% 
      mutate(IA_wo = if_else(enroll_time <= IA_enroll_time, 1, 0), 
             IA_wt = if_else(enroll_time <= IA_cutoff_time, 1, 0)) %>% 
      ungroup() 
    
    #----- IA analysis 
    g0_IA <- subset(d_sim, treat_group == 0 & IA_wo == 1)$y
    g1_IA <- subset(d_sim, treat_group == 1 & IA_wo == 1)$y
    IA_pval <- t.test(x = g0_IA, y = g1_IA, alternative = "two.sided", conf.level = 0.95)$p.value
    
    #----- FA analysis 
    g0_FA <- subset(d_sim, treat_group == 0)$y 
    g1_FA <- subset(d_sim, treat_group == 1)$y 
    FA_pval <- t.test(x = g0_FA, y = g1_FA, alternative = "two.sided", conf.level = 0.95)$p.value 
    
    sim_result$IA_N_wo[i] <- nrow(subset(d_sim, IA_wo == 1))
    sim_result$IA_N_wt[i] <- nrow(subset(d_sim, IA_wt == 1))
    sim_result$IA_stop[i] <- ifelse(IA_pval < 2 * 0.0011, 1, 0)
    sim_result$FA_sig[i]  <- ifelse(FA_pval < 0.05, 1, 0)
    
    sim_result <- sim_result %>% 
      mutate(ASN_wo = if_else(IA_stop == 1, IA_N_wo, N_total), 
             ASN_wt = if_else(IA_stop == 1, IA_N_wt, N_total)) %>% 
      ungroup() 
    
  }
  
  sim_summary$power[k]  <- nrow(subset(sim_result, IA_stop == 1 | FA_sig == 1)) / nrow(sim_result)
  sim_summary$ASN_wo[k] <- mean(sim_result$ASN_wo) 
  sim_summary$ASN_wt[k] <- mean(sim_result$ASN_wt)
  sim_summary$IA_stop_prob[k] <- mean(sim_result$IA_stop)
  
}




