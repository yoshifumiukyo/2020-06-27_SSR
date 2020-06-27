############################################################
# R-project                                                #
# Program      : a6-Plan4-H0.r                             #
# Protocol     :                                           #
# Date         : <2020/06/04 (Thu) 19:16>                  #
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
set.seed(1234)



# delta <- c(1.6, 1.7, 1.8, 1.9, 2.0) 
delta <- 0.0
sigma <- seq(from = 5.0, to = 10.0, by = 0.1)  
IA_N  <- 208 
FA_N  <- 442 
N_max <- 884 

sim_summary <- expand.grid(delta = delta, sigma = sigma, 
                          IA_N = IA_N, FA_N = FA_N, N_max = N_max, 
                          Pr_fav = 0, Pr_prom = 0, Pr_unfav = 0, 
                          Sig_fixed = 0, 
                          Sig_fixed_fav = 0, Sig_fixed_prom = 0, Sig_fixed_unfav = 0, 
                          Sig_ada = 0, 
                          Sig_ada_fav = 0, Sig_ada_prom = 0, Sig_ada_unfav = 0, 
                          ASN = 0, 
                          ASN_fav = 0, ASN_prom = 0, ASN_unfav = 0)

for (k in 1:nrow(sim_summary)) { 
  
  nIter <- 10000 
  d_sim <- data.frame(nIter = 1:nIter, Zone = NA, CP = 0, SS = 0, Sig_fixed = 0, Sig_ada = 0) 
  
  for (i in 1:nrow(d_sim)) { 
    
    mu_g0 <- 0 
    mu_g1 <- sim_summary$delta[k] 
    sd_g0 <- sim_summary$sigma[k] 
    sd_g1 <- sim_summary$sigma[k]
    IA_N  <- sim_summary$IA_N[k] 
    FA_N  <- sim_summary$FA_N[k] 
    N_max <- sim_summary$N_max[k] 
    
    rand_vec <- ifelse(runif(n = N_max) < 0.5, 1, 0)
    y <- rep(0, N_max) 
    for (j in 1:N_max) { 
      if (rand_vec[j] == 0) {
        y[j] <- rnorm(n = 1, mean = mu_g0, sd = sd_g0)
      } else if (rand_vec[j] == 1) {
        y[j] <- rnorm(n = 1, mean = mu_g1, sd = sd_g1)
      }
    }
    
    d     <- data.frame(id = 1:N_max, trt = rand_vec, y = y)
    #----- IA dataset 
    d_IA  <- subset(d, id <= IA_N)
    IA_g0 <- subset(d_IA, trt == 0)$y
    IA_g1 <- subset(d_IA, trt == 1)$y 
    #----- FA dataset (fixed)
    d_FA  <- subset(d, id <= FA_N)
    FA_g0 <- subset(d_FA, trt == 0)$y 
    FA_g1 <- subset(d_FA, trt == 1)$y 
    
    #----- Calc. CP 
    IA_Z1 <- (mean(IA_g1) - mean(IA_g0)) / sqrt(var(IA_g0)/length(IA_g0) + var(IA_g1)/length(IA_g1))
    CP    <- 1 - pnorm(q = (qnorm(p = 0.975)*sqrt(FA_N) - IA_Z1*sqrt(IA_N))/sqrt(FA_N-IA_N) - IA_Z1*sqrt(FA_N-IA_N)/sqrt(IA_N))
    
    pval <- t.test(x = FA_g0, y = FA_g1, alternative = "two.sided", conf.level = 0.95)$p.value
    d_sim$Sig_fixed[i] <- ifelse(pval < 0.05, 1, 0)
    
    if (CP >= 0.8) { 
      
      d_sim$Zone[i]    <- "Favorable"
      d_sim$CP[i]      <- CP 
      d_sim$SS[i]      <- FA_N 
      d_sim$Sig_ada[i] <- ifelse(pval < 0.05, 1, 0)
      
    } else if (0.365 <= CP && CP < 0.8) {
      
      N_max  <- sim_summary$N_max[k]
      N2_adj <- (FA_N + 1):N_max 
      CP_adj <- 1 - pnorm(q = (qnorm(p = 0.975)*sqrt(N2_adj) - IA_Z1*sqrt(IA_N))/sqrt(N2_adj-IA_N) - IA_Z1*sqrt(N2_adj-IA_N)/sqrt(IA_N))
      d_CP   <- data.frame(N2_adj = N2_adj, CP_adj = CP_adj)
      
      if (nrow(subset(d_CP, CP_adj >= 0.8)) == 0) {
        
        d_sim$Zone[i]   <- "Promising" 
        d_sim$CP[i]     <- 1 - pnorm(q = (qnorm(p = 0.975)*sqrt(N_max) - IA_Z1*sqrt(IA_N))/sqrt(N_max-IA_N) - IA_Z1*sqrt(N_max-IA_N)/sqrt(IA_N))
        d_sim$SS[i]     <- N_max 
        g0_ada <- subset(d, id <= N_max & trt == 0)$y
        g1_ada <- subset(d, id <= N_max & trt == 1)$y 
        pval <- t.test(x = g0_ada, y = g1_ada, alternative = "two.sided", conf.level = 0.95)$p.value
        d_sim$Sig_ada[i]<- ifelse(pval < 0.05, 1, 0)
        
      } else {
        d_sim$Zone[i] <- "Promising" 
        d_sim$CP[i]   <- min(subset(d_CP, CP_adj >= 0.8)$CP_adj)
        d_sim$SS[i]   <- min(subset(d_CP, CP_adj >= 0.8)$N2_adj)
        g0_ada <- subset(d, id <= d_sim$SS[i] & trt == 0)$y 
        g1_ada <- subset(d, id <= d_sim$SS[i] & trt == 1)$y 
        pval <- t.test(x = g0_ada, y = g1_ada, alternative = "two.sided", conf.level = 0.95)$p.value
        d_sim$Sig_ada[i]<- ifelse(pval < 0.05, 1, 0)
        
      } 
    } else if (CP < 0.365) {
      d_sim$Zone[i] <- "Unfavorable" 
      d_sim$CP[i]   <- CP 
      d_sim$SS[i]   <- FA_N 
      pval <- t.test(x = FA_g0, y = FA_g1, alternative = "two.sided", conf.level = 0.95)$p.value
      d_sim$Sig_ada[i] <- ifelse(pval < 0.05, 1, 0)
      
    }
      
  }
  
  sim_summary$Pr_fav[k]   <- nrow(subset(d_sim, Zone == "Favorable"))/nrow(d_sim)
  sim_summary$Pr_prom[k]  <- nrow(subset(d_sim, Zone == "Promising"))/nrow(d_sim)
  sim_summary$Pr_unfav[k] <- nrow(subset(d_sim, Zone == "Unfavorable"))/nrow(d_sim) 
  
  sim_summary$Sig_fixed[k]      <- mean(d_sim$Sig_fixed)
  sim_summary$Sig_fixed_fav[k]  <- mean(subset(d_sim, Zone == "Favorable")$Sig_fixed)
  sim_summary$Sig_fixed_prom[k] <- mean(subset(d_sim, Zone == "Promising")$Sig_fixed)
  sim_summary$Sig_fixed_unfav[k]<- mean(subset(d_sim, Zone == "Unfavorable")$Sig_fixed)
  
  sim_summary$Sig_ada[k]      <- mean(d_sim$Sig_ada)
  sim_summary$Sig_ada_fav[k]  <- mean(subset(d_sim, Zone == "Favorable")$Sig_ada)
  sim_summary$Sig_ada_prom[k] <- mean(subset(d_sim, Zone == "Promising")$Sig_ada)
  sim_summary$Sig_ada_unfav[k]<- mean(subset(d_sim, Zone == "Unfavorable")$Sig_ada)
  
  sim_summary$ASN[k]      <- mean(d_sim$SS)
  sim_summary$ASN_fav[k]  <- mean(subset(d_sim, Zone == "Favorable")$SS)
  sim_summary$ASN_prom[k] <- mean(subset(d_sim, Zone == "Promising")$SS)
  sim_summary$ASN_unfav[k]<- mean(subset(d_sim, Zone == "Unfavorable")$SS) 
  
}



#----- read simulation result 
d <- read.csv(file = paste0(base_dir, "/data/a5-Plan4.csv"), sep = ",", header = TRUE)
d_h0 <- read.csv(file = paste0(base_dir, "/data/a5-Plan4-H0.csv"), sep = ",", header = TRUE)




#----- Probability of interim outcome
library(dplyr)
d_ana <- d %>% 
  select(delta, sigma, Pr_fav, Pr_prom, Pr_unfav) %>% 
  mutate(sigma = paste0("SD = ", format(sigma, nsmall = 1))) %>% 
  ungroup() 

library(reshape2)
d_ana <- melt(d_ana, id.vars = c("delta", "sigma"), measure.vars = c("Pr_fav", "Pr_prom", "Pr_unfav"), 
              variable.name = "zone", value.name = "aval")

d_ana <- d_ana %>% 
  mutate(zone = case_when(zone == "Pr_fav" ~  "Favorable", 
                          zone == "Pr_prom" ~ "Promising", 
                          zone == "Pr_unfav" ~ "Unfavorable", 
                          TRUE ~ "")) %>%
  ungroup() 


library(ggplot2)
p <- ggplot(data = d_ana, aes(x = delta, y = aval, colour = zone))
p <- p + facet_wrap( ~ sigma, ncol = 3)
p <- p + geom_line(size = 1.5)
p <- p + scale_x_continuous(limits = c(1.0, 2.5), breaks = seq(from = 1.0, to = 2.5, by = 0.5))
p <- p + scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(from = 0.0, to = 1.0, by = 0.25))
p <- p + xlab(expression(delta))
p <- p + ylab("Probability")
p <- p + theme(strip.text = element_text(size = 14, colour = "black"))
p <- p + theme(axis.title = element_text(size = 14, colour = "black"))
p <- p + theme(axis.text = element_text(size = 14, colour = "black"))
p <- p + theme(legend.title = element_text(size = 14, colour = "black"))
p <- p + theme(legend.text = element_text(size = 14, colour = "black"))
ggsave(file = paste0(base_dir, "/output/a5-1-outcomeProb.png"), plot = p, dpi = 400, w = 8, h = 4)



#----- Overall unconditional power
library(dplyr)
d_ana <- d %>% 
  select(delta, sigma, Sig_fixed, Sig_ada) %>% 
  mutate(sigma = paste0("SD = ", format(sigma, nsmall = 1))) %>% 
  ungroup() 

library(reshape2)
d_ana <- melt(d_ana, id.vars = c("delta", "sigma"), measure.vars = c("Sig_fixed", "Sig_ada"), 
              variable.name = "design", value.name = "aval")

d_ana <- d_ana %>% 
  mutate(design = case_when(design == "Sig_fixed" ~  "Fixed", 
                            design == "Sig_ada" ~ "Adaptive", 
                            TRUE ~ "")) %>%
  ungroup() 


library(ggplot2)
p <- ggplot(data = d_ana, aes(x = delta, y = aval, colour = design))
p <- p + facet_wrap( ~ sigma, ncol = 3)
p <- p + geom_line(size = 1.5)
p <- p + scale_x_continuous(limits = c(1.0, 2.5), breaks = seq(from = 1.0, to = 2.5, by = 0.5))
p <- p + scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(from = 0.0, to = 1.0, by = 0.25))
p <- p + xlab(expression(delta))
p <- p + ylab("Probability")
p <- p + theme(strip.text = element_text(size = 14, colour = "black"))
p <- p + theme(axis.title = element_text(size = 14, colour = "black"))
p <- p + theme(axis.text = element_text(size = 14, colour = "black"))
p <- p + theme(legend.title = element_text(size = 14, colour = "black"))
p <- p + theme(legend.text = element_text(size = 14, colour = "black"))
ggsave(file = paste0(base_dir, "/output/a5-2-unconditionalPower.png"), plot = p, dpi = 400, w = 8, h = 4)




#----- Average sample size
library(dplyr)
d_ana <- d %>% 
  select(delta, sigma, ASN) %>% 
  mutate(sigma = paste0("SD = ", format(sigma, nsmall = 1))) %>% 
  ungroup() 


library(ggplot2)
p <- ggplot(data = d_ana, aes(x = delta, y = ASN, colour = sigma))
# p <- p + facet_wrap( ~ sigma, ncol = 3)
p <- p + geom_line(size = 1.5)
p <- p + geom_hline(yintercept = 442, linetype = "dotted", colour = "black", size = 1.5)
p <- p + scale_x_continuous(limits = c(1.0, 2.5), breaks = seq(from = 1.0, to = 2.5, by = 0.5))
p <- p + scale_y_continuous(limits = c(400, 520), breaks = seq(from = 400, to = 500, by = 25))
p <- p + xlab(expression(delta))
p <- p + ylab("Average sample size")
p <- p + theme(strip.text = element_text(size = 14, colour = "black"))
p <- p + theme(axis.title = element_text(size = 14, colour = "black"))
p <- p + theme(axis.text = element_text(size = 14, colour = "black"))
p <- p + theme(legend.title = element_text(size = 14, colour = "black"))
p <- p + theme(legend.text = element_text(size = 14, colour = "black"))
ggsave(file = paste0(base_dir, "/output/a5-3-ASN.png"), plot = p, dpi = 400, w = 8, h = 4)



#----- type I error rate
library(reshape2)
d_ana <- melt(d_h0, id.vars = c("delta", "sigma"), measure.vars = c("Sig_fixed", "Sig_ada"), 
              variable.name = "design", value.name = "aval")

library(dplyr)
d_ana <- d_ana %>% 
  mutate(design = case_when(design == "Sig_fixed" ~ "Fixed", 
                            design == "Sig_ada" ~ "Adaptive", 
                            TRUE ~ "")) %>% 
 ungroup() 

library(ggplot2)
p <- ggplot(data = d_ana, aes(x = sigma, y = aval, colour = design))
p <- p + geom_line(size = 1.5)
p <- p + geom_hline(yintercept = 0.05, linetype = "dotted", colour = "red", size = 1.5)
p <- p + scale_x_continuous(limits = c(5.0, 10.0), breaks = seq(from = 5.0, to = 10.0, by = 1))
p <- p + scale_y_continuous(limits = c(0.03, 0.06), breaks = seq(from = 0.03, to = 0.06, by = 0.01))
p <- p + xlab(expression(sigma))
p <- p + ylab("Type I error rate")
p <- p + theme(axis.title = element_text(size = 14, colour = "black"))
p <- p + theme(axis.text = element_text(size = 14, colour = "black"))
p <- p + theme(legend.title = element_text(size = 14, colour = "black"))
p <- p + theme(legend.text = element_text(size = 14, colour = "black"))
ggsave(file = paste0(base_dir, "/output/a5-4-TypeIError.png"), plot = p, dpi = 400, w = 8, h = 4)







