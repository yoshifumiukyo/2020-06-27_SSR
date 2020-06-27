############################################################
# R-project                                                #
# Program      : a1-ConventionalCalc.r                     #
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


alpha <- 0.05 
beta  <- 0.2  # power = 0.8
delta <- c(1.6, 1.7, 1.8, 1.9, 2.0)  
sigma <- 7.5 

N_group <- 2 * (qnorm(p = 1 - alpha/2) + qnorm(p = 1 - beta))^2 * sigma^2 / delta^2
N_total <- 2 * N_group 

#----- Plan 1
N_group <- 442/2
Zbeta <- (delta/sigma) * sqrt(N_group/2) - qnorm(p = 1 - alpha/2)
pnorm(q = Zbeta)


#----- Plan 2
N_group <- 690/2
Zbeta <- (delta/sigma) * sqrt(N_group/2) - qnorm(p = 1 - alpha/2)
pnorm(q = Zbeta)




