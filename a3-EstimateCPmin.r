############################################################
# R-project                                                #
# Program      : a3-EstimateCPmin.r                        #
# Protocol     :                                           #
# Date         : <2020/06/02 (Tue) 23:07>                  #
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




N1 <- 208
N2 <- 442 
Nmax <- 884 
N2_tild <- (N2 - N1):(Nmax - N1)


delta <- 2.0 
delta_hat <- seq(from = 1.0, to = 2.0, by = 0.01) 
sigma <- 7.5 


d_CP_min <- data.frame(id = 1:length(delta_hat), CP_n2 = 0, b_value = 0)

for (i in 1:length(delta_hat)){ 
  
  Z1 <- delta_hat[i] / sqrt(4 * sigma^2 / N1)
  CP_n2 <- 1 - pnorm(q = (qnorm(p = 0.975) * sqrt(N2) - Z1 * sqrt(N1))/sqrt(N2-N1) - Z1 * sqrt(N2-N1) / sqrt(N1))
  
  n2_dash <- (N1/(Z1^2)) * ((qnorm(p = 0.975)*sqrt(N2)-Z1*sqrt(N1))/sqrt(N2-N1) + qnorm(p = 0.8))^2
  N2_ast <- min(n2_dash, Nmax-N1)
  
  b_value <- (1/sqrt(N1 + N2_ast)) * (sqrt(N2_ast/(N2-N1))*(qnorm(p = 0.975)*sqrt(N2) - Z1*sqrt(N1))+Z1*sqrt(N1))
  
  d_CP_min$CP_n2[i]   <- CP_n2 
  d_CP_min$b_value[i] <- b_value  
  
}


library(ggplot2)
p <- ggplot(data = d_CP_min, aes(x = CP_n2, y = b_value))
p <- p + theme_bw() 
p <- p + geom_line(size = 1.5)
p <- p + geom_hline(yintercept = 1.96, colour = "red", linetype = "dotted", size = 1.5)
p <- p + geom_vline(xintercept = 0.365, colour = "blue", linetype = "dotted", size = 1.5)
p <- p + xlab("CP")
p <- p + ylab("b")
p <- p + scale_x_continuous(limits = c(0.2, 0.9), breaks = seq(from = 0.2, to = 0.9, by = 0.1))
p <- p + scale_y_continuous(limits = c(1.90, 2.05), breaks = seq(from = 1.90, to = 2.05, by = 0.05))
p <- p + theme(axis.title = element_text(size = 14, colour = "black"))
p <- p + theme(axis.text = element_text(size = 14, colour = "black"))
ggsave(file = paste0(base_dir, "/output/a3-CPmin.png"), plot = p, dpi = 400, w = 6, h = 4)




