#Eleanor Sanderson
#25th July 2023
#Code to replicate SNP effects in Tian and Burgess IJE 2023


library(dplyr)
library(tidyverse)
library(ggplot2)

rm(list = ls(all.names = TRUE))
set.seed(11)

reps = 1000
results <- data.frame()


for(r in 1:reps){

A1 <- rnorm(30,0,0.05)
A2 <- rnorm(30,0,0.15)
A3 <- rnorm(30,0.10,0.01)
A4 <- rnorm(30,0,1)

t <- 1:50
alpha = matrix( , nrow = 50, ncol = 30)
for (i in 1:50){
  for(j in 1:30){
alpha[i,j] <- A1[j] +  A2[j]*cos(A3[j]*i-A4[j])
  }
}

alpha <- data.frame(alpha)


alpha$sum <- alpha[,1]

for(i in 2:30){
  alpha$sum <- alpha$sum + alpha[,i]
}

for(t in 1:50){
  results[r,t] <- alpha[t,"sum"]
}


a <- cor(t(alpha[,1:30]))

#15, 30 with 10 and 50

results[r,"10_15"] <- a[10,15]
results[r,"10_50"] <- a[10,50]
results[r,"30_50"] <- a[30,50]
results[r,"15_50"] <- a[15,50]
results[r,"10_30"] <- a[10,30]


}
##plot mean across the simulation


plotdat <- data.frame(colMeans(results[,1:50]))
plotdat <- plotdat %>% 
  rowid_to_column(var='time') %>%
  rename("snpeffect" = "colMeans.results...1.50..")


ggplot(data = plotdat, aes(x=time)) + 
  geom_line(aes(y=snpeffect)) +
  geom_hline(aes(yintercept = 0), col='black', linetype=2) +
  theme_bw() +
  labs(x="Time", y = "Total SNP-exposure effect across 30 SNPs", title = "(a) Mean effects across 1000 repetitions")

ggsave("mean_plot_changeparams.pdf", width = 7, height = 5)


##plot for one repetition

randno <- round(runif(1,1,reps))

plotdat <- data.frame(t(results[randno,1:50]))
plotdat <- plotdat %>% 
  rowid_to_column(var='time') 
  
colnames(plotdat)[2] = "snpeffect"


ggplot(data = plotdat, aes(x=time)) + 
  geom_line(aes(y=snpeffect)) +
  geom_hline(aes(yintercept = 0), col='black', linetype=2) +
  theme_bw() +
  labs(x="Time", y = "Total SNP-exposure effect across 30 SNPs", title = "(b) Effects for one randomly selected simulation")

ggsave("singlerep_plot_changeparams.pdf", width = 7, height = 5)



###looking at the correlations between the time periods

correlations <- results[,51:ncol(results)]


#estimated effect at time 15 scenario D for each repetition
correlations$D_15 <- correlations[,"10_15"]*0.4
correlations$D_50 <- -0.8 + correlations[,"10_50"]*0.4



