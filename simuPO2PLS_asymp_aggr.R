## ORGANIZATION
## Plots for increasing sample size, increasing dimensionality and increasing power are made
## How to run this script:
## If all went well in simuPO2PLS.R, you've generated the SIMUasymp***.RData files
## Source this file to generate the plots in one SIMUasympCombined.pdf 

library(parallel)
library(tidyverse)
library(magrittr)
library(OmicsPLS)
library(PO2PLS)
library(Cairo)

############ Sample size ########
Dat_asympN <- lapply(list.files(pattern="SIMUasympH0_N"), function(ee) {
  load(ee)
  ret = sapply(outp.H0, function(e) (e$pvals.B$asymp)) %>% unlist
  cbind(SampleSize = outp.H0[[1]]$samp, Dim_p = outp.H0[[1]]$dims, asymp = ret<0.05+0) 
}) %>% Reduce(f=rbind) %>% as_tibble %>% arrange(SampleSize)
#colnames(Dat_asymp) <- c("N_1000", "N_50")
cat("\nResults for asymp only \n")
cat(paste0("Replicates: ",length(outp.H0)," \n"))
Dat_asympN %>% group_by(SampleSize) %>% summarise(Asymptotic = mean(asymp)) %>% print
K_repl <- Dat_asympN %>% group_by(SampleSize) %>% summarise(l = length(asymp)) %>% pull(l) %>% `[[`(1)
Dat_asympN$SampleSize %<>% paste0("N=",.) %>% factor(., levels = unique(.))

# try(dev.new(), silent=TRUE)
# CairoPDF(file = "SIMUasympSampSize.pdf", height=7, width=10)
# Dat_asympN %>% group_by(SampleSize) %>% summarise(TypeI = mean(asymp)) %>% #`[`(c(2:4,1),) %>%
#   mutate(Up = TypeI + 2*sqrt(TypeI*(1-TypeI)/K_repl), Lo = TypeI - 2*sqrt(TypeI*(1-TypeI)/K_repl)) %>%
#   ggplot(aes(x=SampleSize, y=TypeI, col="red", group=1)) + geom_line() + geom_point() +
#   geom_point(aes(x=SampleSize, y=Up, col="red", group=1)) + geom_point(aes(x=SampleSize, y=Lo,col="red", group=1)) +
#   geom_hline(yintercept=0.05, lty=2) + ylim(c(0,0.2)) +
#   theme_bw() + labs(x="Sample Size", y="Type I error") + theme(text=element_text(size=16), legend.position="none")
# dev.off()

############ Dimensionality ########
Dat_asympD <- lapply(list.files(pattern="SIMUasympH0_p"), function(ee) {
  load(ee)
  ret = sapply(outp.H0, function(e) (e$pvals.B$asymp)) %>% unlist
  cbind(SampleSize = outp.H0[[1]]$samp, Dim_p = outp.H0[[1]]$dims, asymp = ret<0.05+0) 
}) %>% Reduce(f=rbind) %>% as_tibble
#colnames(Dat_asymp) <- c("N_1000", "N_50")
cat("\nResults for asymp only \n")
cat(paste0("Replicates: ",length(outp.H0)," \n"))
Dat_asympD %>% group_by(Dim_p) %>% summarise(Asymptotic = mean(asymp)) %>% print
K_repl <- Dat_asympD %>% group_by(Dim_p) %>% summarise(l = length(asymp)) %>% pull(l) %>% `[[`(1)
Dat_asympD$Dim_p %<>% paste0("p=",.) 

# try(dev.new(), silent=TRUE)
# CairoPDF(file = "SIMUasympDim.pdf", height=7, width=10)
# Dat_asympD %>% group_by(Dim_p) %>% summarise(TypeI = mean(asymp)) %>% #`[`(c(2:4,1),) %>%
#   mutate(Up = TypeI + 2*sqrt(TypeI*(1-TypeI)/K_repl), Lo = TypeI - 2*sqrt(TypeI*(1-TypeI)/K_repl)) %>%
#   ggplot(aes(x=Dim_p, y=TypeI, col="red", group=1)) + geom_line() + geom_point() +
#   geom_point(aes(x=Dim_p, y=Up, col="red", group=1)) + geom_point(aes(x=Dim_p, y=Lo,col="red", group=1)) +
#   geom_hline(yintercept=0.05, lty=2) + ylim(c(0,0.2)) +
#   theme_bw() + labs(x="Dimensionality", y="Type I error") + theme(text=element_text(size=16), legend.position="none")
# dev.off()



############ Power ########
Dat <- lapply(list.files(pattern = "SIMUasympAll_N"), function(ee) {
  load(ee) 
  ret0 <- sapply(outp.H0, function(e) (e$pvals.B)) %>% unlist %>% matrix(nrow=4) %>% t
  ret <- cbind(SampleSize = outp.H0[[1]]$samp, 
    Dim_p = outp.H0[[1]]$dims,
    pwr = str_split(ee,"_")[[1]][4] %>% str_sub(start=4) %>% as.numeric,
    adj_pValue = (ret0 < 0.05)+0)
    colnames(ret)[-(1:3)] <- names(outp.H0[[1]]$pvals.B)
  return(ret)
}) %>% Reduce(f=rbind) %>% as_tibble
cat("\nResults for all scenarios \n")
cat(paste0("Replicates: ",length(outp.H0)," \n"))
Dat %>% group_by(SampleSize, Dim_p, pwr) %>% summarise(Asymp = mean(asymp), Boot = mean(boot), nonBoot = mean(nonboot), perm = mean(perm)) %>% print
K_repl <- Dat %>% group_by(SampleSize, Dim_p, pwr) %>% summarise(l = length(asymp)) %>% pull(l) %>% `[[`(1)
Dat$pwr %<>% as.factor
Dat$SampleSize <- ifelse(Dat$SampleSize == 50, "N=50, p=200", "N=500, p=20")

# try(dev.new(), silent=TRUE)
# CairoPDF(file = "SIMUasympAll.pdf", height=7, width=10)
# Dat %>% group_by(SampleSize, Dim_p, pwr) %>% summarise(Asymptotic = mean(asymp), `Param. Boot.` = mean(boot), `Non-param. Boot.` = mean(nonboot), Permutation = mean(perm)) %>% 
#   gather(key="Approach", value="Power", -(1:3)) %>% mutate(Up = Power + 2*sqrt(Power*(1-Power)/K_repl), Lo = Power - 2*sqrt(Power*(1-Power)/K_repl)) %>%
#   ggplot(aes(x=pwr, y=Power, col=Approach, shape=Approach, group=Approach)) + geom_line() + #geom_point() + 
#   geom_point(aes(x=pwr, y=Up, col=Approach, shape=Approach, group=Approach)) + geom_point(aes(x=pwr, y=Lo, col=Approach, shape=Approach, group=Approach)) +
#   facet_grid(SampleSize~.) + 
#   geom_hline(yintercept=0.05, lty=2) + 
#   theme_bw() + labs(x="Effect Size") + theme(text=element_text(size=18))
# dev.off()


K_repl <- Dat_asympN %>% group_by(SampleSize) %>% summarise(l = length(asymp)) %>% pull(l) %>% `[[`(1)
Dat_asympN %>% group_by(SampleSize) %>% summarise(TypeI = mean(asymp)) %>% #`[`(c(2:4,1),) %>%
  mutate(Up = TypeI + 2*sqrt(TypeI*(1-TypeI)/K_repl), Lo = TypeI - 2*sqrt(TypeI*(1-TypeI)/K_repl)) %>%
  ggplot(aes(x=SampleSize, y=TypeI, col="red", group=1)) + geom_line() + geom_point() +
  geom_point(aes(x=SampleSize, y=Up, col="red", group=1)) + geom_point(aes(x=SampleSize, y=Lo,col="red", group=1)) +
  geom_hline(yintercept=0.05, lty=2) + ylim(c(0,0.2)) +
  theme_bw() + labs(x="Sample Size", y="Type I error") + theme(text=element_text(size=14), legend.position="none") -> p1

K_repl <- Dat_asympD %>% group_by(Dim_p) %>% summarise(l = length(asymp)) %>% pull(l) %>% `[[`(1)
Dat_asympD %>% group_by(Dim_p) %>% summarise(TypeI = mean(asymp)) %>% #`[`(c(2:4,1),) %>%
  mutate(Up = TypeI + 2*sqrt(TypeI*(1-TypeI)/K_repl), Lo = TypeI - 2*sqrt(TypeI*(1-TypeI)/K_repl)) %>%
  ggplot(aes(x=Dim_p, y=TypeI, col="red", group=1)) + geom_line() + geom_point() +
  geom_point(aes(x=Dim_p, y=Up, col="red", group=1)) + geom_point(aes(x=Dim_p, y=Lo,col="red", group=1)) +
  geom_hline(yintercept=0.05, lty=2) + ylim(c(0,0.2)) +
  theme_bw() + labs(x="Dimensionality", y="Type I error") + theme(text=element_text(size=14), legend.position="none") -> p2

K_repl <- Dat %>% group_by(SampleSize, Dim_p, pwr) %>% summarise(l = length(asymp)) %>% pull(l) %>% `[[`(1)
Dat %>% group_by(SampleSize, Dim_p, pwr) %>% summarise(Asymptotic = mean(asymp), `Param. Boot.` = mean(boot), `Non-param. Boot.` = mean(nonboot), Permutation = mean(perm)) %>% 
  gather(key="Approach", value="Power", -(1:3)) %>% mutate(Up = Power + 2*sqrt(Power*(1-Power)/K_repl), Lo = Power - 2*sqrt(Power*(1-Power)/K_repl)) %>%
  ggplot(aes(x=pwr, y=Power, col=Approach, shape=Approach, group=Approach)) + geom_line() + #geom_point() + 
  geom_point(aes(x=pwr, y=Up, col=Approach, shape=Approach, group=Approach)) + geom_point(aes(x=pwr, y=Lo, col=Approach, shape=Approach, group=Approach)) +
  facet_grid(SampleSize~.) + 
  geom_hline(yintercept=0.05, lty=2) + 
  theme_bw() + labs(x="Effect Size") + theme(text=element_text(size=14), legend.position="top") -> p3

try(dev.new(), silent=TRUE)
CairoPDF(file = "SIMUasympCombined.pdf", height=7, width=14)
plot.layout <- rbind(c(1, 1, 3, 3, 3), c(2, 2, 3, 3, 3))
gridExtra::grid.arrange(grobs = list(p1, p2, p3), layout_matrix = plot.layout)
dev.off()


