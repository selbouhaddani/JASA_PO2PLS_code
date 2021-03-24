## Organization
## First plots for the normal study are made, then rank, then distribution, then case-control outcome
## How to run this script:
## If all went well in simuPO2PLS.R, you've generated the SIMU***.RData files
## Source this file to generate the plots, in SIMU***.pdf format


library(parallel)
library(tidyverse)
library(magrittr)
library(OmicsPLS)
library(PO2PLS)
library(pryr)
library(gridExtra)
library(Cairo)

######## Normal plots ########
print("Normal plots")
#### Read in simulation results
out.dat <- sapply(simplify = FALSE, 
  list.files(pattern = paste0("SIMU_N")), # lapply over each RData file with output
  function(ee) { # prepare dataframe with output
    load(ee)
    outpt <- sapply(simplify=F,1:length(outp), function(ii) outp[[ii]])
    dat <- list()
    dat$Pred <- outpt %>% sapply(function(ee) unlist(ee$Pred)) %>% t
    dat$trainPred <- outpt %>% sapply(function(ee) unlist(ee$trainPred)) %>% t
    dat$TPR <- outpt %>% sapply(function(ee) colMeans(ee$TPR)) %>% t
    scenario = ee %>% str_split("_") %>% `[[`(1)
    dat$SampSz = ifelse(scenario[2] == "Nlow", "N=100", "N=1000")
    dat$Dims = ifelse(scenario[3] == "plow", "p=2e3", "p=2e4")
    dat$Noise = ifelse(scenario[4] == "nlow", "Low", "High")
    dat$Heter = ifelse(scenario[5] == "heterogen.RData", "Heterogen", "Homogen")
    return(dat)
})

#### Prepare plotting data
plotdat <- list()
plotdat$TPR <- out.dat %>% 
  lapply(function(ee) data.frame(ee$TPR) %>% 
           bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                     Heter = ee$Heter)) %>% Reduce(f=bind_rows)
plotdat$Pred <- out.dat %>% 
  lapply(function(ee) data.frame(ee$Pred) %>% 
           bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                     Heter = ee$Heter)) %>% Reduce(f=bind_rows)
plotdat$trainPred <- out.dat %>% 
  lapply(function(ee) data.frame(ee$trainPred) %>% 
           bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                     Heter = ee$Heter)) %>% Reduce(f=bind_rows)

fctname <- c("PO2PLS", "O2PLS", "PLS", "PPLS", "SIFA")
lay <- rbind(c(1,1,1,3,3,3,3),c(2,2,2,3,3,3,3))

try(dev.off(),silent = T)
CairoPDF(file = paste0("SIMU_plot"), height=10, width=18)
p1a <- plotdat$TPR %>% gather(key="Method", value="TPR", all_of(fctname)) %>% 
  mutate(Method=factor(Method,labels=fctname,levels=fctname), Noise=Noise%>%factor%>%fct_rev) %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method,y=TPR)) + geom_boxplot(aes(col=Noise)) + theme_bw() + 
  facet_grid(. ~ SampSz) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), 
        legend.title=element_text(face="bold"))
p1b <- plotdat$TPR %>% mutate(PPLS=PPLS-PO2PLS, PLS=PLS-PO2PLS, O2PLS=O2PLS-PO2PLS,SIFA=SIFA-PO2PLS) %>%
  gather(key="Method", value="TPR", all_of(fctname)) %>% 
  mutate(Method=factor(Method,labels=fctname,levels=fctname), Noise=Noise%>%factor%>%fct_rev) %>% 
  filter(Method!="PO2PLS") %>%
  filter(Method != "SIFA") %>% 
  #filter(Method != "PPLS") %>%
  ggplot(aes(x=Method,y=TPR)) + geom_boxplot(aes(col=Noise)) + theme_bw() + ylab("Diff in TPR") + geom_hline(yintercept=0) + 
  facet_grid(. ~ SampSz) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
p2 <- plotdat$Pred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Test") %>% 
  bind_rows(plotdat$trainPred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Train")) %>% 
  mutate(Type = factor(Type, labels = c("train", "test"), 
                       levels = c("Train", "Test")),
         Noise=Noise%>%factor%>%fct_rev) %>%
  gather(key="Method",value="RMSEP",1:6) %>% mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
  filter(Method != "TRuE") %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method, y=RMSEP)) + geom_boxplot(aes(col=Type)) + theme_bw() + 
  facet_grid(Noise*Heter ~ SampSz) + 
  geom_hline(yintercept=median(plotdat$Pred$TRuE %>% sqrt), lwd=1) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
grid.arrange(grobs=list(p1a,p1b,p2), layout_matrix=lay)
dev.off()


######### Rank misspecification Plots ####
print("Rank misspecification plots")
#### Read in simulation results
out.dat <- sapply(simplify = FALSE, 
                  list.files(pattern = paste0("SIMUrank_N")), # lapply over each RData file with output
                  function(ee) { # prepare dataframe with output
                    load(ee)
                    outpt <- sapply(simplify=F,1:length(outp), function(ii) outp[[ii]])
                    dat <- list()
                    dat$Pred <- outpt %>% sapply(function(ee) unlist(ee$Pred)) %>% t
                    dat$trainPred <- outpt %>% sapply(function(ee) unlist(ee$trainPred)) %>% t
                    dat$TPR <- outpt %>% sapply(function(ee) colMeans(ee$TPR)) %>% t
                    scenario = ee %>% str_split("_") %>% `[[`(1)
                    dat$SampSz = ifelse(scenario[2] == "Nlow", "N=100", "N=1000")
                    dat$Dims = ifelse(scenario[3] == "plow", "p=2e3", "p=2e4")
                    dat$Noise = ifelse(scenario[4] == "nlow", "Low", "High")
                    dat$Heter = ifelse(scenario[5] == "heterogen.RData", "Heterogen", "Homogen")
                    return(dat)
                  })

#### Prepare plotting data
plotdat <- list()
plotdat$TPR <- out.dat %>% 
  lapply(function(ee) data.frame(ee$TPR) %>% 
           bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                     Heter = ee$Heter)) %>% Reduce(f=bind_rows)
plotdat$Pred <- out.dat %>% 
  lapply(function(ee) data.frame(ee$Pred) %>% 
           bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                     Heter = ee$Heter)) %>% Reduce(f=bind_rows)
plotdat$trainPred <- out.dat %>% 
  lapply(function(ee) data.frame(ee$trainPred) %>% 
           bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                     Heter = ee$Heter)) %>% Reduce(f=bind_rows)

fctname <- c("PO2PLS", "O2PLS", "PLS", "PPLS", "SIFA")
lay <- rbind(c(1,1,1,3,3,3,3),c(2,2,2,3,3,3,3))

try(dev.off(),silent = T)
CairoPDF(file = paste0("SIMUrank_plot"), height=10, width=18)
p1a <- plotdat$TPR %>% gather(key="Method", value="TPR", all_of(fctname)) %>% 
  mutate(Method=factor(Method,labels=fctname,levels=fctname), Noise=Noise%>%factor%>%fct_rev) %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method,y=TPR)) + geom_boxplot(aes(col=Noise)) + theme_bw() + 
  facet_grid(. ~ SampSz) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), 
        legend.title=element_text(face="bold"))
p1b <- plotdat$TPR %>% mutate(PPLS=PPLS-PO2PLS, PLS=PLS-PO2PLS, O2PLS=O2PLS-PO2PLS,SIFA=SIFA-PO2PLS) %>%
  gather(key="Method", value="TPR", all_of(fctname)) %>% 
  mutate(Method=factor(Method,labels=fctname,levels=fctname), Noise=Noise%>%factor%>%fct_rev) %>% 
  filter(Method!="PO2PLS") %>%
  filter(Method != "SIFA") %>% 
  #filter(Method != "PPLS") %>%
  ggplot(aes(x=Method,y=TPR)) + geom_boxplot(aes(col=Noise)) + theme_bw() + ylab("Diff in TPR") + geom_hline(yintercept=0) + 
  facet_grid(. ~ SampSz) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
p2 <- plotdat$Pred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Test") %>% 
  bind_rows(plotdat$trainPred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Train")) %>% 
  mutate(Type = factor(Type, labels = c("train", "test"), 
                       levels = c("Train", "Test")),
         Noise=Noise%>%factor%>%fct_rev) %>%
  gather(key="Method",value="RMSEP",1:6) %>% mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
  filter(Method != "TRuE") %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method, y=RMSEP)) + geom_boxplot(aes(col=Type)) + theme_bw() + 
  facet_grid(Noise*Heter ~ SampSz) + 
  geom_hline(yintercept=median(plotdat$Pred$TRuE %>% sqrt), lwd=1) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
grid.arrange(grobs=list(p1a,p1b,p2), layout_matrix=lay)
dev.off()


###### Distribution misspecification ########
print("Distribution plots")
distribs <- c("student_t", "poiss", "binom")
for(distr_i in 1:3){
  try(dev.off(),silent = T)
  CairoPDF(file = paste0("SIMUdistr",distribs[distr_i],"_plot"), height=10, width=18)
  
  out.dat <- sapply(simplify = FALSE, 
                    list.files(pattern = paste0("SIMUdistr_N")), # lapply over each RData file with output
                    function(ee) { # prepare dataframe with output
                      load(ee)
                      outpt <- sapply(simplify=F,1:length(outp), function(ii) outp[[ii]][[distr_i]])
                      dat <- list()
                      dat$Pred <- outpt %>% sapply(function(ee) unlist(ee$Pred)) %>% t
                      dat$trainPred <- outpt %>% sapply(function(ee) unlist(ee$trainPred)) %>% t
                      dat$TPR <- outpt %>% sapply(function(ee) colMeans(ee$TPR)) %>% t
                      scenario = ee %>% str_split("_") %>% `[[`(1)
                      dat$SampSz = ifelse(scenario[2] == "Nlow", "N=100", "N=1000")
                      dat$Dims = ifelse(scenario[3] == "plow", "p=2e3", "p=2e4")
                      dat$Noise = ifelse(scenario[4] == "nlow", "Low", "High")
                      dat$Heter = ifelse(scenario[5] == "heterogen.RData", "Heterogen", "Homogen")
                      return(dat)
                    })
  
  ## Prepare plotting data
  plotdat <- list()
  plotdat$TPR <- out.dat %>% 
    lapply(function(ee) data.frame(ee$TPR) %>% 
             bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                       Heter = ee$Heter)) %>% Reduce(f=bind_rows)
  plotdat$Pred <- out.dat %>% 
    lapply(function(ee) data.frame(ee$Pred) %>% 
             bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                       Heter = ee$Heter)) %>% Reduce(f=bind_rows)
  plotdat$trainPred <- out.dat %>% 
    lapply(function(ee) data.frame(ee$trainPred) %>% 
             bind_cols(SampSz = ee$SampSz, Dims = ee$Dims, Noise = ee$Noise, 
                       Heter = ee$Heter)) %>% Reduce(f=bind_rows)
  
  fctname <- c("PO2PLS", "O2PLS", "PLS", "PPLS", "SIFA")
  lay <- rbind(c(1,1,1,3,3,3,3),c(2,2,2,3,3,3,3))
  
  p1a <- plotdat$TPR %>% gather(key="Method", value="TPR", all_of(fctname)) %>% 
    mutate(Method=factor(Method,labels=fctname,levels=fctname), Noise=Noise%>%factor%>%fct_rev) %>% 
    filter(Method != "SIFA") %>%
    ggplot(aes(x=Method,y=TPR)) + geom_boxplot(aes(col=Noise)) + theme_bw() + 
    facet_grid(. ~ SampSz) + 
    theme(text=element_text(size=14), axis.title=element_text(face="bold"), 
          legend.title=element_text(face="bold"))
  p1b <- plotdat$TPR %>% mutate(PPLS=PPLS-PO2PLS, PLS=PLS-PO2PLS, O2PLS=O2PLS-PO2PLS,SIFA=SIFA-PO2PLS) %>%
    gather(key="Method", value="TPR", all_of(fctname)) %>% 
    mutate(Method=factor(Method,labels=fctname,levels=fctname), Noise=Noise%>%factor%>%fct_rev) %>% 
    filter(Method!="PO2PLS") %>%
    filter(Method != "SIFA") %>% 
    #filter(Method != "PPLS") %>%
    ggplot(aes(x=Method,y=TPR)) + geom_boxplot(aes(col=Noise)) + theme_bw() + ylab("Diff in TPR") + geom_hline(yintercept=0) + 
    facet_grid(. ~ SampSz) + 
    theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
  p2 <- plotdat$Pred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Test") %>% 
    bind_rows(plotdat$trainPred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Train")) %>% 
    mutate(Type = factor(Type, labels = c("train", "test"), 
                         levels = c("Train", "Test")),
           Noise=Noise%>%factor%>%fct_rev) %>%
    gather(key="Method",value="RMSEP",1:6) %>% mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
    filter(Method != "TRuE") %>% 
    filter(Method != "SIFA") %>%
    ggplot(aes(x=Method, y=RMSEP)) + geom_boxplot(aes(col=Type)) + theme_bw() + 
    facet_grid(Noise*Heter ~ SampSz) + 
    geom_hline(yintercept=median(plotdat$Pred$TRuE %>% sqrt), lwd=1) + 
    theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
  grid.arrange(grobs=list(p1a,p1b,p2), layout_matrix=lay)
  dev.off()
}



######## Outcome plots ######
print("Outcome plots")
out.dat <- sapply(simplify = FALSE, 
                  list.files(pattern = paste0("SIMUoutcome.RData")), # lapply over each RData file with output
                  function(ee) { # prepare dataframe with output
                    load(ee)
                    outpt <- sapply(simplify=F,1:length(outp), function(ii) outp[[ii]])
                    dat <- list()
                    dat$Pred <- outpt %>% sapply(function(ee) unlist(ee$Pred)) %>% t
                    dat$trainPred <- outpt %>% sapply(function(ee) unlist(ee$trainPred)) %>% t
                    dat$TPR <- outpt %>% sapply(function(ee) colMeans(ee$TPR)) %>% t
                    dat$Outc <- outpt %>% sapply(function(ee) unlist(ee$Outc)) %>% t
                    dat$trainOutc <- outpt %>% sapply(function(ee) unlist(ee$trainOutc)) %>% t
                    return(dat)
                  })

plotdat <- list()
plotdat$TPR <- out.dat %>% 
  lapply(function(ee) data.frame(ee$TPR)) %>% Reduce(f=bind_rows)
plotdat$Pred <- out.dat %>% 
  lapply(function(ee) data.frame(ee$Pred)) %>% Reduce(f=bind_rows)
plotdat$trainPred <- out.dat %>% 
  lapply(function(ee) data.frame(ee$trainPred)) %>% Reduce(f=bind_rows)
plotdat$Outc <- out.dat %>% 
  lapply(function(ee) data.frame(ee$Outc)) %>% Reduce(f=bind_rows)
plotdat$trainOutc <- out.dat %>% 
  lapply(function(ee) data.frame(ee$trainOutc)) %>% Reduce(f=bind_rows)


fctname <- c("PO2PLS", "O2PLS", "PLS", "PPLS", "SIFA")
lay <- rbind(c(1,1,1,3,3,3,3),c(2,2,2,4,4,4,4))

try(dev.off(),silent = T)
CairoPDF(file = "SIMUoutcome_plot", height=10, width=18)
p1a <- plotdat$TPR %>% gather(key="Method", value="TPR", all_of(fctname)) %>% 
  mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method,y=TPR)) + geom_boxplot() + theme_bw() + 
#  facet_grid(. ~ SampSz) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), 
        legend.title=element_text(face="bold"))
p1b <- plotdat$TPR %>% mutate(PPLS=PPLS-PO2PLS, PLS=PLS-PO2PLS, O2PLS=O2PLS-PO2PLS,SIFA=SIFA-PO2PLS) %>%
  gather(key="Method", value="TPR", all_of(fctname)) %>% 
  mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
  filter(Method!="PO2PLS") %>%
  filter(Method != "SIFA") %>% 
  #filter(Method != "PPLS") %>%
  ggplot(aes(x=Method,y=TPR)) + geom_boxplot() + theme_bw() + ylab("Diff in TPR") + geom_hline(yintercept=0) + 
#  facet_grid(. ~ SampSz) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
p2a <- plotdat$Pred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Test") %>% 
  bind_rows(plotdat$trainPred %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Train")) %>% 
  mutate(Type = factor(Type, labels = c("train", "test"), 
                       levels = c("Train", "Test"))) %>%
  gather(key="Method",value="RMSEP",1:6) %>% mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
  filter(Method != "TRuE") %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method, y=RMSEP)) + geom_boxplot(aes(col=Type)) + theme_bw() + 
#  facet_grid(Noise*Heter ~ SampSz) + 
  geom_hline(yintercept=median(plotdat$Pred$TRuE %>% sqrt), lwd=1) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))
p2b <- plotdat$Outc %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Test") %>% 
  bind_rows(plotdat$trainOutc %>% mutate_at(.vars = 1:6, .funs = sqrt) %>% mutate(Type="Train")) %>% 
  mutate(Type = factor(Type, labels = c("train", "test"), 
                       levels = c("Train", "Test"))) %>%
  gather(key="Method",value="RMSEP",1:6) %>% mutate(Method=factor(Method,labels=fctname,levels=fctname)) %>% 
  filter(Method != "TRuE") %>% 
  filter(Method != "SIFA") %>%
  ggplot(aes(x=Method, y=RMSEP)) + geom_boxplot(aes(col=Type)) + theme_bw() + 
  #  facet_grid(Noise*Heter ~ SampSz) + 
  geom_hline(yintercept=median(plotdat$Outc$TRuE %>% sqrt), lwd=1) + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), legend.title=element_text(face="bold"))

grid.arrange(grobs=list(p1a,p1b,p2a,p2b), layout_matrix=lay)
dev.off()
