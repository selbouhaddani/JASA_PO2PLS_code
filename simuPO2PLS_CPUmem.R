

## Organization: first memory and cpu times are calculated, then the results are plotted
## HOW TO run this script:
## Source the script, which generates **.out files and a SIMUcpumem.pdf plot
## Please clean up old **.out files before runnning the script
## TIME: about 20 minutes


library(tidyverse)
library(profmem)
library(pryr)
library(OmicsPLS)
library(PO2PLS)

#### About 20 min

out_memtime <- lapply(2:4, function(p_i){

  print(Sys.time())
  cat("\nNow starting power:", p_i, "\n")
  N <- 100
  p <- 1*10^p_i
  q <- p/2
  
  r <- 5
  rx = ry = 5
  
  
  X <- matrix(rnorm(N*p),nrow=N)
  Y <- matrix(rnorm(N*q),nrow=N)
  
  cat("\ninitial fits", "\n")
  
  fit0 = o2m(X,Y,r,0,0,p_thresh = 1,q_thresh = 1)
  fit1 = o2m(X,Y,r,rx,ry,p_thresh = 1,q_thresh = 1)
  fit2 = PO2PLS(X,Y,r,0,0,steps=1e0, init_param="r")
  fit3 = PO2PLS(X,Y,r,rx,ry,steps=1e0, init_param="r")
  
  cat("\nstarting out_mem", "\n")
  out_mem <- sapply(1:20, function(ii){
  
    cat("\n  -- starting PLS iteration", ii, "\n")
    intv <- ifelse(p_i == 2, 0.0005, 0.002)
    #o2m(X,Y,r,0,0)
    Rprof(filename="PLS.out", memory.profiling = TRUE, interval = intv)
    #p0 <- profmem({
    fit0 = suppressMessages(o2m(X,Y,r,0,0,p_thresh = 1,q_thresh = 1, max_iterations = 1))
    #})
    Rprof(NULL)
    
    cat("\n  -- starting O2PLS iteration", ii, "\n")
    #o2m(X,Y,r,rx,ry)
    Rprof(filename="O2PLS.out", memory.profiling = TRUE, interval = intv)
    #p1 <- profmem({
    fit1 = suppressMessages(o2m(X,Y,r,rx,ry,p_thresh = 1,q_thresh = 1, max_iterations = 1))
    #})
    Rprof(NULL)
    
    cat("\n  -- starting PPLS iteration", ii, "\n")
    #fit2 = PO2PLS(X,Y,r,0,0,steps=1e1)
    Rprof(filename="PPLS.out", memory.profiling = TRUE, interval = intv)
    #p2 <- profmem({
    fit2 = suppressMessages(PO2PLS(X,Y,r,0,0,steps=1e0, init_param=fit2$par))
    #})
    Rprof(NULL)
    
    cat("\n  -- starting PO2PLS iteration", ii, "\n")
    #fit3 = PO2PLS(X,Y,r,rx,ry,steps=1e1)
    Rprof(filename="PO2PLS.out", memory.profiling = TRUE, interval = intv)
    #p3 <- profmem({
    fit3 = suppressMessages(PO2PLS(X,Y,r,rx,ry,steps=1e0, init_param=fit3$par))
    #})
    Rprof(NULL)
    
    
    #fit4 = PO2PLS(X,Y,r,rx,ry,steps=1e1,homogen_joint = T)
    #Rprof(filename="SIFA.out", memory.profiling = TRUE, interval = 0.002)
    #p4 <- profmem({
    #fit4 = suppressMessages(PO2PLS(X,Y,r,rx,ry,steps=1e0, init_param=fit4$par, homogen_joint = T))
    #})
    #Rprof(NULL)
    
    cat("\n  -- collecting mem usage of iteration", ii, "\n")
    return(sapply(list.files(pattern=".out"), function(e) summaryRprof(e, memory="both")$by.total[1,3]))
    #cat("\n  -- ending iteration", ii)
  
  })
  
  cat("\nstarting out_time", "\n")
  out_time <- sapply(1:20, function(ii){
    
    cat("\n  -- starting PLS iteration", ii, "\n")
    #o2m(X,Y,r,0,0)
    Rprof(filename="PLS.out", memory.profiling = F, interval = 0.002)
    #p0 <- profmem({
    fit0 = suppressMessages(o2m(X,Y,r,0,0,p_thresh = 1,q_thresh = 1))
    #})
    Rprof(NULL)
    
    cat("\n  -- starting O2PLS iteration", ii, "\n")
    #o2m(X,Y,r,rx,ry)
    Rprof(filename="O2PLS.out", memory.profiling = F, interval = 0.002)
    #p1 <- profmem({
    fit1 = suppressMessages(o2m(X,Y,r,rx,ry,p_thresh = 1,q_thresh = 1))
    #})
    Rprof(NULL)
    
    cat("\n  -- starting PPLS iteration", ii, "\n")
    #fit2 = PO2PLS(X,Y,r,0,0,steps=1e1)
    Rprof(filename="PPLS.out", memory.profiling = F, interval = 0.02)
    #p2 <- profmem({
    fit2 = suppressMessages(PO2PLS(X,Y,r,0,0,steps=1e2, init_param=fit2$par,verbose = F))
    #})
    Rprof(NULL)
    
    cat("\n  -- starting PO2PLS iteration", ii, "\n")
    #fit3 = PO2PLS(X,Y,r,rx,ry,steps=1e1)
    Rprof(filename="PO2PLS.out", memory.profiling = F, interval = 0.02)
    #p3 <- profmem({
    fit3 = suppressMessages(PO2PLS(X,Y,r,rx,ry,steps=1e2, init_param=fit3$par,verbose = F))
    #})
    Rprof(NULL)
    
    
    #fit4 = PO2PLS(X,Y,r,rx,ry,steps=1e1,homogen_joint = T)
    #Rprof(filename="SIFA.out", memory.profiling = TRUE, interval = 0.02)
    #p4 <- profmem({
    #fit4 = suppressMessages(PO2PLS(X,Y,r,rx,ry,steps=1e2, init_param=fit4$par,verbose = F, homogen_joint = T))
    #})
    #Rprof(NULL)
    
    cat("\n  -- collecting time of iteration", ii, "\n")
    return(sapply(list.files(pattern=".out"), function(e) summaryRprof(e)$by.total[1,1]))
    #cat("\n  -- ending iteration", ii, "\n")
    
  })
  
  cat("\nend: success", "\n")
  print(Sys.time())
  return(list(mem=out_mem, time=out_time))

})

save(out_memtime, file = "out_memtime2.RData")

p1 <- lapply(out_memtime, function(e) {
  e$mem %>% t %>% as_tibble %>% gather(key="Method", value="Memory") %>% 
    mutate(Method = str_sub(Method, end=-5))
}) %>% Reduce(f=bind_rows) %>% bind_cols(Dim=rep(paste0("p=10^",2:4), each=4*20)) %>% 
  mutate(Method = factor(Method, levels=c("PO2PLS", "PPLS", "PLS", "O2PLS"), 
                         labels=c("PO2PLS", "PPLS", "PLS", "O2PLS"))) %>% 
  ggplot(aes(x=Dim, y=log10(Memory), col=Method)) + geom_boxplot() + 
  theme_bw() + xlab(NA) + ylab("log10 Memory usage (Mb)") + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), 
        axis.title.x = element_blank(), legend.title=element_text(face="bold"))

p2 <- lapply(out_memtime, function(e) {
  e$time %>% t %>% as_tibble %>% gather(key="Method", value="Time") %>% 
    mutate(Method = str_sub(Method, end=-5))
}) %>% Reduce(f=bind_rows) %>% bind_cols(Dim=rep(paste0("p=10^",2:4), each=4*20)) %>% 
  mutate(Method = factor(Method, levels=c("PO2PLS", "PPLS", "PLS", "O2PLS"), 
                         labels=c("PO2PLS", "PPLS", "PLS", "O2PLS"))) %>% 
  ggplot(aes(x=Dim, y=log10(Time), col=Method)) + geom_boxplot() + 
  #coord_cartesian(ylim = c(-1.65, 1.5)) + 
  theme_bw() + xlab("Dimensionality") + ylab("log10 Time (s)") + 
  theme(text=element_text(size=14), axis.title=element_text(face="bold"), 
        legend.title=element_text(face="bold"))

try(dev.off(),silent = T)
CairoPDF(file = "SIMUcpumem.pdf", height=7, width=9)
gridExtra::grid.arrange(p1, p2, nrow=2)
dev.off()
