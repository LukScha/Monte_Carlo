---
title: '\begin{center} Monte Carlo-Simulation Study: CPPI-Strategy \\ Case study teaching: Cheat sheet style \\ (Gruppe 24) \end{center}'
author: 
- David Pitz (11840550) 
- Paul Gassner (11814819) 
- Maximilian Lang (11844866) 
- Marlene Bala (11817168)  
- Lukas Schafferhofer (11771039) 
date: "12.11.2020"
output: 
  pdf_document:
    toc: yes
    number_sections: TRUE
    #keep_tex: TRUE
    citation_package: natbib
bibliography: literature.bib
biblio-style: humannat
include-before:
  \pagebreak
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
W_0 <- 1000000
F_T <- 0.65*W_0
R_f <- 0.01
r_f <- log(1+R_f)
T_T <- 1
months <- 12
m = 1
# m = multiplier
b = 0.1
# b = maximale Risikoanteil
Voestalpine.v<-c(27.769463,26.715246,29.577293,28.870247,31.685719,32.066765,31.609516,
                 33.091343,29.139774,26.704212,28.648314,26.791195,24.655727,21.085018,
                 23.224834,25.577765,27.382690,26.778147,26.112719,27.417482,26.702175,
                 27.719145,29.006104,32.178509,33.559963,35.278912,35.751396,33.208977,
                 34.518433,36.273380,36.718868,38.572811,40.275539,39.864376,43.643383,
                 45.098625,46.054928,48.304775,43.980625,39.333096,40.386414,42.465328,
                 36.450333,39.610283,36.968590,37.695717,30.032194,28.080437,24.971020,
                 26.674023,26.042572,25.899059,27.391581,22.359108,25.994736,23.000128,
                 20.950001,21.080000,22.450001,23.719999)


Voestalpine_last_twelve.v <- c(28.080437,24.971020,26.674023,26.042572,25.899059,
                               27.391581,22.359108,25.994736,23.000128,20.950001,
                               21.080000,22.450001,23.719999)


#Unterprüfung
#W_0 <- 1000
#F_T <- W_0
#R_f <- 0.05
#r_f <- log(1+R_f)
#T_T <- 1
#months <- 12
#m <- 15
# m = multiplier
#b <- 0.3
# b = maximale Risikoanteil
#PROBE.v <- c(1000, 968.51, 982.92, 940.53, 1035.58, 1059.87, 1015.06, 1048.39, 
             #1098.61, 1140.49, 1125.24, 1232.99, 1266.32)
```


```{r}
set.seed(1) # Ensuring path consistency in simulation studies
CPPI.m<-matrix(0,13,11, byrow=TRUE)
colnames(CPPI.m) <- c("t", "T_t,T", "F_t", "C_t", "X_r,t", "X_f,t", "S_t",
                      "TSR_t", "W_t", "delta%W_t", "delta%W_t/TSR_t")
```


```{r}
func_F_t <- function(arg_F_T, arg_r_f, arg_T_t) {
	return(arg_F_T * exp((-1)*arg_r_f*arg_T_t))}

func_C_t <- function(arg_W_t, arg_F_t) {
	return(max((arg_W_t - arg_F_t),0)) 
}

func_X_rt <- function(arg_m, arg_C_t, arg_b, arg_W_t) {
  if (arg_m*arg_C_t > arg_b*arg_W_t) return (arg_b*arg_W_t) else return (arg_m*arg_C_t)
}

func_X_ft <- function(arg_W_t, arg_X_rt) {
	return (arg_W_t - arg_X_rt)
}

func_TSR_t <- function(arg_S_t, arg_S_t1) {
	return (log(arg_S_t/arg_S_t1))
}

func_W_t <- function(arg_X_rt, arg_TSR, arg_X_ft, arg_r_f) {
	return ((arg_X_rt*exp(arg_TSR)) + (arg_X_ft*exp(arg_r_f/12)))
}

func_delta_Wt <- function(arg_W_t1, arg_W_t0) {
	return ((arg_W_t1 - arg_W_t0)/arg_W_t0)
}
for (j in 0:12) {
  #t
  CPPI.m[j+1,1] <- j
  #T,t,T
  CPPI.m[j+1,2] <- 1-(j/months)
  #F,t
  CPPI.m[j+1,3] <- func_F_t(F_T, r_f, CPPI.m[j+1,2])
  #S,t
  CPPI.m[j+1,7] <- Voestalpine_last_twelve.v[j+1]
}

CPPI.m[1,4] <- max(W_0-CPPI.m[1,3], 0)
CPPI.m[1,5] <- min(m*CPPI.m[1,4], b*W_0)
CPPI.m[1,6] <- W_0-CPPI.m[1,5]

for (j in 1:12){
  # TSR
  CPPI.m[j+1,8] <- func_TSR_t(CPPI.m[j+1,7],CPPI.m[j,7])
  # W_t
  CPPI.m[j+1,9] <- func_W_t(CPPI.m[j,5], CPPI.m[j+1,8], CPPI.m[j,6], r_f)
  #C_t
  CPPI.m[j+1,4] <- func_C_t(CPPI.m[j+1,9], CPPI.m[j+1,3])
  # X_r,t
  CPPI.m[j+1,5] <- func_X_rt(m, CPPI.m[j+1,4], b, CPPI.m[j+1,9])
  #X_f,t
  CPPI.m[j+1,6] <- func_X_ft(CPPI.m[j+1,9], CPPI.m[j+1,5])
  #delta%W_t
  CPPI.m[j+1,10] <- func_delta_Wt(CPPI.m[j+1,9], CPPI.m[j,9])
  #delta%W_t/TSR_t
  CPPI.m[j,11] <- (CPPI.m[j,10]/CPPI.m[j,8])
}

#delta%W_0
CPPI.m[1,10] <- 0
#delta%W_1
CPPI.m[2,10] <- (CPPI.m[2,9]-W_0)/W_0
#delta%W_t/TSR_1
CPPI.m[2,11] <- (CPPI.m[2,10]/CPPI.m[2,8])
#delta%W_t/TSR_12
CPPI.m[13,11] <- (CPPI.m[13,10]/CPPI.m[13,8])


#monatlichen Momente des RENDITEVEKTORS
mu <- 0
sig <- 0

#Erstellung des Renditevektors
Renditenvektor.v <- c()
for (i in 1:length(Voestalpine.v)-1){
 Renditenvektor.v[i] <- func_TSR_t(Voestalpine.v[i+1],Voestalpine.v[i])
}

Renditenvektor.v
mu <- mean(Renditenvektor.v)
Varianz <-0
for (i in 1:length(Renditenvektor.v)){
  Varianz <- Varianz + (1/(length(Renditenvektor.v-1)))*(Renditenvektor.v[i]-mu)^2
}
sig <- Varianz^(1/2)

cat("Erwartungswert der Rendite: ", mu, "\n")

cat("Standardabweichung der Rendite: ", sig, "\n")
```
```{r}
#begin Monte_Carlo

#simulation runs
sim = 1000
#dimensions
subp = 12

set.seed(1) # So wird jedesmal versichert dass immer die gleichen Zahlen gewählt werden
zufällige_Rendite.m<-matrix(rnorm(sim*subp),sim,subp, byrow=TRUE)


#Normalverteilte Renditen
r.m <- (mu-(sig^2))/2+sig*zufällige_Rendite.m
#Restlaufzeit
T_t_T.v <- c()
#Floor
F_t.v <- c()
  
for (j in 0:13) {
  #T,t,T #Restlaufzeit
  T_t_T.v[j] <- 1-((j-1)/months)
  
  #F,t
  F_t.v[j] <- func_F_t(F_T, r_f, T_t_T.v[j])
}

#Risikoparamter m und b deklarieren
multi.v <- c(1,3,5)
b_risk.v <- c(0.05, 0.1, 0.2, 0.25, 0.5)
Erwartungswert_Portfolio <- matrix(0,length(multi.v),length(b_risk.v))
colnames(Erwartungswert_Portfolio) <- c("0.05", "0.1", "0.2", "0.25", "0.5")
rownames(Erwartungswert_Portfolio) <- c("m=1 ","m=3 ","m=5 ")
Standardabweichung <- matrix(0,length(multi.v),length(b_risk.v))
colnames(Standardabweichung) <- c("0.05", "0.1", "0.2", "0.25", "0.5")
rownames(Standardabweichung) <- c("m=1 ","m=3 ","m=5 ")


#Iteration über alle Werte m
for(j in 1:length(multi.v)){
  
  #Iteration über alle Werte b
  for(i in 1:length(b_risk.v)){
    
    #notwendige Variablen
    C_t_sim <- matrix(0,sim, subp+1)
    X_r_t_sim <- matrix(0,sim,subp+1)
    X_f_t_sim <- matrix(0,sim,subp+1)
    W_t_sim <- matrix(0,sim,subp+1)
    #Summe aller Portfoliowerte am Zeitpunkt t=12
    sum_of_results <- 0
  
    #Iteration durch alle Simulationsläufe
    for(k in 1:sim){
      #Startwert W_0_sim
      W_t_sim[k,1] <- W_0
      C_t_sim[k,1] <- func_C_t(W_t_sim[k,1], F_t.v[1])
      X_r_t_sim[k,1] <- func_X_rt(multi.v[j], C_t_sim[k,1], b_risk.v[i],  W_t_sim[k,1])
      X_f_t_sim[k,1] <- func_X_ft(W_t_sim[k,1], X_r_t_sim[k,1])
      
      #Iteration durch alle Zeitpunkte
      for(l in 1:subp){
        
        W_t_sim[k,l+1] <- func_W_t(X_r_t_sim[k,l], r.m[k,l], X_f_t_sim[k,l], r_f)
        C_t_sim[k,l+1] <- func_C_t(W_t_sim[k,l+1], F_t.v[l+1])
        X_r_t_sim[k,l+1] <- func_X_rt(multi.v[j], C_t_sim[k,l+1], b_risk.v[i],  
                                      W_t_sim[k,l+1])
        X_f_t_sim[k,l+1] <- func_X_ft(W_t_sim[k,l+1], X_r_t_sim[k,l+1])
      }
      sum_of_results <- sum_of_results + W_t_sim[k,subp+1]
    }
    
    Mittwelwert_Gesamtwert <- sum_of_results/sim
    
   
    Erwartungswert_Portfolio[j,i] <- sum_of_results/sim
    
    Varianz_Gesamtwert <- 0
    for (x in 1:1000){
      Varianz_Gesamtwert <-Varianz_Gesamtwert + 
        ((W_t_sim[x,subp+1]-Mittwelwert_Gesamtwert)^2)
    }
    Varianz_Gesamtwert <- Varianz_Gesamtwert/1000
    sig_Gesamtwert <- Varianz_Gesamtwert^(1/2)

    Standardabweichung[j,i] <- sig_Gesamtwert

    #Teilen durch Million um Lesbarkeit zu verbessern
    hist((W_t_sim[,13]/1000000), 
     main = paste("m: " , multi.v[j], "b: ", b_risk.v[i]),
     xlab="Gesamtwert des Portfolios / 1e+06", 
     border="Black", 
     col="orange", 
     breaks=30
     )
    
  }
}

print(Erwartungswert_Portfolio)
print(Standardabweichung)
```










