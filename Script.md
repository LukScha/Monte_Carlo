---
title: "MCS-Case Study: CPPI-Strategy (Gruppe 24)"
author: "David Pitz (11840550), Paul Gassner (11814819), Maximilian Lang (11844866), Marlene Bala (11817168), Lukas Schafferhofer (11771039)"
date: "12 11 2020"
output: 
  pdf_document:
    toc: yes
    number_sections: TRUE
    citation_package: natbib
    #keep_tex: true
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
m <- 1/0.2
# m = multiplier
b <- 0.3
# b = maximale Risikoanteil
Voestalpine.v<-c(27.769463,26.715246,29.577293,28.870247,31.685719,32.066765,31.609516,33.091343,29.139774,26.704212,28.648314,26.791195,24.655727,21.085018,23.224834,25.577765,27.382690,26.778147,26.112719,27.417482,26.702175,27.719145,29.006104,32.178509,33.559963,35.278912,35.751396,33.208977,34.518433,36.273380,36.718868,38.572811,40.275539,39.864376,43.643383,45.098625,46.054928,48.304775,43.980625,39.333096,40.386414,42.465328,36.450333,39.610283,36.968590,37.695717,30.032194,28.080437,24.971020,26.674023,26.042572,25.899059,27.391581,22.359108,25.994736,23.000128,20.950001,21.080000,22.450001,23.719999)
length(Voestalpine.v)

Voestalpine_last_twelve.v <- c(28.080437,24.971020,26.674023,26.042572,25.899059,27.391581,22.359108,25.994736,23.000128,20.950001,21.080000,22.450001,23.719999)
length(Voestalpine_last_twelve.v)

Erwartung <- mean(Voestalpine.v)
Standardabweichung <- sd(Voestalpine.v)
print(Erwartung)
print(Standardabweichung)
```


```{r}
set.seed(1) # Ensuring path consistency in simulation studies
(CPPI.m<-matrix(0,13,11, byrow=TRUE))
colnames(CPPI.m) <- c("t", "T_t,T", "F_t", "C_t", "X_r,t", "X_f,t", "S_t", "TSR_t", "W_t", "delta%W_t", "delta%W_t/TSR_t")

CPPI.m
```



```{r}
for (j in 0:12) {
  CPPI.m[j+1,1] <- j
  CPPI.m[j+1,2] <- 1-(j/months)
  CPPI.m[j+1,3] <- F_T*exp(-r_f*CPPI.m[j+1,2])
  CPPI.m[j+1,7] <- Voestalpine_last_twelve.v[j+1]
  
}

CPPI.m[1,4] <- max(W_0-F_T, 0)
CPPI.m[1,5] <- min(m*CPPI.m[1,4], b*W_0)
CPPI.m[1,6] <- W_0-CPPI.m[1,5]


for (j in 1:12){
  # TSR
  CPPI.m[j+1,8] <- log(CPPI.m[j+1,7]/CPPI.m[j-1+1,7])
  # W_t
  CPPI.m[j+1,9] <- CPPI.m[j,5]*exp(CPPI.m[j+1,8]) + CPPI.m[j,6]*exp((1/12)*r_f)
  #C_t
  CPPI.m[j+1,4] <- max(CPPI.m[j+1,9]-F_T, 0)
  # X_r,t
  CPPI.m[j+1,5] <- min(m*CPPI.m[j+1,4], b*CPPI.m[j+1,9])
  #X_f,t
  CPPI.m[j+1,6] <- CPPI.m[j+1,9]-CPPI.m[j+1,5]
  #delta%W_t
  CPPI.m[j,10] <- (CPPI.m[j+1,9]-CPPI.m[j,9])/CPPI.m[j,9]
  #delta%W_t/TSR_t
  CPPI.m[j,11] <- (CPPI.m[j+1,10]/CPPI.m[j+1,8])
}


CPPI.m
```
