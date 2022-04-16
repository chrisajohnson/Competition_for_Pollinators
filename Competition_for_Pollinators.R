###############################################################################################
#### This R script estimates the parameters of the annual plant population model (Eq. 1);  ####
#### quantifies the niche and average fitness differences and propagates their associated  ####
#### uncertainty (Fig. 2); compares parameter estimates between the pollination treatments ####
#### (Fig 3), and evaluates the effects of the pollinator decline experiment (Fig. 4)      ####
###############################################################################################

# Load packages and set working directory
library(ggplot2)
library(plyr)
library(broom)
library(dotwhisker)
library(propagate)
library(nlme)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# READ DATA
data <- read.csv("Competition data.csv", header=TRUE)

# COMPETITION FUNCTIONS
Scomp <- function(S,B,P,C,N,r,aSS,aSB,aSP,aSC,aSN){r/(1+aSS*S+aSB*B+aSP*P+aSC*C+aSN*N)}
Bcomp <- function(S,B,P,C,N,r,aBS,aBB,aBP,aBC,aBN){r/(1+aBS*S+aBB*B+aBP*P+aBC*C+aBN*N)}
Pcomp <- function(S,B,P,C,N,r,aPS,aPB,aPP,aPC,aPN){r/(1+aPS*S+aPB*B+aPP*P+aPC*C+aPN*N)}
Ccomp <- function(S,B,P,C,N,r,aCS,aCB,aCP,aCC,aCN){r/(1+aCS*S+aCB*B+aCP*P+aCC*C+aCN*N)}
Ncomp <- function(S,B,P,C,N,r,aNS,aNB,aNP,aNC,aNN){r/(1+aNS*S+aNB*B+aNP*P+aNC*C+aNN*N)}