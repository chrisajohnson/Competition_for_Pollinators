###############################################################################################
#### This R script estimates the parameters of the annual plant population model (Eq. 1);  ####
#### quantifies the niche and average fitness differences and propagates their associated  ####
#### uncertainty (Fig. 2); compares parameter estimates between the pollination treatments ####
#### (Fig 3), and evaluates the effects of the pollinator decline experiment (Fig. 4)      ####
###############################################################################################

# IMPORTANT: nls works on R v. 3.3 and v. 3.6.3, but does not work on R v. 4.0

# INSTALL AND LOAD PACKAGES AND SET WORKING DIRECTORY
#install.packages("tmvtnorm", dependencies = FALSE) # install if necessary
#install.packages("minpack.lm", dependencies = FALSE) # install if necessary
library(propagate)
library(broom)
library(tidyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# READ DATA AND SORT BY FOCAL SPECIES
data <- read.csv("Competition data.csv", header=TRUE)
s.data <- data[data$Focal == 'S',]
b.data <- data[data$Focal == 'B',]
p.data <- data[data$Focal == 'P',]
c.data <- data[data$Focal == 'C',]
n.data <- data[data$Focal == 'N',]


###############################################################################################################################################
# INSTRUMENT VARIABLE ANALYSIS  (Extended Data Table 7)
# Step 1: regress neighbor density on sowing density
# Sinapis
summary(lm(Neighbors ~ 0 + S + B + P + C + N, data = s.data))
# Buglossoides
summary(lm(Neighbors ~ 0 + S + B + P + C + N, data = b.data))
# Papaver
summary(lm(Neighbors ~ 0 + S + B + P + C + N, data = p.data))
# Centaurea
summary(lm(Neighbors ~ 0 + S + B + P + C + N, data = c.data))
# Nigella
summary(lm(Neighbors ~ 0 + S + B + P + C + N, data = n.data))


# Step 2: regress per capita seed production on the *predicted* neighbor density from step 1
# Sinapis
# Full model with separate parameter estimates for each treatment
#          (parameter 1 (e.g., aSS1) is for control ambient pollination treatment, parameter 2 (e.g., aSS2) is for hand pollination treatment)
S.full <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aSS[Tr]*S_predict+aSB[Tr]*B_predict+aSP[Tr]*P_predict+aSC[Tr]*C_predict+aSN[Tr]*N_predict)+1),
                 data=s.data, start=list(lambda=c(30000,30000), aSS=c(1,1), aSB=c(1,1), aSP=c(1,1), aSC=c(1,1), aSN=c(0.1,0.1)))
summary(S.full)
# Reduced model with one estimate of each parameter across both treatments
S.reduced <- nls(log(Seeds+1)~log(lambda/(1+aSS*S_predict+aSB*B_predict+aSP*P_predict+aSC*C_predict+aSN*N_predict)+1),
                      data=s.data, start=list(lambda=30000, aSS=1, aSB=1, aSP=1, aSC=1, aSN=0.1))
summary(S.reduced)

# Buglossoides
# Full model with separate parameter estimates for each treatment
#          (parameter 1 (e.g., aBS1) is for control ambient pollination treatment, parameter 2 (e.g., aBS2) is for hand pollination treatment)
B.full <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aBS[Tr]*S_predict+aBB[Tr]*B_predict+aBP[Tr]*P_predict+aBC[Tr]*C_predict+aBN[Tr]*N_predict)+1),
                 data=b.data, start=list(lambda=c(3000,3000), aBS=c(1,1), aBB=c(1,1), aBP=c(1,1), aBC=c(1,1), aBN=c(0.01,0.01)))
summary(B.full)
# Reduced model with one estimate of each parameter across both treatments
B.reduced <- nls(log(Seeds+1)~log(lambda/(1+aBS*S_predict+aBB*B_predict+aBP*P_predict+aBC*C_predict+aBN*N_predict)+1),
                      data=b.data, start=list(lambda=3000, aBS=1, aBB=1, aBP=1, aBC=1, aBN=0.01))
summary(B.reduced)

# Papaver
# Full model with separate parameter estimates for each treatment
#          (parameter 1 (e.g., aPS1) is for control ambient pollination treatment, parameter 2 (e.g., aPS2) is for hand pollination treatment)
P.full <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aPS[Tr]*S_predict+aPB[Tr]*B_predict+aPP[Tr]*P_predict+aPC[Tr]*C_predict+aPN[Tr]*N_predict)+1),
                 data=p.data, start=list(lambda=c(30000,30000), aPS=c(1,1), aPB=c(1,1), aPP=c(1,1), aPC=c(1,1), aPN=c(0.1,0.1)))
summary(P.full)
# Reduced model with one estimate of each parameter across both treatments
P.reduced <- nls(log(Seeds+1)~log(lambda/(1+aPS*S_predict+aPB*B_predict+aPP*P_predict+aPC*C_predict+aPN*N_predict)+1),
                      data=p.data, start=list(lambda=30000, aPS=1, aPB=1, aPP=1, aPC=1, aPN=0.1))
summary(P.reduced)

# Centaurea
# Full model with separate parameter estimates for each treatment
#          (parameter 1 (e.g., aCS1) is for control ambient pollination treatment, parameter 2 (e.g., aCS2) is for hand pollination treatment)
C.full <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aCS[Tr]*S_predict+aCB[Tr]*B_predict+aCP[Tr]*P_predict+aCC[Tr]*C_predict+aCN[Tr]*N_predict)+1),
                 data=c.data, start=list(lambda=c(1000,1000), aCS=c(0.1,0.1), aCB=c(1,1), aCP=c(0.1,0.1), aCC=c(0.1,0.1), aCN=c(0.1,0.1)))
summary(C.full)
# Reduced model with one estimate of each parameter across both treatments
C.reduced <- nls(log(Seeds+1)~log(lambda/(1+aCS*S_predict+aCB*B_predict+aCP*P_predict+aCC*C_predict+aCN*N_predict)+1),
                      data=c.data, start=list(lambda=1000, aCS=0.1, aCB=1, aCP=0.1, aCC=0.1, aCN=0.1))
summary(C.reduced)

# Nigella
# Full model with separate parameter estimates for each treatment
#          (parameter 1 (e.g., aNS1) is for control ambient pollination treatment, parameter 2 (e.g., aNS2) is for hand pollination treatment)
N.full <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aNS[Tr]*S_predict+aNB[Tr]*B_predict+aNP[Tr]*P_predict+aNC[Tr]*C_predict+aNN[Tr]*N_predict)+1),
                 data=n.data, start=list(lambda=c(1000,1000), aNS=c(0.1,0.1), aNB=c(0.1,0.1), aNP=c(0.1,0.1), aNC=c(0.1,0.1), aNN=c(0.01,0.01)))
summary(N.full)
# Reduced model with one estimate of each parameter across both treatments
N.reduced <- nls(log(Seeds+1)~log(lambda/(1+aNS*S_predict+aNB*B_predict+aNP*P_predict+aNC*C_predict+aNN*N_predict)+1),
                      data=n.data, start=list(lambda=3000, aNS=0.1, aNB=0.1, aNP=0.1, aNC=0.1, aNN=0.01))
summary(N.reduced)



###############################################################################################################################################
# LIKELIHOOD-RATIO TESTS EVALUATING WHETHER PARAMETERS DIFFER BETWEEN POLLINATION TREATMENTS (Figure 3 and Extended Data Table 8)
# Sinapis
# Test if separate parameters are required for lambda
S.lam <- nls(log(Seeds+1)~log(lambda/(1+aSS[Tr]*S_predict+aSB[Tr]*B_predict+aSP[Tr]*P_predict+aSC[Tr]*C_predict+aSN[Tr]*N_predict)+1),
             data=s.data, start=list(lambda=30000, aSS=c(1,1), aSB=c(1,1), aSP=c(1,1), aSC=c(1,1), aSN=c(0.01,0.01)))
Q <- -2 * (logLik(S.lam) - logLik(S.full))
Q
df.Q <- df.residual(S.lam) - df.residual(S.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aSS 
S.ass <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aSS*S_predict+aSB[Tr]*B_predict+aSP[Tr]*P_predict+aSC[Tr]*C_predict+aSN[Tr]*N_predict)+1),
             data=s.data, start=list(lambda=c(30000,30000), aSS=1, aSB=c(1,1), aSP=c(1,1), aSC=c(1,1), aSN=c(0.01,0.01)))
Q <- -2 * (logLik(S.ass) - logLik(S.full))
Q
df.Q <- df.residual(S.ass) - df.residual(S.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aSB 
S.asb <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aSS[Tr]*S_predict+aSB*B_predict+aSP[Tr]*P_predict+aSC[Tr]*C_predict+aSN[Tr]*N_predict)+1),
             data=s.data, start=list(lambda=c(30000,30000), aSS=c(1,1), aSB=1, aSP=c(1,1), aSC=c(1,1), aSN=c(0.01,0.01)))
Q <- -2 * (logLik(S.asb) - logLik(S.full))
Q
df.Q <- df.residual(S.asb) - df.residual(S.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aSP 
S.asp <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aSS[Tr]*S_predict+aSB[Tr]*B_predict+aSP*P_predict+aSC[Tr]*C_predict+aSN[Tr]*N_predict)+1),
             data=s.data, start=list(lambda=c(30000,30000), aSS=c(1,1), aSB=c(1,1), aSP=1, aSC=c(1,1), aSN=c(0.01,0.01)))
Q <- -2 * (logLik(S.asp) - logLik(S.full))
Q
df.Q <- df.residual(S.asp) - df.residual(S.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aSC 
S.asc <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aSS[Tr]*S_predict+aSB[Tr]*B_predict+aSP[Tr]*P_predict+aSC*C_predict+aSN[Tr]*N_predict)+1),
             data=s.data, start=list(lambda=c(30000,30000), aSS=c(1,1), aSB=c(1,1), aSP=c(1,1), aSC=1, aSN=c(0.01,0.01)))
Q <- -2 * (logLik(S.asc) - logLik(S.full))
Q
df.Q <- df.residual(S.asc) - df.residual(S.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aSN
S.asn <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aSS[Tr]*S_predict+aSB[Tr]*B_predict+aSP[Tr]*P_predict+aSC[Tr]*C_predict+aSN*N_predict)+1),
             data=s.data, start=list(lambda=c(30000,30000), aSS=c(1,1), aSB=c(1,1), aSP=c(1,1), aSC=c(1,1), aSN=0.01))
Q <- -2 * (logLik(S.asn) - logLik(S.full))
Q
df.Q <- df.residual(S.asn) - df.residual(S.full)
1 - pchisq(Q, df.Q)

# Buglossoides
# Test if separate parameters are required for lambda 
B.lam <- nls(log(Seeds+1)~log(lambda/(1+aBS[Tr]*S_predict+aBB[Tr]*B_predict+aBP[Tr]*P_predict+aBC[Tr]*C_predict+aBN[Tr]*N_predict)+1),
             data=b.data, start=list(lambda=3000, aBS=c(1,1), aBB=c(1,1), aBP=c(1,1), aBC=c(1,1), aBN=c(0.01,0.01)))
Q <- -2 * (logLik(B.lam) - logLik(B.full))
Q
df.Q <- df.residual(B.lam) - df.residual(B.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aBS 
B.abs <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aBS*S_predict+aBB[Tr]*B_predict+aBP[Tr]*P_predict+aBC[Tr]*C_predict+aBN[Tr]*N_predict)+1),
             data=b.data, start=list(lambda=c(3000,3000), aBS=1, aBB=c(1,1), aBP=c(1,1), aBC=c(1,1), aBN=c(0.01,0.01)))
Q <- -2 * (logLik(B.abs) - logLik(B.full))
Q
df.Q <- df.residual(B.abs) - df.residual(B.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aBB 
B.abb <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aBS[Tr]*S_predict+aBB*B_predict+aBP[Tr]*P_predict+aBC[Tr]*C_predict+aBN[Tr]*N_predict)+1),
             data=b.data, start=list(lambda=c(3000,3000), aBS=c(1,1), aBB=1, aBP=c(1,1), aBC=c(1,1), aBN=c(0.01,0.01)))
Q <- -2 * (logLik(B.abb) - logLik(B.full))
Q
df.Q <- df.residual(B.abb) - df.residual(B.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aBP 
B.abp <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aBS[Tr]*S_predict+aBB[Tr]*B_predict+aBP*P_predict+aBC[Tr]*C_predict+aBN[Tr]*N_predict)+1),
             data=b.data, start=list(lambda=c(3000,3000), aBS=c(1,1), aBB=c(1,1), aBP=1, aBC=c(1,1), aBN=c(0.01,0.01)))
Q <- -2 * (logLik(B.abp) - logLik(B.full))
Q
df.Q <- df.residual(B.abp) - df.residual(B.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aBC 
B.abc <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aBS[Tr]*S_predict+aBB[Tr]*B_predict+aBP[Tr]*P_predict+aBC*C_predict+aBN[Tr]*N_predict)+1),
             data=b.data, start=list(lambda=c(4000,4000), aBS=c(0.1,0.1), aBB=c(1,1), aBP=c(1,1), aBC=1, aBN=c(0.01,0.01)))
Q <- -2 * (logLik(B.abc) - logLik(B.full))
Q
df.Q <- df.residual(B.abc) - df.residual(B.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aBN
B.abn <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aBS[Tr]*S_predict+aBB[Tr]*B_predict+aBP[Tr]*P_predict+aBC[Tr]*C_predict+aBN*N_predict)+1),
             data=b.data, start=list(lambda=c(3000,3000), aBS=c(1,1), aBB=c(1,1), aBP=c(1,1), aBC=c(1,1), aBN=0.01))
Q <- -2 * (logLik(B.abn) - logLik(B.full))
Q
df.Q <- df.residual(B.abn) - df.residual(B.full)
1 - pchisq(Q, df.Q)

# Papaver
# Test if separate parameters are required for lambda 
P.lam <- nls(log(Seeds+1)~log(lambda/(1+aPS[Tr]*S_predict+aPB[Tr]*B_predict+aPP[Tr]*P_predict+aPC[Tr]*C_predict+aPN[Tr]*N_predict)+1),
             data=p.data, start=list(lambda=50000, aPS=c(0.1,0.1), aPB=c(0.1,0.1), aPP=c(0.1,0.1), aPC=c(0.1,0.1), aPN=c(0.1,0.1)))
Q <- -2 * (logLik(P.lam) - logLik(P.full))
Q
df.Q <- df.residual(P.lam) - df.residual(P.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aPS 
P.aps <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aPS*S_predict+aPB[Tr]*B_predict+aPP[Tr]*P_predict+aPC[Tr]*C_predict+aPN[Tr]*N_predict)+1),
             data=p.data, start=list(lambda=c(50000,50000), aPS=1, aPB=c(1,1), aPP=c(1,1), aPC=c(1,1), aPN=c(0.1,0.1)))
Q <- -2 * (logLik(P.aps) - logLik(P.full))
Q
df.Q <- df.residual(P.aps) - df.residual(P.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aPB 
P.apb <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aPS[Tr]*S_predict+aPB*B_predict+aPP[Tr]*P_predict+aPC[Tr]*C_predict+aPN[Tr]*N_predict)+1),
             data=p.data, start=list(lambda=c(50000,50000), aPS=c(1,1), aPB=1, aPP=c(1,1), aPC=c(1,1), aPN=c(0.1,0.1)))
Q <- -2 * (logLik(P.apb) - logLik(P.full))
Q
df.Q <- df.residual(P.apb) - df.residual(P.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aPP 
P.app <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aPS[Tr]*S_predict+aPB[Tr]*B_predict+aPP*P_predict+aPC[Tr]*C_predict+aPN[Tr]*N_predict)+1),
             data=p.data, start=list(lambda=c(50000,50000), aPS=c(0.1,0.1), aPB=c(0.1,0.1), aPP=1, aPC=c(0.1,0.1), aPN=c(0.1,0.1)))
Q <- -2 * (logLik(P.app) - logLik(P.full))
Q
df.Q <- df.residual(P.app) - df.residual(P.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aPC 
P.apc <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aPS[Tr]*S_predict+aPB[Tr]*B_predict+aPP[Tr]*P_predict+aPC*C_predict+aPN[Tr]*N_predict)+1),
             data=p.data, start=list(lambda=c(50000,50000), aPS=c(0.1,0.1), aPB=c(1,1), aPP=c(1,1), aPC=1, aPN=c(0.1,0.1)))
Q <- -2 * (logLik(P.apc) - logLik(P.full))
Q
df.Q <- df.residual(P.apc) - df.residual(P.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aPN
P.apn <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aPS[Tr]*S_predict+aPB[Tr]*B_predict+aPP[Tr]*P_predict+aPC[Tr]*C_predict+aPN*N_predict)+1),
             data=p.data, start=list(lambda=c(50000,50000), aPS=c(1,1), aPB=c(1,1), aPP=c(1,1), aPC=c(1,1), aPN=0.1))
Q <- -2 * (logLik(P.apn) - logLik(P.full))
Q
df.Q <- df.residual(P.apn) - df.residual(P.full)
1 - pchisq(Q, df.Q)

# Centaurea
# Test if separate parameters are required for lambda 
C.lam <- nls(log(Seeds+1)~log(lambda/(1+aCS[Tr]*S_predict+aCB[Tr]*B_predict+aCP[Tr]*P_predict+aCC[Tr]*C_predict+aCN[Tr]*N_predict)+1),
             data=c.data, start=list(lambda=3000, aCS=c(0.1,0.1), aCB=c(0.1,0.1), aCP=c(0.1,0.1), aCC=c(0.1,0.1), aCN=c(0.1,0.1)))
Q <- -2 * (logLik(C.lam) - logLik(C.full))
Q
df.Q <- df.residual(C.lam) - df.residual(C.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aCS 
C.acs <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aCS*S_predict+aCB[Tr]*B_predict+aCP[Tr]*P_predict+aCC[Tr]*C_predict+aCN[Tr]*N_predict)+1),
             data=c.data, start=list(lambda=c(3000,3000), aCS=0.1, aCB=c(0.1,0.1), aCP=c(0.1,0.1), aCC=c(0.1,0.1), aCN=c(0.1,0.1)))
Q <- -2 * (logLik(C.acs) - logLik(C.full))
Q
df.Q <- df.residual(C.acs) - df.residual(C.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aCB 
C.acb <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aCS[Tr]*S_predict+aCB*B_predict+aCP[Tr]*P_predict+aCC[Tr]*C_predict+aCN[Tr]*N_predict)+1),
             data=c.data, start=list(lambda=c(3000,3000), aCS=c(0.1,0.1), aCB=0.1, aCP=c(0.1,0.1), aCC=c(0.1,0.1), aCN=c(0.1,0.1)))
Q <- -2 * (logLik(C.acb) - logLik(C.full))
Q
df.Q <- df.residual(C.acb) - df.residual(C.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aCP 
C.acp <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aCS[Tr]*S_predict+aCB[Tr]*B_predict+aCP*P_predict+aCC[Tr]*C_predict+aCN[Tr]*N_predict)+1),
             data=c.data, start=list(lambda=c(3000,3000), aCS=c(0.1,0.1), aCB=c(0.1,0.1), aCP=0.1, aCC=c(0.1,0.1), aCN=c(0.1,0.1)))
Q <- -2 * (logLik(C.acp) - logLik(C.full))
Q
df.Q <- df.residual(C.acp) - df.residual(C.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aCC 
C.acc <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aCS[Tr]*S_predict+aCB[Tr]*B_predict+aCP[Tr]*P_predict+aCC*C_predict+aCN[Tr]*N_predict)+1),
             data=c.data, start=list(lambda=c(1000,1000), aCS=c(0.1,0.1), aCB=c(0.1,0.1), aCP=c(0.1,0.1), aCC=0.1, aCN=c(0.1,0.1)))
Q <- -2 * (logLik(C.acc) - logLik(C.full))
Q
df.Q <- df.residual(C.acc) - df.residual(C.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aCN
C.acn <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aCS[Tr]*S_predict+aCB[Tr]*B_predict+aCP[Tr]*P_predict+aCC[Tr]*C_predict+aCN*N_predict)+1),
             data=c.data, start=list(lambda=c(3000,3000), aCS=c(0.1,0.1), aCB=c(0.1,0.1), aCP=c(0.1,0.1), aCC=c(0.1,0.1), aCN=0.1))
Q <- -2 * (logLik(C.acn) - logLik(C.full))
Q
df.Q <- df.residual(C.acn) - df.residual(C.full)
1 - pchisq(Q, df.Q)

# Nigella
# Test if separate parameters are required for lambda 
N.lam <- nls(log(Seeds+1)~log(lambda/(1+aNS[Tr]*S_predict+aNB[Tr]*B_predict+aNP[Tr]*P_predict+aNC[Tr]*C_predict+aNN[Tr]*N_predict)+1),
             data=n.data, start=list(lambda=3000, aNS=c(0.1,0.1), aNB=c(1,1), aNP=c(0.1,0.1), aNC=c(0.1,0.1), aNN=c(0.1,0.1)))
Q <- -2 * (logLik(N.lam) - logLik(N.full))
Q
df.Q <- df.residual(N.lam) - df.residual(N.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aNS 
N.ans <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aNS*S_predict+aNB[Tr]*B_predict+aNP[Tr]*P_predict+aNC[Tr]*C_predict+aNN[Tr]*N_predict)+1),
             data=n.data, start=list(lambda=c(1000,1000), aNS=0.1, aNB=c(1,1), aNP=c(0.1,0.1), aNC=c(0.1,0.1), aNN=c(0.1,0.1)))
Q <- -2 * (logLik(N.ans) - logLik(N.full))
Q
df.Q <- df.residual(N.ans) - df.residual(N.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aNB 
N.anb <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aNS[Tr]*S_predict+aNB*B_predict+aNP[Tr]*P_predict+aNC[Tr]*C_predict+aNN[Tr]*N_predict)+1),
             data=n.data, start=list(lambda=c(1000,1000), aNS=c(0.1,0.1), aNB=1, aNP=c(0.1,0.1), aNC=c(0.1,0.1), aNN=c(0.1,0.1)))
Q <- -2 * (logLik(N.anb) - logLik(N.full))
Q
df.Q <- df.residual(N.anb) - df.residual(N.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aNP 
N.anp <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aNS[Tr]*S_predict+aNB[Tr]*B_predict+aNP*P_predict+aNC[Tr]*C_predict+aNN[Tr]*N_predict)+1),
             data=n.data, start=list(lambda=c(1000,1000), aNS=c(0.1,0.1), aNB=c(1,1), aNP=0.1, aNC=c(0.1,0.1), aNN=c(0.1,0.1)))
Q <- -2 * (logLik(N.anp) - logLik(N.full))
Q
df.Q <- df.residual(N.anp) - df.residual(N.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aNC 
N.anc <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aNS[Tr]*S_predict+aNB[Tr]*B_predict+aNP[Tr]*P_predict+aNC*C_predict+aNN[Tr]*N_predict)+1),
             data=n.data, start=list(lambda=c(1000,1000), aNS=c(0.1,0.1), aNB=c(1,1), aNP=c(0.1,0.1), aNC=0.1, aNN=c(0.1,0.1)))
Q <- -2 * (logLik(N.anc) - logLik(N.full))
Q
df.Q <- df.residual(N.anc) - df.residual(N.full)
1 - pchisq(Q, df.Q)

# Test if separate parameters are required for aNN
N.ann <- nls(log(Seeds+1)~log(lambda[Tr]/(1+aNS[Tr]*S_predict+aNB[Tr]*B_predict+aNP[Tr]*P_predict+aNC[Tr]*C_predict+aNN*N_predict)+1),
             data=n.data, start=list(lambda=c(1000,1000), aNS=c(0.1,0.1), aNB=c(1,1), aNP=c(0.1,0.1), aNC=c(0.1,0.1), aNN=0.1))
Q <- -2 * (logLik(N.ann) - logLik(N.full))
Q
df.Q <- df.residual(N.ann) - df.residual(N.full)
1 - pchisq(Q, df.Q)



###############################################################################################################################################
# QUANTIFY NICHE AND FITNESS DIFFERENCES (Figure 2)
# Set Parameters
# Sinapis
lambda.s.h<-coef(S.full)[2]
lambda.s.c<-coef(S.full)[1]
g.s<-0.021 # Germination fraction (Extended Data Table 4)
s.s<-0.774 # Seed bank survival (Extended Data Table 4)
ass.h<-coef(S.full)[4]
ass.c<-coef(S.full)[3]
asb.h<-coef(S.full)[6]
asb.c<-coef(S.full)[5]
asp.h<-coef(S.full)[8]
asp.c<-coef(S.full)[7]
asc.h<-coef(S.full)[10]
asc.c<-coef(S.full)[9]
asn.h<-coef(S.full)[12]
asn.c<-coef(S.full)[11]
s.par<-coef(summary(S.full))
# Buglossoides
lambda.b.h<-coef(B.full)[2]
lambda.b.c<-coef(B.full)[1]
g.b<-0.038 # Germination fraction (Extended Data Table 4)
s.b<-0.231 # Seed bank survival (Extended Data Table 4)
abs.h<-coef(B.full)[4]
abs.c<-coef(B.full)[3]
abb.h<-coef(B.full)[6]
abb.c<-coef(B.full)[5]
abp.h<-coef(B.full)[8]
abp.c<-coef(B.full)[7]
abc.h<-coef(B.full)[10]
abc.c<-coef(B.full)[9]
abn.h<-coef(B.full)[12]
abn.c<-coef(B.full)[11]
b.par<-coef(summary(B.full))
# Papaver
lambda.p.h<-coef(P.full)[2]
lambda.p.c<-coef(P.full)[1]
g.p<-0.005 # Germination fraction (Extended Data Table 4)
s.p<-0.166 # Seed bank survival (Extended Data Table 4)
aps.h<-coef(P.full)[4]
aps.c<-coef(P.full)[3]
apb.h<-coef(P.full)[6]
apb.c<-coef(P.full)[5]
app.h<-coef(P.full)[8]
app.c<-coef(P.full)[7]
apc.h<-coef(P.full)[10]
apc.c<-coef(P.full)[9]
apn.h<-coef(P.full)[12]
apn.c<-coef(P.full)[11]
p.par<-coef(summary(P.full))
# Centaurea
lambda.c.h<-coef(C.full)[2]
lambda.c.c<-coef(C.full)[1]
g.c<-0.021 # Germination fraction (Extended Data Table 4)
s.c<-0.413 # Seed bank survival (Extended Data Table 4)
acs.h<-coef(C.full)[4]
acs.c<-coef(C.full)[3]
acb.h<-coef(C.full)[6]
acb.c<-coef(C.full)[5]
acp.h<-coef(C.full)[8]
acp.c<-coef(C.full)[7]
acc.h<-coef(C.full)[10]
acc.c<-coef(C.full)[9]
acn.h<-coef(C.full)[12]
acn.c<-coef(C.full)[11]
c.par<-coef(summary(C.full))
# Nigella
lambda.n.h<-coef(N.full)[2]
lambda.n.c<-coef(N.full)[1]
g.n<-0.008 # Germination fraction (Extended Data Table 4)
s.n<-0.084 # Seed bank survival (Extended Data Table 4)
ans.h<-coef(N.full)[4]
ans.c<-coef(N.full)[3]
anb.h<-coef(N.full)[6]
anb.c<-coef(N.full)[5]
anp.h<-coef(N.full)[8]
anp.c<-coef(N.full)[7]
anc.h<-coef(N.full)[10]
anc.c<-coef(N.full)[9]
ann.h<-coef(N.full)[12]
ann.c<-coef(N.full)[11]
n.par<-coef(summary(N.full))

# Calculate niche overlap (rho)
rho.sb.h <- as.numeric(sqrt((asb.h*abs.h)/(ass.h*abb.h)))
rho.sb.c <- as.numeric(sqrt((asb.c*abs.c)/(ass.c*abb.c)))
rho.sp.h <- as.numeric(sqrt((asp.h*aps.h)/(ass.h*app.h)))
rho.sp.c <- as.numeric(sqrt((asp.c*aps.c)/(ass.c*app.c)))
rho.sc.h <- as.numeric(sqrt((asc.h*acs.h)/(ass.h*acc.h)))
rho.sc.c <- as.numeric(sqrt((asc.c*acs.c)/(ass.c*acc.c)))
rho.sn.h <- as.numeric(sqrt((asn.h*ans.h)/(ass.h*ann.h)))
rho.sn.c <- as.numeric(sqrt((asn.c*ans.c)/(ass.c*ann.c)))
rho.bp.h <- as.numeric(sqrt((abp.h*apb.h)/(abb.h*app.h)))
rho.bp.c <- as.numeric(sqrt((abp.c*apb.c)/(abb.c*app.c)))
rho.bc.h <- as.numeric(sqrt((abc.h*acb.h)/(abb.h*acc.h)))
rho.bc.c <- as.numeric(sqrt((abc.c*acb.c)/(abb.c*acc.c)))
rho.bn.h <- as.numeric(sqrt((abn.h*anb.h)/(abb.h*ann.h)))
rho.bn.c <- as.numeric(sqrt((abn.c*anb.c)/(abb.c*ann.c)))
rho.pc.h <- as.numeric(sqrt((apc.h*acp.h)/(app.h*acc.h)))
rho.pc.c <- as.numeric(sqrt((apc.c*acp.c)/(app.c*acc.c)))
rho.pn.h <- as.numeric(sqrt((apn.h*anp.h)/(app.h*ann.h)))
rho.pn.c <- as.numeric(sqrt((apn.c*anp.c)/(app.c*ann.c)))
rho.cn.h <- as.numeric(sqrt((acn.h*anc.h)/(acc.h*ann.h)))
rho.cn.c <- as.numeric(sqrt((acn.c*anc.c)/(acc.c*ann.c)))

# Calculate niche difference (-ln[rho])
-log(rho.sb.h)
-log(rho.sb.c)
-log(rho.sp.h)
-log(rho.sp.c)
-log(rho.sc.h)
-log(rho.sc.c)
-log(rho.sn.h)
-log(rho.sn.c)
-log(rho.bp.h)
-log(rho.bp.c)
-log(rho.bc.h)
-log(rho.bc.c)
-log(rho.bn.h)
-log(rho.bn.c)
-log(rho.pc.h)
-log(rho.pc.c)
-log(rho.pn.h)
-log(rho.pn.c)
-log(rho.cn.h)
-log(rho.cn.c)

# Calculate fitness differences
k.sb.h <- as.numeric((lambda.b.h*g.b/(1-(1-g.b)*s.b)-1)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asb.h*ass.h)/(abs.h*abb.h)))
k.sb.c <- as.numeric((lambda.b.c*g.b/(1-(1-g.b)*s.b)-1)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asb.c*ass.c)/(abs.c*abb.c)))
k.sp.h <- as.numeric((lambda.p.h*g.p/(1-(1-g.p)*s.p)-1)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asp.h*ass.h)/(aps.h*app.h)))
k.sp.c <- as.numeric((lambda.p.c*g.p/(1-(1-g.p)*s.p)-1)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asp.c*ass.c)/(aps.c*app.c)))
k.sc.h <- as.numeric((lambda.c.h*g.c/(1-(1-g.c)*s.c)-1)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asc.h*ass.h)/(acs.h*acc.h)))
k.sc.c <- as.numeric((lambda.c.c*g.c/(1-(1-g.c)*s.c)-1)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asc.c*ass.c)/(acs.c*acc.c)))
k.sn.h <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asn.h*ass.h)/(ans.h*ann.h)))
k.sn.c <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1)*sqrt((asn.c*ass.c)/(ans.c*ann.c)))
k.bp.h <- as.numeric((lambda.p.h*g.p/(1-(1-g.p)*s.p)-1)/(lambda.b.h*g.b/(1-(1-g.b)*s.b)-1)*sqrt((abp.h*abb.h)/(apb.h*app.h)))
k.bp.c <- as.numeric((lambda.p.c*g.p/(1-(1-g.p)*s.p)-1)/(lambda.b.c*g.b/(1-(1-g.b)*s.b)-1)*sqrt((abp.c*abb.c)/(apb.c*app.c)))
k.bc.h <- as.numeric((lambda.c.h*g.c/(1-(1-g.c)*s.c)-1)/(lambda.b.h*g.b/(1-(1-g.b)*s.b)-1)*sqrt((abc.h*abb.h)/(acb.h*acc.h)))
k.bc.c <- as.numeric((lambda.c.c*g.c/(1-(1-g.c)*s.c)-1)/(lambda.b.c*g.b/(1-(1-g.b)*s.b)-1)*sqrt((abc.c*abb.c)/(acb.c*acc.c)))
k.bn.h <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1)/(lambda.b.h*g.b/(1-(1-g.b)*s.b)-1)*sqrt((abn.h*abb.h)/(anb.h*ann.h)))
k.bn.c <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1)/(lambda.b.c*g.b/(1-(1-g.b)*s.b)-1)*sqrt((abn.c*abb.c)/(anb.c*ann.c)))
k.pc.h <- as.numeric((lambda.c.h*g.c/(1-(1-g.c)*s.c)-1)/(lambda.p.h*g.p/(1-(1-g.p)*s.p)-1)*sqrt((apc.h*app.h)/(acp.h*acc.h)))
k.pc.c <- as.numeric((lambda.c.c*g.c/(1-(1-g.c)*s.c)-1)/(lambda.p.c*g.p/(1-(1-g.p)*s.p)-1)*sqrt((apc.c*app.c)/(acp.c*acc.c)))
k.pn.h <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1)/(lambda.p.h*g.p/(1-(1-g.p)*s.p)-1)*sqrt((apn.h*app.h)/(anp.h*ann.h)))
k.pn.c <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1)/(lambda.p.c*g.p/(1-(1-g.p)*s.p)-1)*sqrt((apn.c*app.c)/(anp.c*ann.c)))
k.cn.h <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1)/(lambda.c.h*g.c/(1-(1-g.c)*s.c)-1)*sqrt((acn.h*acc.h)/(anc.h*ann.h)))
k.cn.c <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1)/(lambda.c.c*g.c/(1-(1-g.c)*s.c)-1)*sqrt((acn.c*acc.c)/(anc.c*ann.c)))

# Calculate fitness differences (ln[k2/k1] where species 2 is dominant competitor)
log(1/k.sb.h)
log(1/k.sb.c)
log(1/k.sp.h)
log(1/k.sp.c)
log(1/k.sc.h)
log(1/k.sc.c)
log(1/k.sn.h)
log(1/k.sn.c)
log(k.bp.h)
log(k.bp.c)
log(k.bc.h)
log(k.bc.c)
log(1/k.bn.h)
log(1/k.bn.c)
log(k.pc.h)
log(k.pc.c)
log(1/k.pn.h)
log(1/k.pn.c)
log(1/k.cn.h)
log(1/k.cn.c)


# STATISTICAL ANALYSES REPORTED IN MAIN TEXT
# Evaluate whether niche differences (niche overlap) are significantly different between hand and control treatments
rho.h <- c(rho.sb.h,rho.sp.h,rho.sc.h,rho.sn.h,rho.bp.h,rho.bc.h,rho.bn.h,rho.pc.h,rho.pn.h,rho.cn.h)
rho.c <- c(rho.sb.c,rho.sp.c,rho.sc.c,rho.sn.c,rho.bp.c,rho.bc.c,rho.bn.c,rho.pc.c,rho.pn.c,rho.cn.c)
ratio <- (1/rho.h)/(1/rho.c)
t.test(ratio, mu=1)

# Evaluate whether fitness differences are significantly different between hand and control treatments
k.h <- c(1/k.sb.h,1/k.sp.h,1/k.sc.h,1/k.sn.h,k.bp.h,k.bc.h,1/k.bn.h,k.pc.h,1/k.pn.h,1/k.cn.h)
k.c <- c(1/k.sb.c,1/k.sp.c,1/k.sc.c,1/k.sn.c,k.bp.c,k.bc.c,1/k.bn.c,k.pc.c,1/k.pn.c,1/k.cn.c)
ratio <- k.h/k.c
t.test(ratio, mu=1)

# Binomial test for probability of 0/10 pairs coexisting in control treatment given 3/10 pairs coexisting in hand pollination treatment
binom.test(0, 10, p=0.3, alternative = "two.sided")


# ANALYSES FOR SUPPLEMENTARY DISCUSSION
# Calculated sigma_iB B terms
sigmaSB <- lambda.s.h / lambda.s.c - 1
sigmaBB <- lambda.b.h / lambda.b.c - 1
sigmaPB <- lambda.p.h / lambda.p.c - 1
sigmaCB <- lambda.c.h / lambda.c.c - 1
sigmaNB <- lambda.n.h / lambda.n.c - 1

# Calculate fitness differences without competition with background competitors in other plots (Eq. S9)
k.sb.h.new <- as.numeric((lambda.b.h*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asb.h*ass.h)/(abs.h*abb.h)))
k.sb.c.new <- as.numeric((lambda.b.c*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asb.c*ass.c)/(abs.c*abb.c)))
k.sp.h.new <- as.numeric((lambda.p.h*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asp.h*ass.h)/(aps.h*app.h)))
k.sp.c.new <- as.numeric((lambda.p.c*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asp.c*ass.c)/(aps.c*app.c)))
k.sc.h.new <- as.numeric((lambda.c.h*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asc.h*ass.h)/(acs.h*acc.h)))
k.sc.c.new <- as.numeric((lambda.c.c*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asc.c*ass.c)/(acs.c*acc.c)))
k.sn.h.new <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.s.h*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asn.h*ass.h)/(ans.h*ann.h)))
k.sn.c.new <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.s.c*g.s/(1-(1-g.s)*s.s)-1-sigmaSB)*sqrt((asn.c*ass.c)/(ans.c*ann.c)))
k.bp.h.new <- as.numeric((lambda.p.h*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)/(lambda.b.h*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)*sqrt((abp.h*abb.h)/(apb.h*app.h)))
k.bp.c.new <- as.numeric((lambda.p.c*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)/(lambda.b.c*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)*sqrt((abp.c*abb.c)/(apb.c*app.c)))
k.bc.h.new <- as.numeric((lambda.c.h*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)/(lambda.b.h*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)*sqrt((abc.h*abb.h)/(acb.h*acc.h)))
k.bc.c.new <- as.numeric((lambda.c.c*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)/(lambda.b.c*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)*sqrt((abc.c*abb.c)/(acb.c*acc.c)))
k.bn.h.new <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.b.h*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)*sqrt((abn.h*abb.h)/(anb.h*ann.h)))
k.bn.c.new <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.b.c*g.b/(1-(1-g.b)*s.b)-1-sigmaBB)*sqrt((abn.c*abb.c)/(anb.c*ann.c)))
k.pc.h.new <- as.numeric((lambda.c.h*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)/(lambda.p.h*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)*sqrt((apc.h*app.h)/(acp.h*acc.h)))
k.pc.c.new <- as.numeric((lambda.c.c*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)/(lambda.p.c*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)*sqrt((apc.c*app.c)/(acp.c*acc.c)))
k.pn.h.new <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.p.h*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)*sqrt((apn.h*app.h)/(anp.h*ann.h)))
k.pn.c.new <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.p.c*g.p/(1-(1-g.p)*s.p)-1-sigmaPB)*sqrt((apn.c*app.c)/(anp.c*ann.c)))
k.cn.h.new <- as.numeric((lambda.n.h*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.c.h*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)*sqrt((acn.h*acc.h)/(anc.h*ann.h)))
k.cn.c.new <- as.numeric((lambda.n.c*g.n/(1-(1-g.n)*s.n)-1-sigmaNB)/(lambda.c.c*g.c/(1-(1-g.c)*s.c)-1-sigmaCB)*sqrt((acn.c*acc.c)/(anc.c*ann.c)))

# Calculate fitness differences (ln[k2/k1] where species 2 is dominant competitor)
log(1/k.sb.h.new)
log(1/k.sb.c.new)
log(1/k.sp.h.new)
log(1/k.sp.c.new)
log(1/k.sc.h.new)
log(1/k.sc.c.new)
log(1/k.sn.h.new)
log(1/k.sn.c.new)
log(k.bp.h.new)
log(k.bp.c.new)
log(k.bc.h.new)
log(k.bc.c.new)
log(1/k.bn.h.new)
log(1/k.bn.c.new)
log(1/k.pc.h.new)
log(1/k.pc.c.new)
log(1/k.pn.h.new)
log(1/k.pn.c.new)
log(1/k.cn.h.new)
log(1/k.cn.c.new)



###############################################################################################################################################
# EFFECTS OF POLLINATOR DECLINE (Figure 4)
# data (Sinapis, Buglossoides, Papaver, Centaurea, Nigella)
data2 <- data.frame(vis = c(-0.89,-1.47,-2.25,-0.34,1.26),
                    pgr = c(-0.33,-0.25,-0.70,-0.17,0.20))
summary(lm(pgr ~ vis, data = data2))



###############################################################################################################################################
# ERROR PROPAGATION FOR NICHE AND FITNESS DIFFERENCES (Figure 2)
# parameter estimate and standard errors (=sd) as is required for 'propagate' package
# Hand pollination treatments
# Sinapis
p.lambda.s.h<-as.vector(s.par[2,c(1,2)]*g.s)/(1-(1-g.s)*s.s)
p.ass.h<-as.vector(s.par[4,c(1,2)])
p.asb.h<-as.vector(s.par[6,c(1,2)])
p.asp.h<-as.vector(s.par[8,c(1,2)])
p.asc.h<-as.vector(s.par[10,c(1,2)])
p.asn.h<-as.vector(s.par[12,c(1,2)])
# Buglossoides
p.lambda.b.h<-as.vector(b.par[2,c(1,2)]*g.b)/(1-(1-g.b)*s.b)
p.abs.h<-as.vector(b.par[4,c(1,2)])
p.abb.h<-as.vector(b.par[6,c(1,2)])
p.abp.h<-as.vector(b.par[8,c(1,2)])
p.abc.h<-as.vector(b.par[10,c(1,2)])
p.abn.h<-as.vector(b.par[12,c(1,2)])
# Papaver
p.lambda.p.h<-as.vector(p.par[2,c(1,2)]*g.p)/(1-(1-g.p)*s.p)
p.aps.h<-as.vector(p.par[4,c(1,2)])
p.apb.h<-as.vector(p.par[6,c(1,2)])
p.app.h<-as.vector(p.par[8,c(1,2)])
p.apc.h<-as.vector(p.par[10,c(1,2)])
p.apn.h<-as.vector(p.par[12,c(1,2)])
# Centraurea
p.lambda.c.h<-as.vector(c.par[2,c(1,2)]*g.c)/(1-(1-g.c)*s.c)
p.acs.h<-as.vector(c.par[4,c(1,2)])
p.acb.h<-as.vector(c.par[6,c(1,2)])
p.acp.h<-as.vector(c.par[8,c(1,2)])
p.acc.h<-as.vector(c.par[10,c(1,2)])
p.acn.h<-as.vector(c.par[12,c(1,2)])
# Nigella
p.lambda.n.h<-as.vector(n.par[2,c(1,2)]*g.n)/(1-(1-g.n)*s.n)
p.ans.h<-as.vector(n.par[4,c(1,2)])
p.anb.h<-as.vector(n.par[6,c(1,2)])
p.anp.h<-as.vector(n.par[8,c(1,2)])
p.anc.h<-as.vector(n.par[10,c(1,2)])
p.ann.h<-as.vector(n.par[12,c(1,2)])

# Control Treatments
# Sinapis
p.lambda.s.c<-as.vector(s.par[1,c(1,2)]*g.s)/(1-(1-g.s)*s.s)
p.ass.c<-as.vector(s.par[3,c(1,2)])
p.asb.c<-as.vector(s.par[5,c(1,2)])
p.asp.c<-as.vector(s.par[7,c(1,2)])
p.asc.c<-as.vector(s.par[9,c(1,2)])
p.asn.c<-as.vector(s.par[11,c(1,2)])
# Buglossoides
p.lambda.b.c<-as.vector(b.par[1,c(1,2)]*g.b)/(1-(1-g.b)*s.b)
p.abs.c<-as.vector(b.par[3,c(1,2)])
p.abb.c<-as.vector(b.par[5,c(1,2)])
p.abp.c<-as.vector(b.par[7,c(1,2)])
p.abc.c<-as.vector(b.par[9,c(1,2)])
p.abn.c<-as.vector(b.par[11,c(1,2)])
# Papaver
p.lambda.p.c<-as.vector(p.par[1,c(1,2)]*g.p)/(1-(1-g.p)*s.p)
p.aps.c<-as.vector(p.par[3,c(1,2)])
p.apb.c<-as.vector(p.par[5,c(1,2)])
p.app.c<-as.vector(p.par[7,c(1,2)])
p.apc.c<-as.vector(p.par[9,c(1,2)])
p.apn.c<-as.vector(p.par[11,c(1,2)])
# Centaurea
p.lambda.c.c<-as.vector(c.par[1,c(1,2)]*g.c)/(1-(1-g.c)*s.c)
p.acs.c<-as.vector(c.par[3,c(1,2)])
p.acb.c<-as.vector(c.par[5,c(1,2)])
p.acp.c<-as.vector(c.par[7,c(1,2)])
p.acc.c<-as.vector(c.par[9,c(1,2)])
p.acn.c<-as.vector(c.par[11,c(1,2)])
# Nigella
p.lambda.n.c<-as.vector(n.par[1,c(1,2)]*g.n)/(1-(1-g.n)*s.n)
p.ans.c<-as.vector(n.par[3,c(1,2)])
p.anb.c<-as.vector(n.par[5,c(1,2)])
p.anp.c<-as.vector(n.par[7,c(1,2)])
p.anc.c<-as.vector(n.par[9,c(1,2)])
p.ann.c<-as.vector(n.par[11,c(1,2)])


# NICHE DIFFERENCE IN HAND POLLINATION TREATMENT
# make dataframe with parameter estimates and standard errors
data.niche.sb.h <- cbind(p.ass.h,p.asb.h,p.abb.h,p.abs.h)
data.niche.sp.h <- cbind(p.ass.h,p.asp.h,p.app.h,p.aps.h)
data.niche.sc.h <- cbind(p.ass.h,p.asc.h,p.acc.h,p.acs.h)
data.niche.sn.h <- cbind(p.ass.h,p.asn.h,p.ann.h,p.ans.h)
data.niche.bp.h <- cbind(p.abb.h,p.abp.h,p.app.h,p.apb.h)
data.niche.bc.h <- cbind(p.abb.h,p.abc.h,p.acc.h,p.acb.h)
data.niche.bn.h <- cbind(p.abb.h,p.abn.h,p.ann.h,p.anb.h)
data.niche.pc.h <- cbind(p.app.h,p.apc.h,p.acc.h,p.acp.h)
data.niche.pn.h <- cbind(p.app.h,p.apn.h,p.ann.h,p.anp.h)
data.niche.cn.h <- cbind(p.acc.h,p.acn.h,p.ann.h,p.anc.h)

# covariance matrices
# Sinapis
niche.sb.h.vcov<-matrix(c(vcov(S.full)[4,4],vcov(S.full)[6,4],vcov(S.full)[4,6],vcov(S.full)[6,6],0,0,0,0),nrow=2,byrow=F)
niche.sp.h.vcov<-matrix(c(vcov(S.full)[4,4],vcov(S.full)[8,4],vcov(S.full)[4,8],vcov(S.full)[8,8],0,0,0,0),nrow=2,byrow=F)
niche.sc.h.vcov<-matrix(c(vcov(S.full)[4,4],vcov(S.full)[10,4],vcov(S.full)[4,10],vcov(S.full)[10,10],0,0,0,0),nrow=2,byrow=F)
niche.sn.h.vcov<-matrix(c(vcov(S.full)[4,4],vcov(S.full)[12,4],vcov(S.full)[4,12],vcov(S.full)[12,12],0,0,0,0),nrow=2,byrow=F)
# Buglossoides
niche.bs.h.vcov<-matrix(c(0,0,0,0,vcov(B.full)[6,6],vcov(B.full)[4,6],vcov(B.full)[6,4],vcov(B.full)[4,4]),nrow=2,byrow=F)
niche.bp.h.vcov<-matrix(c(vcov(B.full)[6,6],vcov(B.full)[8,6],vcov(B.full)[6,8],vcov(B.full)[8,8],0,0,0,0),nrow=2,byrow=F)
niche.bc.h.vcov<-matrix(c(vcov(B.full)[6,6],vcov(B.full)[10,6],vcov(B.full)[6,10],vcov(B.full)[10,10],0,0,0,0),nrow=2,byrow=F)
niche.bn.h.vcov<-matrix(c(vcov(B.full)[6,6],vcov(B.full)[12,6],vcov(B.full)[6,12],vcov(B.full)[12,12],0,0,0,0),nrow=2,byrow=F)
# Papaver
niche.ps.h.vcov<-matrix(c(0,0,0,0,vcov(P.full)[8,8],vcov(P.full)[4,8],vcov(P.full)[8,4],vcov(P.full)[4,4]),nrow=2,byrow=F)
niche.pb.h.vcov<-matrix(c(0,0,0,0,vcov(P.full)[8,8],vcov(P.full)[6,8],vcov(P.full)[8,6],vcov(P.full)[6,6]),nrow=2,byrow=F)
niche.pc.h.vcov<-matrix(c(vcov(P.full)[8,8],vcov(P.full)[10,8],vcov(P.full)[8,10],vcov(P.full)[10,10],0,0,0,0),nrow=2,byrow=F)
niche.pn.h.vcov<-matrix(c(vcov(P.full)[8,8],vcov(P.full)[12,8],vcov(P.full)[8,12],vcov(P.full)[12,12],0,0,0,0),nrow=2,byrow=F)
# Centaurea
niche.cs.h.vcov<-matrix(c(0,0,0,0,vcov(C.full)[10,10],vcov(C.full)[4,10],vcov(C.full)[10,4],vcov(C.full)[4,4]),nrow=2,byrow=F)
niche.cb.h.vcov<-matrix(c(0,0,0,0,vcov(C.full)[10,10],vcov(C.full)[6,10],vcov(C.full)[10,6],vcov(C.full)[6,6]),nrow=2,byrow=F)
niche.cp.h.vcov<-matrix(c(0,0,0,0,vcov(C.full)[10,10],vcov(C.full)[8,10],vcov(C.full)[10,8],vcov(C.full)[8,8]),nrow=2,byrow=F)
niche.cn.h.vcov<-matrix(c(vcov(C.full)[10,10],vcov(C.full)[12,10],vcov(C.full)[10,12],vcov(C.full)[12,12],0,0,0,0),nrow=2,byrow=F)
# Nigella
niche.ns.h.vcov<-matrix(c(0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[4,12],vcov(N.full)[12,4],vcov(N.full)[4,4]),nrow=2,byrow=F)
niche.nb.h.vcov<-matrix(c(0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[6,12],vcov(N.full)[12,6],vcov(N.full)[6,6]),nrow=2,byrow=F)
niche.np.h.vcov<-matrix(c(0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[8,12],vcov(N.full)[12,8],vcov(N.full)[8,8]),nrow=2,byrow=F)
niche.nc.h.vcov<-matrix(c(0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[10,12],vcov(N.full)[12,10],vcov(N.full)[10,10]),nrow=2,byrow=F)

# combine the covariance matrices, assuming zero covariance between parameters estimated separately
niche.SiBu.h.vcov<-rbind(niche.sb.h.vcov,niche.bs.h.vcov)
colnames(niche.SiBu.h.vcov)<-c("p.ass.h","p.asb.h","p.abb.h","p.abs.h")
rownames(niche.SiBu.h.vcov)<-c("p.ass.h","p.asb.h","p.abb.h","p.abs.h")
niche.SiPa.h.vcov<-rbind(niche.sp.h.vcov,niche.ps.h.vcov)
colnames(niche.SiPa.h.vcov)<-c("p.ass.h","p.asp.h","p.app.h","p.aps.h")
rownames(niche.SiPa.h.vcov)<-c("p.ass.h","p.asp.h","p.app.h","p.aps.h")
niche.SiCe.h.vcov<-rbind(niche.sc.h.vcov,niche.cs.h.vcov)
colnames(niche.SiCe.h.vcov)<-c("p.ass.h","p.asc.h","p.acc.h","p.acs.h")
rownames(niche.SiCe.h.vcov)<-c("p.ass.h","p.asc.h","p.acc.h","p.acs.h")
niche.SiNi.h.vcov<-rbind(niche.sn.h.vcov,niche.ns.h.vcov)
colnames(niche.SiNi.h.vcov)<-c("p.ass.h","p.asn.h","p.ann.h","p.ans.h")
rownames(niche.SiNi.h.vcov)<-c("p.ass.h","p.asn.h","p.ann.h","p.ans.h")
niche.BuPa.h.vcov<-rbind(niche.bp.h.vcov,niche.pb.h.vcov)
colnames(niche.BuPa.h.vcov)<-c("p.abb.h","p.abp.h","p.app.h","p.apb.h")
rownames(niche.BuPa.h.vcov)<-c("p.abb.h","p.abp.h","p.app.h","p.apb.h")
niche.BuCe.h.vcov<-rbind(niche.bc.h.vcov,niche.cb.h.vcov)
colnames(niche.BuCe.h.vcov)<-c("p.abb.h","p.abc.h","p.acc.h","p.acb.h")
rownames(niche.BuCe.h.vcov)<-c("p.abb.h","p.abc.h","p.acc.h","p.acb.h")
niche.BuNi.h.vcov<-rbind(niche.bn.h.vcov,niche.nb.h.vcov)
colnames(niche.BuNi.h.vcov)<-c("p.abb.h","p.abn.h","p.ann.h","p.anb.h")
rownames(niche.BuNi.h.vcov)<-c("p.abb.h","p.abn.h","p.ann.h","p.anb.h")
niche.PaCe.h.vcov<-rbind(niche.pc.h.vcov,niche.cp.h.vcov)
colnames(niche.PaCe.h.vcov)<-c("p.app.h","p.apc.h","p.acc.h","p.acp.h")
rownames(niche.PaCe.h.vcov)<-c("p.app.h","p.apc.h","p.acc.h","p.acp.h")
niche.PaNi.h.vcov<-rbind(niche.pn.h.vcov,niche.np.h.vcov)
colnames(niche.PaNi.h.vcov)<-c("p.app.h","p.apn.h","p.ann.h","p.anp.h")
rownames(niche.PaNi.h.vcov)<-c("p.app.h","p.apn.h","p.ann.h","p.anp.h")
niche.CeNi.h.vcov<-rbind(niche.cn.h.vcov,niche.nc.h.vcov)
colnames(niche.CeNi.h.vcov)<-c("p.acc.h","p.acn.h","p.ann.h","p.anc.h")
rownames(niche.CeNi.h.vcov)<-c("p.acc.h","p.acn.h","p.ann.h","p.anc.h")

# an expression for niche overlap (rho)
niche.expr.sb.h<-expression(sqrt((p.asb.h*p.abs.h)/(p.ass.h*p.abb.h)))
niche.expr.sp.h<-expression(sqrt((p.asp.h*p.aps.h)/(p.ass.h*p.app.h)))
niche.expr.sc.h<-expression(sqrt((p.asc.h*p.acs.h)/(p.ass.h*p.acc.h)))
niche.expr.sn.h<-expression(sqrt((p.asn.h*p.ans.h)/(p.ass.h*p.ann.h)))
niche.expr.bp.h<-expression(sqrt((p.abp.h*p.apb.h)/(p.abb.h*p.app.h)))
niche.expr.bc.h<-expression(sqrt((p.abc.h*p.acb.h)/(p.abb.h*p.acc.h)))
niche.expr.bn.h<-expression(sqrt((p.abn.h*p.anb.h)/(p.abb.h*p.ann.h)))
niche.expr.pc.h<-expression(sqrt((p.apc.h*p.acp.h)/(p.app.h*p.acc.h)))
niche.expr.pn.h<-expression(sqrt((p.apn.h*p.anp.h)/(p.app.h*p.ann.h)))
niche.expr.cn.h<-expression(sqrt((p.acn.h*p.anc.h)/(p.acc.h*p.ann.h)))

# propagation of uncertainty, including the covariance matrix
niche.sb.h.error <- propagate(niche.expr.sb.h,data=data.niche.sb.h,cov=niche.SiBu.h.vcov)
niche.sp.h.error <- propagate(niche.expr.sp.h,data=data.niche.sp.h,cov=niche.SiPa.h.vcov)
niche.sc.h.error <- propagate(niche.expr.sc.h,data=data.niche.sc.h,cov=niche.SiCe.h.vcov)
niche.sn.h.error <- propagate(niche.expr.sn.h,data=data.niche.sn.h,cov=niche.SiNi.h.vcov)
niche.bp.h.error <- propagate(niche.expr.bp.h,data=data.niche.bp.h,cov=niche.BuPa.h.vcov)
niche.bc.h.error <- propagate(niche.expr.bc.h,data=data.niche.bc.h,cov=niche.BuCe.h.vcov)
niche.bn.h.error <- propagate(niche.expr.bn.h,data=data.niche.bn.h,cov=niche.BuNi.h.vcov)
niche.pc.h.error <- propagate(niche.expr.pc.h,data=data.niche.pc.h,cov=niche.PaCe.h.vcov)
niche.pn.h.error <- propagate(niche.expr.pn.h,data=data.niche.pn.h,cov=niche.PaNi.h.vcov)
niche.cn.h.error <- propagate(niche.expr.cn.h,data=data.niche.cn.h,cov=niche.CeNi.h.vcov)

# mean values (-ln[rho])
-log(as.numeric(niche.sb.h.error$`prop`[1]))
-log(as.numeric(niche.sp.h.error$`prop`[1]))
-log(as.numeric(niche.sc.h.error$`prop`[1]))
-log(as.numeric(niche.sn.h.error$`prop`[1]))
-log(as.numeric(niche.bp.h.error$`prop`[1]))
-log(as.numeric(niche.bc.h.error$`prop`[1]))
-log(as.numeric(niche.bn.h.error$`prop`[1]))
-log(as.numeric(niche.pc.h.error$`prop`[1]))
-log(as.numeric(niche.pn.h.error$`prop`[1]))
-log(as.numeric(niche.cn.h.error$`prop`[1]))
# standard deviation (natural log transformed)
as.numeric(niche.sb.h.error$`prop`[3])/as.numeric(niche.sb.h.error$`prop`[1])
as.numeric(niche.sp.h.error$`prop`[3])/as.numeric(niche.sp.h.error$`prop`[1])
as.numeric(niche.sc.h.error$`prop`[3])/as.numeric(niche.sc.h.error$`prop`[1])
as.numeric(niche.sn.h.error$`prop`[3])/as.numeric(niche.sn.h.error$`prop`[1])
as.numeric(niche.bp.h.error$`prop`[3])/as.numeric(niche.bp.h.error$`prop`[1])
as.numeric(niche.bc.h.error$`prop`[3])/as.numeric(niche.bc.h.error$`prop`[1])
as.numeric(niche.bn.h.error$`prop`[3])/as.numeric(niche.bn.h.error$`prop`[1])
as.numeric(niche.pc.h.error$`prop`[3])/as.numeric(niche.pc.h.error$`prop`[1])
as.numeric(niche.pn.h.error$`prop`[3])/as.numeric(niche.pn.h.error$`prop`[1])
as.numeric(niche.cn.h.error$`prop`[3])/as.numeric(niche.cn.h.error$`prop`[1])


# NICHE DIFFERENCES IN CONTROL AMBIENT POLLINATION TREATMENT
# make dataframe with parameter estimates and standard errors
data.niche.sb.c <- cbind(p.ass.c,p.asb.c,p.abb.c,p.abs.c)
data.niche.sp.c <- cbind(p.ass.c,p.asp.c,p.app.c,p.aps.c)
data.niche.sc.c <- cbind(p.ass.c,p.asc.c,p.acc.c,p.acs.c)
data.niche.sn.c <- cbind(p.ass.c,p.asn.c,p.ann.c,p.ans.c)
data.niche.bp.c <- cbind(p.abb.c,p.abp.c,p.app.c,p.apb.c)
data.niche.bc.c <- cbind(p.abb.c,p.abc.c,p.acc.c,p.acb.c)
data.niche.bn.c <- cbind(p.abb.c,p.abn.c,p.ann.c,p.anb.c)
data.niche.pc.c <- cbind(p.app.c,p.apc.c,p.acc.c,p.acp.c)
data.niche.pn.c <- cbind(p.app.c,p.apn.c,p.ann.c,p.anp.c)
data.niche.cn.c <- cbind(p.acc.c,p.acn.c,p.ann.c,p.anc.c)

# covariance matrices
# Sinapis
niche.sb.c.vcov<-matrix(c(vcov(S.full)[3,3],vcov(S.full)[5,3],vcov(S.full)[3,5],vcov(S.full)[5,5],0,0,0,0),nrow=2,byrow=F)
niche.sp.c.vcov<-matrix(c(vcov(S.full)[3,3],vcov(S.full)[7,3],vcov(S.full)[3,7],vcov(S.full)[7,7],0,0,0,0),nrow=2,byrow=F)
niche.sc.c.vcov<-matrix(c(vcov(S.full)[3,3],vcov(S.full)[9,3],vcov(S.full)[3,9],vcov(S.full)[9,9],0,0,0,0),nrow=2,byrow=F)
niche.sn.c.vcov<-matrix(c(vcov(S.full)[3,3],vcov(S.full)[11,3],vcov(S.full)[3,11],vcov(S.full)[11,11],0,0,0,0),nrow=2,byrow=F)
# Buglossoides
niche.bs.c.vcov<-matrix(c(0,0,0,0,vcov(B.full)[5,5],vcov(B.full)[3,5],vcov(B.full)[5,3],vcov(B.full)[3,3]),nrow=2,byrow=F)
niche.bp.c.vcov<-matrix(c(vcov(B.full)[5,5],vcov(B.full)[7,5],vcov(B.full)[5,7],vcov(B.full)[7,7],0,0,0,0),nrow=2,byrow=F)
niche.bc.c.vcov<-matrix(c(vcov(B.full)[5,5],vcov(B.full)[9,5],vcov(B.full)[5,9],vcov(B.full)[9,9],0,0,0,0),nrow=2,byrow=F)
niche.bn.c.vcov<-matrix(c(vcov(B.full)[5,5],vcov(B.full)[11,5],vcov(B.full)[5,11],vcov(B.full)[11,11],0,0,0,0),nrow=2,byrow=F)
# Papaver
niche.ps.c.vcov<-matrix(c(0,0,0,0,vcov(P.full)[7,7],vcov(P.full)[3,7],vcov(P.full)[7,3],vcov(P.full)[3,3]),nrow=2,byrow=F)
niche.pb.c.vcov<-matrix(c(0,0,0,0,vcov(P.full)[7,7],vcov(P.full)[5,7],vcov(P.full)[7,5],vcov(P.full)[5,5]),nrow=2,byrow=F)
niche.pc.c.vcov<-matrix(c(vcov(P.full)[7,7],vcov(P.full)[9,7],vcov(P.full)[7,9],vcov(P.full)[9,9],0,0,0,0),nrow=2,byrow=F)
niche.pn.c.vcov<-matrix(c(vcov(P.full)[7,7],vcov(P.full)[11,7],vcov(P.full)[7,11],vcov(P.full)[11,11],0,0,0,0),nrow=2,byrow=F)
# Centaurea
niche.cs.c.vcov<-matrix(c(0,0,0,0,vcov(C.full)[9,9],vcov(C.full)[3,9],vcov(C.full)[9,3],vcov(C.full)[3,3]),nrow=2,byrow=F)
niche.cb.c.vcov<-matrix(c(0,0,0,0,vcov(C.full)[9,9],vcov(C.full)[5,9],vcov(C.full)[9,5],vcov(C.full)[5,5]),nrow=2,byrow=F)
niche.cp.c.vcov<-matrix(c(0,0,0,0,vcov(C.full)[9,9],vcov(C.full)[7,9],vcov(C.full)[9,7],vcov(C.full)[7,7]),nrow=2,byrow=F)
niche.cn.c.vcov<-matrix(c(vcov(C.full)[9,9],vcov(C.full)[11,9],vcov(C.full)[9,11],vcov(C.full)[11,11],0,0,0,0),nrow=2,byrow=F)
# Nigella
niche.ns.c.vcov<-matrix(c(0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[3,11],vcov(N.full)[11,3],vcov(N.full)[3,3]),nrow=2,byrow=F)
niche.nb.c.vcov<-matrix(c(0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[5,11],vcov(N.full)[11,5],vcov(N.full)[5,5]),nrow=2,byrow=F)
niche.np.c.vcov<-matrix(c(0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[7,11],vcov(N.full)[11,7],vcov(N.full)[7,7]),nrow=2,byrow=F)
niche.nc.c.vcov<-matrix(c(0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[9,11],vcov(N.full)[11,9],vcov(N.full)[9,9]),nrow=2,byrow=F)

# combine the covariance matrices, assuming zero covariance between parameters estimated separately
niche.SiBu.c.vcov<-rbind(niche.sb.c.vcov,niche.bs.c.vcov)
colnames(niche.SiBu.c.vcov)<-c("p.ass.c","p.asb.c","p.abb.c","p.abs.c")
rownames(niche.SiBu.c.vcov)<-c("p.ass.c","p.asb.c","p.abb.c","p.abs.c")
niche.SiPa.c.vcov<-rbind(niche.sp.c.vcov,niche.ps.c.vcov)
colnames(niche.SiPa.c.vcov)<-c("p.ass.c","p.asp.c","p.app.c","p.aps.c")
rownames(niche.SiPa.c.vcov)<-c("p.ass.c","p.asp.c","p.app.c","p.aps.c")
niche.SiCe.c.vcov<-rbind(niche.sc.c.vcov,niche.cs.c.vcov)
colnames(niche.SiCe.c.vcov)<-c("p.ass.c","p.asc.c","p.acc.c","p.acs.c")
rownames(niche.SiCe.c.vcov)<-c("p.ass.c","p.asc.c","p.acc.c","p.acs.c")
niche.SiNi.c.vcov<-rbind(niche.sn.c.vcov,niche.ns.c.vcov)
colnames(niche.SiNi.c.vcov)<-c("p.ass.c","p.asn.c","p.ann.c","p.ans.c")
rownames(niche.SiNi.c.vcov)<-c("p.ass.c","p.asn.c","p.ann.c","p.ans.c")
niche.BuPa.c.vcov<-rbind(niche.bp.c.vcov,niche.pb.c.vcov)
colnames(niche.BuPa.c.vcov)<-c("p.abb.c","p.abp.c","p.app.c","p.apb.c")
rownames(niche.BuPa.c.vcov)<-c("p.abb.c","p.abp.c","p.app.c","p.apb.c")
niche.BuCe.c.vcov<-rbind(niche.bc.c.vcov,niche.cb.c.vcov)
colnames(niche.BuCe.c.vcov)<-c("p.abb.c","p.abc.c","p.acc.c","p.acb.c")
rownames(niche.BuCe.c.vcov)<-c("p.abb.c","p.abc.c","p.acc.c","p.acb.c")
niche.BuNi.c.vcov<-rbind(niche.bn.c.vcov,niche.nb.c.vcov)
colnames(niche.BuNi.c.vcov)<-c("p.abb.c","p.abn.c","p.ann.c","p.anb.c")
rownames(niche.BuNi.c.vcov)<-c("p.abb.c","p.abn.c","p.ann.c","p.anb.c")
niche.PaCe.c.vcov<-rbind(niche.pc.c.vcov,niche.cp.c.vcov)
colnames(niche.PaCe.c.vcov)<-c("p.app.c","p.apc.c","p.acc.c","p.acp.c")
rownames(niche.PaCe.c.vcov)<-c("p.app.c","p.apc.c","p.acc.c","p.acp.c")
niche.PaNi.c.vcov<-rbind(niche.pn.c.vcov,niche.np.c.vcov)
colnames(niche.PaNi.c.vcov)<-c("p.app.c","p.apn.c","p.ann.c","p.anp.c")
rownames(niche.PaNi.c.vcov)<-c("p.app.c","p.apn.c","p.ann.c","p.anp.c")
niche.CeNi.c.vcov<-rbind(niche.cn.c.vcov,niche.nc.c.vcov)
colnames(niche.CeNi.c.vcov)<-c("p.acc.c","p.acn.c","p.ann.c","p.anc.c")
rownames(niche.CeNi.c.vcov)<-c("p.acc.c","p.acn.c","p.ann.c","p.anc.c")

# an expression for niche overlap (rho)
niche.expr.sb.c<-expression(sqrt((p.asb.c*p.abs.c)/(p.ass.c*p.abb.c)))
niche.expr.sp.c<-expression(sqrt((p.asp.c*p.aps.c)/(p.ass.c*p.app.c)))
niche.expr.sc.c<-expression(sqrt((p.asc.c*p.acs.c)/(p.ass.c*p.acc.c)))
niche.expr.sn.c<-expression(sqrt((p.asn.c*p.ans.c)/(p.ass.c*p.ann.c)))
niche.expr.bp.c<-expression(sqrt((p.abp.c*p.apb.c)/(p.abb.c*p.app.c)))
niche.expr.bc.c<-expression(sqrt((p.abc.c*p.acb.c)/(p.abb.c*p.acc.c)))
niche.expr.bn.c<-expression(sqrt((p.abn.c*p.anb.c)/(p.abb.c*p.ann.c)))
niche.expr.pc.c<-expression(sqrt((p.apc.c*p.acp.c)/(p.app.c*p.acc.c)))
niche.expr.pn.c<-expression(sqrt((p.apn.c*p.anp.c)/(p.app.c*p.ann.c)))
niche.expr.cn.c<-expression(sqrt((p.acn.c*p.anc.c)/(p.acc.c*p.ann.c)))

# propagation of uncertainty, including the covariance matrix
niche.sb.c.error<-propagate(niche.expr.sb.c,data=data.niche.sb.c,cov=niche.SiBu.c.vcov)
niche.sp.c.error<-propagate(niche.expr.sp.c,data=data.niche.sp.c,cov=niche.SiPa.c.vcov)
niche.sc.c.error<-propagate(niche.expr.sc.c,data=data.niche.sc.c,cov=niche.SiCe.c.vcov)
niche.sn.c.error<-propagate(niche.expr.sn.c,data=data.niche.sn.c,cov=niche.SiNi.c.vcov)
niche.bp.c.error<-propagate(niche.expr.bp.c,data=data.niche.bp.c,cov=niche.BuPa.c.vcov)
niche.bc.c.error<-propagate(niche.expr.bc.c,data=data.niche.bc.c,cov=niche.BuCe.c.vcov)
niche.bn.c.error<-propagate(niche.expr.bn.c,data=data.niche.bn.c,cov=niche.BuNi.c.vcov)
niche.pc.c.error<-propagate(niche.expr.pc.c,data=data.niche.pc.c,cov=niche.PaCe.c.vcov)
niche.pn.c.error<-propagate(niche.expr.pn.c,data=data.niche.pn.c,cov=niche.PaNi.c.vcov)
niche.cn.c.error<-propagate(niche.expr.cn.c,data=data.niche.cn.c,cov=niche.CeNi.c.vcov)

# mean values (-ln[rho])
-log(as.numeric(niche.sb.c.error$`prop`[1]))
-log(as.numeric(niche.sp.c.error$`prop`[1]))
-log(as.numeric(niche.sc.c.error$`prop`[1]))
-log(as.numeric(niche.sn.c.error$`prop`[1]))
-log(as.numeric(niche.bp.c.error$`prop`[1]))
-log(as.numeric(niche.bc.c.error$`prop`[1]))
-log(as.numeric(niche.bn.c.error$`prop`[1]))
-log(as.numeric(niche.pc.c.error$`prop`[1]))
-log(as.numeric(niche.pn.c.error$`prop`[1]))
-log(as.numeric(niche.cn.c.error$`prop`[1]))
# standard deviation (natural log transformed)
as.numeric(niche.sb.c.error$`prop`[3])/as.numeric(niche.sb.c.error$`prop`[1])
as.numeric(niche.sp.c.error$`prop`[3])/as.numeric(niche.sp.c.error$`prop`[1])
as.numeric(niche.sc.c.error$`prop`[3])/as.numeric(niche.sc.c.error$`prop`[1])
as.numeric(niche.sn.c.error$`prop`[3])/as.numeric(niche.sn.c.error$`prop`[1])
as.numeric(niche.bp.c.error$`prop`[3])/as.numeric(niche.bp.c.error$`prop`[1])
as.numeric(niche.bc.c.error$`prop`[3])/as.numeric(niche.bc.c.error$`prop`[1])
as.numeric(niche.bn.c.error$`prop`[3])/as.numeric(niche.bn.c.error$`prop`[1])
as.numeric(niche.pc.c.error$`prop`[3])/as.numeric(niche.pc.c.error$`prop`[1])
as.numeric(niche.pn.c.error$`prop`[3])/as.numeric(niche.pn.c.error$`prop`[1])
as.numeric(niche.cn.c.error$`prop`[3])/as.numeric(niche.cn.c.error$`prop`[1])


# FITNESS DIFFERENCES IN HAND POLLINATION TREATMENTS
# make dataframe with parameter estimates and standard errors
data.fitness.sb.h <- cbind(p.lambda.s.h,p.ass.h,p.asb.h,p.abb.h,p.abs.h,p.lambda.b.h)
data.fitness.sp.h <- cbind(p.lambda.s.h,p.ass.h,p.asp.h,p.app.h,p.aps.h,p.lambda.p.h)
data.fitness.sc.h <- cbind(p.lambda.s.h,p.ass.h,p.asc.h,p.acc.h,p.acs.h,p.lambda.c.h)
data.fitness.sn.h <- cbind(p.lambda.s.h,p.ass.h,p.asn.h,p.ann.h,p.ans.h,p.lambda.n.h)
data.fitness.bp.h <- cbind(p.lambda.b.h,p.abb.h,p.abp.h,p.app.h,p.apb.h,p.lambda.p.h)
data.fitness.bc.h <- cbind(p.lambda.b.h,p.abb.h,p.abc.h,p.acc.h,p.acb.h,p.lambda.c.h)
data.fitness.bn.h <- cbind(p.lambda.b.h,p.abb.h,p.abn.h,p.ann.h,p.anb.h,p.lambda.n.h)
data.fitness.pc.h <- cbind(p.lambda.p.h,p.app.h,p.apc.h,p.acc.h,p.acp.h,p.lambda.c.h)
data.fitness.pn.h <- cbind(p.lambda.p.h,p.app.h,p.apn.h,p.ann.h,p.anp.h,p.lambda.n.h)
data.fitness.cn.h <- cbind(p.lambda.c.h,p.acc.h,p.acn.h,p.ann.h,p.anc.h,p.lambda.n.h)

# covariance matrices
# Sinapis
fitness.sb.h.vcov<-matrix(c(vcov(S.full)[2,2]*g.s^2,vcov(S.full)[4,2]*g.s,vcov(S.full)[6,2]*g.s,vcov(S.full)[2,4]*g.s,vcov(S.full)[4,4],vcov(S.full)[6,4],vcov(S.full)[2,6]*g.s,vcov(S.full)[4,6],vcov(S.full)[6,6],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.sp.h.vcov<-matrix(c(vcov(S.full)[2,2]*g.s^2,vcov(S.full)[4,2]*g.s,vcov(S.full)[8,2]*g.s,vcov(S.full)[2,4]*g.s,vcov(S.full)[4,4],vcov(S.full)[8,4],vcov(S.full)[2,8]*g.s,vcov(S.full)[4,8],vcov(S.full)[8,8],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.sc.h.vcov<-matrix(c(vcov(S.full)[2,2]*g.s^2,vcov(S.full)[4,2]*g.s,vcov(S.full)[10,2]*g.s,vcov(S.full)[2,4]*g.s,vcov(S.full)[4,4],vcov(S.full)[10,4],vcov(S.full)[2,10]*g.s,vcov(S.full)[4,10],vcov(S.full)[10,10],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.sn.h.vcov<-matrix(c(vcov(S.full)[2,2]*g.s^2,vcov(S.full)[4,2]*g.s,vcov(S.full)[12,2]*g.s,vcov(S.full)[2,4]*g.s,vcov(S.full)[4,4],vcov(S.full)[12,4],vcov(S.full)[2,12]*g.s,vcov(S.full)[4,12],vcov(S.full)[12,12],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Buglossoides
fitness.bs.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(B.full)[6,6],vcov(B.full)[4,6],vcov(B.full)[2,6]*g.b,vcov(B.full)[6,4],vcov(B.full)[4,4],vcov(B.full)[2,4]*g.b,vcov(B.full)[6,2]*g.b,vcov(B.full)[4,2]*g.b,vcov(B.full)[2,2]*g.b^2),nrow=3,byrow=F)
fitness.bp.h.vcov<-matrix(c(vcov(B.full)[2,2]*g.b^2,vcov(B.full)[6,2]*g.b,vcov(B.full)[8,2]*g.b,vcov(B.full)[2,6]*g.b,vcov(B.full)[6,6],vcov(B.full)[8,6],vcov(B.full)[2,8]*g.b,vcov(B.full)[6,8],vcov(B.full)[8,8],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.bc.h.vcov<-matrix(c(vcov(B.full)[2,2]*g.b^2,vcov(B.full)[6,2]*g.b,vcov(B.full)[10,2]*g.b,vcov(B.full)[2,6]*g.b,vcov(B.full)[6,6],vcov(B.full)[10,6],vcov(B.full)[2,10]*g.b,vcov(B.full)[6,10],vcov(B.full)[10,10],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.bn.h.vcov<-matrix(c(vcov(B.full)[2,2]*g.b^2,vcov(B.full)[6,2]*g.b,vcov(B.full)[12,2]*g.b,vcov(B.full)[2,6]*g.b,vcov(B.full)[6,6],vcov(B.full)[12,6],vcov(B.full)[2,12]*g.b,vcov(B.full)[6,12],vcov(B.full)[12,12],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Papaver
fitness.ps.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(P.full)[8,8],vcov(P.full)[4,8],vcov(P.full)[2,8]*g.p,vcov(P.full)[8,4],vcov(P.full)[4,4],vcov(P.full)[2,4]*g.p,vcov(P.full)[8,2]*g.p,vcov(P.full)[4,2]*g.p,vcov(P.full)[2,2]*g.p^2),nrow=3,byrow=F)
fitness.pb.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(P.full)[8,8],vcov(P.full)[6,8],vcov(P.full)[2,8]*g.p,vcov(P.full)[8,6],vcov(P.full)[6,6],vcov(P.full)[2,6]*g.p,vcov(P.full)[8,2]*g.p,vcov(P.full)[6,2]*g.p,vcov(P.full)[2,2]*g.p^2),nrow=3,byrow=F)
fitness.pc.h.vcov<-matrix(c(vcov(P.full)[2,2]*g.p^2,vcov(P.full)[8,2]*g.p,vcov(P.full)[10,2]*g.p,vcov(P.full)[2,8]*g.p,vcov(P.full)[8,8],vcov(P.full)[10,8],vcov(P.full)[2,10]*g.p,vcov(P.full)[8,10],vcov(P.full)[10,10],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.pn.h.vcov<-matrix(c(vcov(P.full)[2,2]*g.p^2,vcov(P.full)[8,2]*g.p,vcov(P.full)[12,2]*g.p,vcov(P.full)[2,8]*g.p,vcov(P.full)[8,8],vcov(P.full)[12,8],vcov(P.full)[2,12]*g.p,vcov(P.full)[8,12],vcov(P.full)[12,12],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Centaurea
fitness.cs.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(C.full)[10,10],vcov(C.full)[4,10],vcov(C.full)[2,10]*g.c,vcov(C.full)[10,4],vcov(C.full)[4,4],vcov(C.full)[2,4]*g.c,vcov(C.full)[10,2]*g.c,vcov(C.full)[4,2]*g.c,vcov(C.full)[2,2]*g.c^2),nrow=3,byrow=F)
fitness.cb.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(C.full)[10,10],vcov(C.full)[6,10],vcov(C.full)[2,10]*g.c,vcov(C.full)[10,6],vcov(C.full)[6,6],vcov(C.full)[2,6]*g.c,vcov(C.full)[10,2]*g.c,vcov(C.full)[6,2]*g.c,vcov(C.full)[2,2]*g.c^2),nrow=3,byrow=F)
fitness.cp.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(C.full)[10,10],vcov(C.full)[8,10],vcov(C.full)[2,10]*g.c,vcov(C.full)[10,8],vcov(C.full)[8,8],vcov(C.full)[2,8]*g.c,vcov(C.full)[10,2]*g.c,vcov(C.full)[8,2]*g.c,vcov(C.full)[2,2]*g.c^2),nrow=3,byrow=F)
fitness.cn.h.vcov<-matrix(c(vcov(C.full)[2,2]*g.c^2,vcov(C.full)[10,2]*g.c,vcov(C.full)[12,2]*g.c,vcov(C.full)[2,10]*g.c,vcov(C.full)[10,10],vcov(C.full)[12,10],vcov(C.full)[2,12]*g.c,vcov(C.full)[10,12],vcov(C.full)[12,12],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Nigella
fitness.ns.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[4,12],vcov(N.full)[2,12]*g.n,vcov(N.full)[12,4],vcov(N.full)[4,4],vcov(N.full)[2,4]*g.n,vcov(N.full)[12,2]*g.n,vcov(N.full)[4,2]*g.n,vcov(N.full)[2,2]*g.n^2),nrow=3,byrow=F)
fitness.nb.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[6,12],vcov(N.full)[2,12]*g.n,vcov(N.full)[12,6],vcov(N.full)[6,6],vcov(N.full)[2,6]*g.n,vcov(N.full)[12,2]*g.n,vcov(N.full)[6,2]*g.n,vcov(N.full)[2,2]*g.n^2),nrow=3,byrow=F)
fitness.np.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[8,12],vcov(N.full)[2,12]*g.n,vcov(N.full)[12,8],vcov(N.full)[8,8],vcov(N.full)[2,8]*g.n,vcov(N.full)[12,2]*g.n,vcov(N.full)[8,2]*g.n,vcov(N.full)[2,2]*g.n^2),nrow=3,byrow=F)
fitness.nc.h.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[12,12],vcov(N.full)[10,12],vcov(N.full)[2,12]*g.n,vcov(N.full)[12,10],vcov(N.full)[10,10],vcov(N.full)[2,10]*g.n,vcov(N.full)[12,2]*g.n,vcov(N.full)[10,2]*g.n,vcov(N.full)[2,2]*g.n^2),nrow=3,byrow=F)

# combine the covariance matrices, assuming zero covariance between parameters estimated separately
fitness.SiBu.h.vcov<-rbind(fitness.sb.h.vcov,fitness.bs.h.vcov)
colnames(fitness.SiBu.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asb.h","p.abb.h","p.abs.h","p.lambda.b.h")
rownames(fitness.SiBu.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asb.h","p.abb.h","p.abs.h","p.lambda.b.h")
fitness.SiPa.h.vcov<-rbind(fitness.sp.h.vcov,fitness.ps.h.vcov)
colnames(fitness.SiPa.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asp.h","p.app.h","p.aps.h","p.lambda.p.h")
rownames(fitness.SiPa.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asp.h","p.app.h","p.aps.h","p.lambda.p.h")
fitness.SiCe.h.vcov<-rbind(fitness.sc.h.vcov,fitness.cs.h.vcov)
colnames(fitness.SiCe.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asc.h","p.acc.h","p.acs.h","p.lambda.c.h")
rownames(fitness.SiCe.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asc.h","p.acc.h","p.acs.h","p.lambda.c.h")
fitness.SiNi.h.vcov<-rbind(fitness.sn.h.vcov,fitness.ns.h.vcov)
colnames(fitness.SiNi.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asn.h","p.ann.h","p.ans.h","p.lambda.n.h")
rownames(fitness.SiNi.h.vcov)<-c("p.lambda.s.h","p.ass.h","p.asn.h","p.ann.h","p.ans.h","p.lambda.n.h")
fitness.BuPa.h.vcov<-rbind(fitness.bp.h.vcov,fitness.pb.h.vcov)
colnames(fitness.BuPa.h.vcov)<-c("p.lambda.b.h","p.abb.h","p.abp.h","p.app.h","p.apb.h","p.lambda.p.h")
rownames(fitness.BuPa.h.vcov)<-c("p.lambda.b.h","p.abb.h","p.abp.h","p.app.h","p.apb.h","p.lambda.p.h")
fitness.BuCe.h.vcov<-rbind(fitness.bc.h.vcov,fitness.cb.h.vcov)
colnames(fitness.BuCe.h.vcov)<-c("p.lambda.b.h","p.abb.h","p.abc.h","p.acc.h","p.acb.h","p.lambda.c.h")
rownames(fitness.BuCe.h.vcov)<-c("p.lambda.b.h","p.abb.h","p.abc.h","p.acc.h","p.acb.h","p.lambda.c.h")
fitness.BuNi.h.vcov<-rbind(fitness.bn.h.vcov,fitness.nb.h.vcov)
colnames(fitness.BuNi.h.vcov)<-c("p.lambda.b.h","p.abb.h","p.abn.h","p.ann.h","p.anb.h","p.lambda.n.h")
rownames(fitness.BuNi.h.vcov)<-c("p.lambda.b.h","p.abb.h","p.abn.h","p.ann.h","p.anb.h","p.lambda.n.h")
fitness.PaCe.h.vcov<-rbind(fitness.pc.h.vcov,fitness.cp.h.vcov)
colnames(fitness.PaCe.h.vcov)<-c("p.lambda.p.h","p.app.h","p.apc.h","p.acc.h","p.acp.h","p.lambda.c.h")
rownames(fitness.PaCe.h.vcov)<-c("p.lambda.p.h","p.app.h","p.apc.h","p.acc.h","p.acp.h","p.lambda.c.h")
fitness.PaNi.h.vcov<-rbind(fitness.pn.h.vcov,fitness.np.h.vcov)
colnames(fitness.PaNi.h.vcov)<-c("p.lambda.p.h","p.app.h","p.apn.h","p.ann.h","p.anp.h","p.lambda.n.h")
rownames(fitness.PaNi.h.vcov)<-c("p.lambda.p.h","p.app.h","p.apn.h","p.ann.h","p.anp.h","p.lambda.n.h")
fitness.CeNi.h.vcov<-rbind(fitness.cn.h.vcov,fitness.nc.h.vcov)
colnames(fitness.CeNi.h.vcov)<-c("p.lambda.c.h","p.acc.h","p.acn.h","p.ann.h","p.anc.h","p.lambda.n.h")
rownames(fitness.CeNi.h.vcov)<-c("p.lambda.c.h","p.acc.h","p.acn.h","p.ann.h","p.anc.h","p.lambda.n.h")

# an expression for fitness differences (k_j/k_i) (Note: g, s incorporated within lambda)
fitness.expr.sb.h<-expression((p.lambda.b.h-1)/(p.lambda.s.h-1)*sqrt((p.asb.h*p.ass.h)/(p.abs.h*p.abb.h)))
fitness.expr.sp.h<-expression((p.lambda.p.h-1)/(p.lambda.s.h-1)*sqrt((p.asp.h*p.ass.h)/(p.aps.h*p.app.h)))
fitness.expr.sc.h<-expression((p.lambda.c.h-1)/(p.lambda.s.h-1)*sqrt((p.asc.h*p.ass.h)/(p.acs.h*p.acc.h)))
fitness.expr.sn.h<-expression((p.lambda.n.h-1)/(p.lambda.s.h-1)*sqrt((p.asn.h*p.ass.h)/(p.ans.h*p.ann.h)))
fitness.expr.bp.h<-expression((p.lambda.p.h-1)/(p.lambda.b.h-1)*sqrt((p.abp.h*p.abb.h)/(p.apb.h*p.app.h)))
fitness.expr.bc.h<-expression((p.lambda.c.h-1)/(p.lambda.b.h-1)*sqrt((p.abc.h*p.abb.h)/(p.acb.h*p.acc.h)))
fitness.expr.bn.h<-expression((p.lambda.n.h-1)/(p.lambda.b.h-1)*sqrt((p.abn.h*p.abb.h)/(p.anb.h*p.ann.h)))
fitness.expr.pc.h<-expression((p.lambda.c.h-1)/(p.lambda.p.h-1)*sqrt((p.apc.h*p.app.h)/(p.acp.h*p.acc.h)))
fitness.expr.pn.h<-expression((p.lambda.n.h-1)/(p.lambda.p.h-1)*sqrt((p.apn.h*p.app.h)/(p.anp.h*p.ann.h)))
fitness.expr.cn.h<-expression((p.lambda.n.h-1)/(p.lambda.c.h-1)*sqrt((p.acn.h*p.acc.h)/(p.anc.h*p.ann.h)))

# propagation of uncertainty, including the covariance matrix
fitness.sb.h.error<-propagate(fitness.expr.sb.h,data=data.fitness.sb.h,cov=fitness.SiBu.h.vcov)
fitness.sp.h.error<-propagate(fitness.expr.sp.h,data=data.fitness.sp.h,cov=fitness.SiPa.h.vcov)
fitness.sc.h.error<-propagate(fitness.expr.sc.h,data=data.fitness.sc.h,cov=fitness.SiCe.h.vcov)
fitness.sn.h.error<-propagate(fitness.expr.sn.h,data=data.fitness.sn.h,cov=fitness.SiNi.h.vcov)
fitness.bp.h.error<-propagate(fitness.expr.bp.h,data=data.fitness.bp.h,cov=fitness.BuPa.h.vcov)
fitness.bc.h.error<-propagate(fitness.expr.bc.h,data=data.fitness.bc.h,cov=fitness.BuCe.h.vcov)
fitness.bn.h.error<-propagate(fitness.expr.bn.h,data=data.fitness.bn.h,cov=fitness.BuNi.h.vcov)
fitness.pc.h.error<-propagate(fitness.expr.pc.h,data=data.fitness.pc.h,cov=fitness.PaCe.h.vcov)
fitness.pn.h.error<-propagate(fitness.expr.pn.h,data=data.fitness.pn.h,cov=fitness.PaNi.h.vcov)
fitness.cn.h.error<-propagate(fitness.expr.cn.h,data=data.fitness.cn.h,cov=fitness.CeNi.h.vcov)

# mean values (ln[k2/k1])
log(1/as.numeric(fitness.sb.h.error$`prop`[1]))
log(1/as.numeric(fitness.sp.h.error$`prop`[1]))
log(1/as.numeric(fitness.sc.h.error$`prop`[1]))
log(1/as.numeric(fitness.sn.h.error$`prop`[1]))
log(as.numeric(fitness.bp.h.error$`prop`[1]))
log(as.numeric(fitness.bc.h.error$`prop`[1]))
log(1/as.numeric(fitness.bn.h.error$`prop`[1]))
log(as.numeric(fitness.pc.h.error$`prop`[1]))
log(1/as.numeric(fitness.pn.h.error$`prop`[1]))
log(1/as.numeric(fitness.cn.h.error$`prop`[1]))
# standard deviation (natural log transformed)
as.numeric(fitness.sb.h.error$`prop`[3])/as.numeric(fitness.sb.h.error$`prop`[1])
as.numeric(fitness.sp.h.error$`prop`[3])/as.numeric(fitness.sp.h.error$`prop`[1])
as.numeric(fitness.sc.h.error$`prop`[3])/as.numeric(fitness.sc.h.error$`prop`[1])
as.numeric(fitness.sn.h.error$`prop`[3])/as.numeric(fitness.sn.h.error$`prop`[1])
as.numeric(fitness.bp.h.error$`prop`[3])/as.numeric(fitness.bp.h.error$`prop`[1])
as.numeric(fitness.bc.h.error$`prop`[3])/as.numeric(fitness.bc.h.error$`prop`[1])
as.numeric(fitness.bn.h.error$`prop`[3])/as.numeric(fitness.bn.h.error$`prop`[1])
as.numeric(fitness.pc.h.error$`prop`[3])/as.numeric(fitness.pc.h.error$`prop`[1])
as.numeric(fitness.pn.h.error$`prop`[3])/as.numeric(fitness.pn.h.error$`prop`[1])
as.numeric(fitness.cn.h.error$`prop`[3])/as.numeric(fitness.cn.h.error$`prop`[1])


# FITNESS DIFFERENCE IN CONTROL AMBIENT POLLINATION TREATMENTS
# make dataframe with parameter estimates and standard errors
data.fitness.sb.c <- cbind(p.lambda.s.c,p.ass.c,p.asb.c,p.abb.c,p.abs.c,p.lambda.b.c)
data.fitness.sp.c <- cbind(p.lambda.s.c,p.ass.c,p.asp.c,p.app.c,p.aps.c,p.lambda.p.c)
data.fitness.sc.c <- cbind(p.lambda.s.c,p.ass.c,p.asc.c,p.acc.c,p.acs.c,p.lambda.c.c)
data.fitness.sn.c <- cbind(p.lambda.s.c,p.ass.c,p.asn.c,p.ann.c,p.ans.c,p.lambda.n.c)
data.fitness.bp.c <- cbind(p.lambda.b.c,p.abb.c,p.abp.c,p.app.c,p.apb.c,p.lambda.p.c)
data.fitness.bc.c <- cbind(p.lambda.b.c,p.abb.c,p.abc.c,p.acc.c,p.acb.c,p.lambda.c.c)
data.fitness.bn.c <- cbind(p.lambda.b.c,p.abb.c,p.abn.c,p.ann.c,p.anb.c,p.lambda.n.c)
data.fitness.pc.c <- cbind(p.lambda.p.c,p.app.c,p.apc.c,p.acc.c,p.acp.c,p.lambda.c.c)
data.fitness.pn.c <- cbind(p.lambda.p.c,p.app.c,p.apn.c,p.ann.c,p.anp.c,p.lambda.n.c)
data.fitness.cn.c <- cbind(p.lambda.c.c,p.acc.c,p.acn.c,p.ann.c,p.anc.c,p.lambda.n.c)

# covariance matrices
# Sinapis
fitness.sb.c.vcov<-matrix(c(vcov(S.full)[1,1]*g.s^2,vcov(S.full)[3,1]*g.s,vcov(S.full)[5,1]*g.s,vcov(S.full)[1,3]*g.s,vcov(S.full)[3,3],vcov(S.full)[5,3],vcov(S.full)[1,5]*g.s,vcov(S.full)[3,5],vcov(S.full)[5,5],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.sp.c.vcov<-matrix(c(vcov(S.full)[1,1]*g.s^2,vcov(S.full)[3,1]*g.s,vcov(S.full)[7,1]*g.s,vcov(S.full)[1,3]*g.s,vcov(S.full)[3,3],vcov(S.full)[7,3],vcov(S.full)[1,7]*g.s,vcov(S.full)[3,7],vcov(S.full)[7,7],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.sc.c.vcov<-matrix(c(vcov(S.full)[1,1]*g.s^2,vcov(S.full)[3,1]*g.s,vcov(S.full)[9,1]*g.s,vcov(S.full)[1,3]*g.s,vcov(S.full)[3,3],vcov(S.full)[9,3],vcov(S.full)[1,9]*g.s,vcov(S.full)[3,9],vcov(S.full)[9,9],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.sn.c.vcov<-matrix(c(vcov(S.full)[1,1]*g.s^2,vcov(S.full)[3,1]*g.s,vcov(S.full)[11,1]*g.s,vcov(S.full)[1,3]*g.s,vcov(S.full)[3,3],vcov(S.full)[11,3],vcov(S.full)[1,11]*g.s,vcov(S.full)[3,11],vcov(S.full)[11,11],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Buglossoides
fitness.bs.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(B.full)[5,5],vcov(B.full)[3,5],vcov(B.full)[1,5]*g.b,vcov(B.full)[5,3],vcov(B.full)[3,3],vcov(B.full)[1,3]*g.b,vcov(B.full)[5,1]*g.b,vcov(B.full)[3,1]*g.b,vcov(B.full)[1,1]*g.b^2),nrow=3,byrow=F)
fitness.bp.c.vcov<-matrix(c(vcov(B.full)[1,1]*g.b^2,vcov(B.full)[5,1]*g.b,vcov(B.full)[7,1]*g.b,vcov(B.full)[1,5]*g.b,vcov(B.full)[5,5],vcov(B.full)[7,5],vcov(B.full)[1,7]*g.b,vcov(B.full)[5,7],vcov(B.full)[7,7],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.bc.c.vcov<-matrix(c(vcov(B.full)[1,1]*g.b^2,vcov(B.full)[5,1]*g.b,vcov(B.full)[9,1]*g.b,vcov(B.full)[1,5]*g.b,vcov(B.full)[5,5],vcov(B.full)[9,5],vcov(B.full)[1,9]*g.b,vcov(B.full)[5,9],vcov(B.full)[9,9],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.bn.c.vcov<-matrix(c(vcov(B.full)[1,1]*g.b^2,vcov(B.full)[5,1]*g.b,vcov(B.full)[11,1]*g.b,vcov(B.full)[1,5]*g.b,vcov(B.full)[5,5],vcov(B.full)[11,5],vcov(B.full)[1,11]*g.b,vcov(B.full)[5,11],vcov(B.full)[11,11],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Papaver
fitness.ps.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(P.full)[7,7],vcov(P.full)[3,7],vcov(P.full)[1,7]*g.p,vcov(P.full)[7,3],vcov(P.full)[3,3],vcov(P.full)[1,3]*g.p,vcov(P.full)[7,1]*g.p,vcov(P.full)[3,1]*g.p,vcov(P.full)[1,1]*g.p^2),nrow=3,byrow=F)
fitness.pb.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(P.full)[7,7],vcov(P.full)[5,7],vcov(P.full)[1,7]*g.p,vcov(P.full)[7,5],vcov(P.full)[5,5],vcov(P.full)[1,5]*g.p,vcov(P.full)[7,1]*g.p,vcov(P.full)[5,1]*g.p,vcov(P.full)[1,1]*g.p^2),nrow=3,byrow=F)
fitness.pc.c.vcov<-matrix(c(vcov(P.full)[1,1]*g.p^2,vcov(P.full)[7,1]*g.p,vcov(P.full)[9,1]*g.p,vcov(P.full)[1,7]*g.p,vcov(P.full)[7,7],vcov(P.full)[9,7],vcov(P.full)[1,9]*g.p,vcov(P.full)[7,9],vcov(P.full)[9,9],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
fitness.pn.c.vcov<-matrix(c(vcov(P.full)[1,1]*g.p^2,vcov(P.full)[7,1]*g.p,vcov(P.full)[11,1]*g.p,vcov(P.full)[1,7]*g.p,vcov(P.full)[7,7],vcov(P.full)[11,7],vcov(P.full)[1,11]*g.p,vcov(P.full)[7,11],vcov(P.full)[11,11],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Centaurea
fitness.cs.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(C.full)[9,9],vcov(C.full)[3,9],vcov(C.full)[1,9]*g.c,vcov(C.full)[9,3],vcov(C.full)[3,3],vcov(C.full)[1,3]*g.c,vcov(C.full)[9,1]*g.c,vcov(C.full)[3,1]*g.c,vcov(C.full)[1,1]*g.c^2),nrow=3,byrow=F)
fitness.cb.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(C.full)[9,9],vcov(C.full)[5,9],vcov(C.full)[1,9]*g.c,vcov(C.full)[9,5],vcov(C.full)[5,5],vcov(C.full)[1,5]*g.c,vcov(C.full)[9,1]*g.c,vcov(C.full)[5,1]*g.c,vcov(C.full)[1,1]*g.c^2),nrow=3,byrow=F)
fitness.cp.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(C.full)[9,9],vcov(C.full)[7,9],vcov(C.full)[1,9]*g.c,vcov(C.full)[9,7],vcov(C.full)[7,7],vcov(C.full)[1,7]*g.c,vcov(C.full)[9,1]*g.c,vcov(C.full)[7,1]*g.c,vcov(C.full)[1,1]*g.c^2),nrow=3,byrow=F)
fitness.cn.c.vcov<-matrix(c(vcov(C.full)[1,1]*g.c^2,vcov(C.full)[9,1]*g.c,vcov(C.full)[11,1]*g.c,vcov(C.full)[1,9]*g.c,vcov(C.full)[9,9],vcov(C.full)[11,9],vcov(C.full)[1,11]*g.c,vcov(C.full)[9,11],vcov(C.full)[11,11],0,0,0,0,0,0,0,0,0),nrow=3,byrow=F)
# Nigella
fitness.ns.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[3,11],vcov(N.full)[1,11]*g.n,vcov(N.full)[11,3],vcov(N.full)[3,3],vcov(N.full)[1,3]*g.n,vcov(N.full)[11,1]*g.n,vcov(N.full)[3,1]*g.n,vcov(N.full)[1,1]*g.n^2),nrow=3,byrow=F)
fitness.nb.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[5,11],vcov(N.full)[1,11]*g.n,vcov(N.full)[11,5],vcov(N.full)[5,5],vcov(N.full)[1,5]*g.n,vcov(N.full)[11,1]*g.n,vcov(N.full)[5,1]*g.n,vcov(N.full)[1,1]*g.n^2),nrow=3,byrow=F)
fitness.np.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[7,11],vcov(N.full)[1,11]*g.n,vcov(N.full)[11,7],vcov(N.full)[7,7],vcov(N.full)[1,7]*g.n,vcov(N.full)[11,1]*g.n,vcov(N.full)[7,1]*g.n,vcov(N.full)[1,1]*g.n^2),nrow=3,byrow=F)
fitness.nc.c.vcov<-matrix(c(0,0,0,0,0,0,0,0,0,vcov(N.full)[11,11],vcov(N.full)[9,11],vcov(N.full)[1,11]*g.n,vcov(N.full)[11,9],vcov(N.full)[9,9],vcov(N.full)[1,9]*g.n,vcov(N.full)[11,1]*g.n,vcov(N.full)[9,1]*g.n,vcov(N.full)[1,1]*g.n^2),nrow=3,byrow=F)

# combine the covariance matrices, assuming zero covariance between parameters estimated separately
fitness.SiBu.c.vcov<-rbind(fitness.sb.c.vcov,fitness.bs.c.vcov)
colnames(fitness.SiBu.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asb.c","p.abb.c","p.abs.c","p.lambda.b.c")
rownames(fitness.SiBu.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asb.c","p.abb.c","p.abs.c","p.lambda.b.c")
fitness.SiPa.c.vcov<-rbind(fitness.sp.c.vcov,fitness.ps.c.vcov)
colnames(fitness.SiPa.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asp.c","p.app.c","p.aps.c","p.lambda.p.c")
rownames(fitness.SiPa.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asp.c","p.app.c","p.aps.c","p.lambda.p.c")
fitness.SiCe.c.vcov<-rbind(fitness.sc.c.vcov,fitness.cs.c.vcov)
colnames(fitness.SiCe.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asc.c","p.acc.c","p.acs.c","p.lambda.c.c")
rownames(fitness.SiCe.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asc.c","p.acc.c","p.acs.c","p.lambda.c.c")
fitness.SiNi.c.vcov<-rbind(fitness.sn.c.vcov,fitness.ns.c.vcov)
colnames(fitness.SiNi.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asn.c","p.ann.c","p.ans.c","p.lambda.n.c")
rownames(fitness.SiNi.c.vcov)<-c("p.lambda.s.c","p.ass.c","p.asn.c","p.ann.c","p.ans.c","p.lambda.n.c")
fitness.BuPa.c.vcov<-rbind(fitness.bp.c.vcov,fitness.pb.c.vcov)
colnames(fitness.BuPa.c.vcov)<-c("p.lambda.b.c","p.abb.c","p.abp.c","p.app.c","p.apb.c","p.lambda.p.c")
rownames(fitness.BuPa.c.vcov)<-c("p.lambda.b.c","p.abb.c","p.abp.c","p.app.c","p.apb.c","p.lambda.p.c")
fitness.BuCe.c.vcov<-rbind(fitness.bc.c.vcov,fitness.cb.c.vcov)
colnames(fitness.BuCe.c.vcov)<-c("p.lambda.b.c","p.abb.c","p.abc.c","p.acc.c","p.acb.c","p.lambda.c.c")
rownames(fitness.BuCe.c.vcov)<-c("p.lambda.b.c","p.abb.c","p.abc.c","p.acc.c","p.acb.c","p.lambda.c.c")
fitness.BuNi.c.vcov<-rbind(fitness.bn.c.vcov,fitness.nb.c.vcov)
colnames(fitness.BuNi.c.vcov)<-c("p.lambda.b.c","p.abb.c","p.abn.c","p.ann.c","p.anb.c","p.lambda.n.c")
rownames(fitness.BuNi.c.vcov)<-c("p.lambda.b.c","p.abb.c","p.abn.c","p.ann.c","p.anb.c","p.lambda.n.c")
fitness.PaCe.c.vcov<-rbind(fitness.pc.c.vcov,fitness.cp.c.vcov)
colnames(fitness.PaCe.c.vcov)<-c("p.lambda.p.c","p.app.c","p.apc.c","p.acc.c","p.acp.c","p.lambda.c.c")
rownames(fitness.PaCe.c.vcov)<-c("p.lambda.p.c","p.app.c","p.apc.c","p.acc.c","p.acp.c","p.lambda.c.c")
fitness.PaNi.c.vcov<-rbind(fitness.pn.c.vcov,fitness.np.c.vcov)
colnames(fitness.PaNi.c.vcov)<-c("p.lambda.p.c","p.app.c","p.apn.c","p.ann.c","p.anp.c","p.lambda.n.c")
rownames(fitness.PaNi.c.vcov)<-c("p.lambda.p.c","p.app.c","p.apn.c","p.ann.c","p.anp.c","p.lambda.n.c")
fitness.CeNi.c.vcov<-rbind(fitness.cn.c.vcov,fitness.nc.c.vcov)
colnames(fitness.CeNi.c.vcov)<-c("p.lambda.c.c","p.acc.c","p.acn.c","p.ann.c","p.anc.c","p.lambda.n.c")
rownames(fitness.CeNi.c.vcov)<-c("p.lambda.c.c","p.acc.c","p.acn.c","p.ann.c","p.anc.c","p.lambda.n.c")

# an expression for fitness differences (k_j/k_i) (Note: g and s incorporated within lambda)
fitness.expr.sb.c<-expression((p.lambda.b.c-1)/(p.lambda.s.c-1)*sqrt((p.asb.c*p.ass.c)/(p.abs.c*p.abb.c)))
fitness.expr.sp.c<-expression((p.lambda.p.c-1)/(p.lambda.s.c-1)*sqrt((p.asp.c*p.ass.c)/(p.aps.c*p.app.c)))
fitness.expr.sc.c<-expression((p.lambda.c.c-1)/(p.lambda.s.c-1)*sqrt((p.asc.c*p.ass.c)/(p.acs.c*p.acc.c)))
fitness.expr.sn.c<-expression((p.lambda.n.c-1)/(p.lambda.s.c-1)*sqrt((p.asn.c*p.ass.c)/(p.ans.c*p.ann.c)))
fitness.expr.bp.c<-expression((p.lambda.p.c-1)/(p.lambda.b.c-1)*sqrt((p.abp.c*p.abb.c)/(p.apb.c*p.app.c)))
fitness.expr.bc.c<-expression((p.lambda.c.c-1)/(p.lambda.b.c-1)*sqrt((p.abc.c*p.abb.c)/(p.acb.c*p.acc.c)))
fitness.expr.bn.c<-expression((p.lambda.n.c-1)/(p.lambda.b.c-1)*sqrt((p.abn.c*p.abb.c)/(p.anb.c*p.ann.c)))
fitness.expr.pc.c<-expression((p.lambda.c.c-1)/(p.lambda.p.c-1)*sqrt((p.apc.c*p.app.c)/(p.acp.c*p.acc.c)))
fitness.expr.pn.c<-expression((p.lambda.n.c-1)/(p.lambda.p.c-1)*sqrt((p.apn.c*p.app.c)/(p.anp.c*p.ann.c)))
fitness.expr.cn.c<-expression((p.lambda.n.c-1)/(p.lambda.c.c-1)*sqrt((p.acn.c*p.acc.c)/(p.anc.c*p.ann.c)))

# propagation of uncertainty, including the covariance matrix (add sigma terms here for Eq. S9 in Supplementary Discussion)
fitness.sb.c.error<-propagate(fitness.expr.sb.c,data=data.fitness.sb.c,cov=fitness.SiBu.c.vcov)
fitness.sp.c.error<-propagate(fitness.expr.sp.c,data=data.fitness.sp.c,cov=fitness.SiPa.c.vcov)
fitness.sc.c.error<-propagate(fitness.expr.sc.c,data=data.fitness.sc.c,cov=fitness.SiCe.c.vcov)
fitness.sn.c.error<-propagate(fitness.expr.sn.c,data=data.fitness.sn.c,cov=fitness.SiNi.c.vcov)
fitness.bp.c.error<-propagate(fitness.expr.bp.c,data=data.fitness.bp.c,cov=fitness.BuPa.c.vcov)
fitness.bc.c.error<-propagate(fitness.expr.bc.c,data=data.fitness.bc.c,cov=fitness.BuCe.c.vcov)
fitness.bn.c.error<-propagate(fitness.expr.bn.c,data=data.fitness.bn.c,cov=fitness.BuNi.c.vcov)
fitness.pc.c.error<-propagate(fitness.expr.pc.c,data=data.fitness.pc.c,cov=fitness.PaCe.c.vcov)
fitness.pn.c.error<-propagate(fitness.expr.pn.c,data=data.fitness.pn.c,cov=fitness.PaNi.c.vcov)
fitness.cn.c.error<-propagate(fitness.expr.cn.c,data=data.fitness.cn.c,cov=fitness.CeNi.c.vcov)

# mean values (ln[k2/k1])
log(1/as.numeric(fitness.sb.c.error$`prop`[1]))
log(1/as.numeric(fitness.sp.c.error$`prop`[1]))
log(1/as.numeric(fitness.sc.c.error$`prop`[1]))
log(1/as.numeric(fitness.sn.c.error$`prop`[1]))
log(as.numeric(fitness.bp.c.error$`prop`[1]))
log(as.numeric(fitness.bc.c.error$`prop`[1]))
log(1/as.numeric(fitness.bn.c.error$`prop`[1]))
log(as.numeric(fitness.pc.c.error$`prop`[1]))
log(1/as.numeric(fitness.pn.c.error$`prop`[1]))
log(1/as.numeric(fitness.cn.c.error$`prop`[1]))
# standard deviation (natural log transformed)
as.numeric(fitness.sb.c.error$`prop`[3])/as.numeric(fitness.sb.c.error$`prop`[1])
as.numeric(fitness.sp.c.error$`prop`[3])/as.numeric(fitness.sp.c.error$`prop`[1])
as.numeric(fitness.sc.c.error$`prop`[3])/as.numeric(fitness.sc.c.error$`prop`[1])
as.numeric(fitness.sn.c.error$`prop`[3])/as.numeric(fitness.sn.c.error$`prop`[1])
as.numeric(fitness.bp.c.error$`prop`[3])/as.numeric(fitness.bp.c.error$`prop`[1])
as.numeric(fitness.bc.c.error$`prop`[3])/as.numeric(fitness.bc.c.error$`prop`[1])
as.numeric(fitness.bn.c.error$`prop`[3])/as.numeric(fitness.bn.c.error$`prop`[1])
as.numeric(fitness.pc.c.error$`prop`[3])/as.numeric(fitness.pc.c.error$`prop`[1])
as.numeric(fitness.pn.c.error$`prop`[3])/as.numeric(fitness.pn.c.error$`prop`[1])
as.numeric(fitness.cn.c.error$`prop`[3])/as.numeric(fitness.cn.c.error$`prop`[1])

