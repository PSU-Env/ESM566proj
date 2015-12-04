#### Bern Romey, 01Dec15, ESM566 Term Project, Westslope Cutthroat Trout (WCT)

####--- GLOBAL OPTIONS ---####

setwd("D:/R/ESM566Proj")
data <- read.csv("wct.csv")
op <- par(mfrow=c(1,1), xpd=NA) # set the default plot grid to 1X1, &set default clip region based on logical value. 
## If TRUE, clipped to figure region, if False, clipped to plot region.
require(dplyr) # require is same as library, loads package
source("cor.matrix.r")

##### --- KNOW YOUR DATA --- ####
str(data) # look at the structure of the data
summary(data)

## Variables to use in model selection
library(dplyr)
ct <- na.omit(data[ ,c(6,9,11:14,19:27)])  # variables to use in logit regression
dt <- dplyr::filter(ct, ct$WID_m > 0, ct$DEP > 0, ct$TH > 0) # IWCT P/A=1/0 & no zero water
ct$IWCT <- factor(ct$IWCT, levels = 0:1, labels = c("No CT", "CT")) # turn IWCT interger into a factor
wct <- dplyr::filter(ct, ct$WID_m > 0, ct$DEP > 0, ct$TH > 0) # P/A = No CT/CT
rm(ct)
library(MVN) # load MVN package
uniNorm(wct[ ,2:15], type = "SW", desc = T) # summary statistics and Shipiro-Wilks normality test

### Cumulative distribution function
attach(wct)
windows(7,5)
par(mfrow=c(2,3))
plot(ecdf(WID_m))
plot(ecdf(DEP_m))
plot(ecdf(ACW_m))
plot(ecdf(FIN))
plot(ecdf(BLD))
plot(ecdf(LWDP))
par(op)

## Cor.matrix plots for predictor groups
cor.matrix(wct[,c(1:6)]) # Habitat
cor.matrix(wct[ ,c(1,7:11)]) # Substrate
cor.matrix(wct[ ,c(1,12:15)]) # Wood

## Boxplots of predictors included in final model
windows(7,7)
par(mfrow=c(2,3))
boxplot(IWCT,ACW_m, xlab="Presence", ylab="Active Channel Width (m)", col="darkgreen")
boxplot(IWCT,WID_m, xlab="Presence", ylab="Wet Width (m)", col="green")
boxplot(IWCT,DEP_m, xlab="Presence", ylab="Channel Depth (m)",col="lightgreen")
boxplot(IWCT,LWDP, xlab="Presence", ylab="Large Wood (count)",col="brown")
boxplot(IWCT,FIN, xlab="Presence", ylab="Fine Substrate (%)",col="tan")
boxplot(IWCT,BLD, xlab="Presence", ylab="Boulder Substrate (%)",col="beige")
par(op)
detach(wct)

summary(wct)


####--- LOGISTIC REGRESSION ---####
## Normality & Equal variance assumption not required; Predict presence of Westslope Cutthroat Trout;
## Response variable (y): IWCT

## Classification Tree Model (CART) - variable classification selection tree (for model interaction terms)
library(rpart) # Loads the Recursive Partitioning & Regression Trees package (rpart), used for variable selection & interaction selection
par(xpd = NA) # set default clip region based on logical value
ct.tree<-rpart(as.factor(IWCT)~.,data=wct, method="class")
plot(ct.tree, branch = .6, compress = T)
text(ct.tree, use.n = T)  # break points, bottom of tree is the number of occurances (No CT/CT).
plotcp(ct.tree) # What splits are significant?
summary(ct.tree) # model summary
print(ct.tree) # a short summary of the summary output

##---- Full Model ----
ctlogit <- glm(IWCT ~ .,data=dt, binomial) # Full model, no interaction
ctlogit # short summary
summary(ctlogit)

## Chi-squared test for differenc in deviance significance. 
anova(ctlogit,test="Chisq") # Used as a rough indicator to find good model, sig. predictors

## Variance Invlation Factor (VIF) 
library(car)
vif(ctlogit) # Test for multicolinearity: > 4-5 suggest a problem, > 10 highly likely

## Variable Selection - Stepwise
step(ctlogit)

##--- Reduced Model ---
library(rpart)
modR <- glm(formula= IWCT~LENG_m+WID_m+DEP_m+ACW_m+GRA+BLD+BDR+LWDP, binomial,data=dt)
modR # varaible selection tree (no interaction) & stepwise selection methods were used
summary(modR) 

pchisq(1166-1166.9, 8) #  the deviance of the 8 predictor modR is higher than the full model (not ideal)

## Variance Invlation Factor (VIF) & AIC 
library(car)
vif(modR)
AIC(ctlogit,modR,k=2) # smaller is better

####--------------------####
## HYPOTHESIS TEST
####--------------------####
mod0 <- glm(formula= IWCT~1, binomial,data=dt) # null model
mod0
anova(ctlogit, mod0, test="Chi") # Is the slope of the full model = 0?

anova(modR,test="Chisq")# Chi-squared test for differenc in deviance significance. 
## Used as a rough indicator to find/check model varaible coefficients

anova(ctlogit,modR, test="Chi") # Is the reduced model (modR) model significantly worse than full (ctlogit) model?
## no p-value if model has interaction term???

## MAXIMUM LOG LIKELIHOOD (ML). it can be thought of as a chi-square value - smallest possible deviance 
## between the observed and predicted values (kind of like finding the best fitting line) 
logLik(modR) # smaller is better
logLik(ctlogit)

## LIKELIHOOD RATIO TEST - G: Is the reduced model an adequate fit?
with(modR, null.deviance - deviance) # difference in deviance between null and modR
with(modR, df.null - df.residual) # df
with(modR, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE)) # p-value
## With a chi-square of 249, 9 degrees of freedom, & an associated p-value of < 0.001 tells us 
## that our model as a whole fits significantly better than a null model. 

## GOODNESS OF FIT: Pseudo R^2
1-(deviance(modR)/deviance(mod0))
1-(deviance(ctlogit)/deviance(mod0))

## ODDS RATIO COEFFICIENTS (ORC) & 95% CI
orc <-exp(cbind(IWCT = coef(modR), confint(modR))) 
orc
## The results of the ORC can be checked by simply taking the natural log of each coefficient and comparing to the
## summary of the model output.
##---- example ----
## For a one unit increase in LWDP (large wood), the odds of WCT presence (vs not present) increase 
## by a factor of 1.39.  Values < 1, reduce the odds, radio = 1, have no influance, and > 1 increase the odds of Y

####-----------------------------------------####
## WCT prediction model (probability of presence)
####-----------------------------------------####
## Replace each habitat variable in the x vector (l,w,d,acw,gr,b,br,lwp) to get predicted presence of WCT
## Where: l=length, w=weted width, d=water depth, acw=active channel width, gr=%gravel, co=%cobble, b=%boulder, br=%bedrock,& lwp=number of large wood pieces
## Make sure to run orc above first.

library(faraway)
x <- c(6,2,.75,2.5,23,60,2,4)
ilogit(orc[1,1]-orc[2,1]*(x[1])-orc[3,1]*(x[2])-orc[4,1]*(x[3])-orc[5,1]*(x[4])-
  orc[6,1]*(x[5])+orc[7,1]*(x[6])-orc[8,1]*(x[7])+orc[9,1]*(x[8]))
rm(x)


###############
## logit Plots
#### ------####
attach(dt)
mod <- glm(IWCT~WID_m, binomial)

par(mfrow=c(1,2))
xv <- seq(0,16,0.01)
yv <- predict(mod, list(WID_m=xv), type="resp")
plot(WID_m,IWCT)
lines(xv,yv,col="red")

####-----------------------------------------####
## Confusion table - Lect 17 CART & Random Forest
# -------------------------------------------####
## Slide 12:13,15
## Lecture 17CART
dim(dt)
s <- sample(2, nrow(dt), replace = TRUE, prob=c(0.8, 0.2))  #split sample in to 80 & 20 %.  Use 20% to validate
table(s)
library(rpart)
ct.t <- rpart(as.factor(IWCT)~LENG_m + WID_m + DEP_m + ACW_m + GRA + BLD + BDR + LWDP,method="class",data=dt[s==1,])
plot(ct.t)
text(ct.t, use.n=TRUE)
summary(ct.t)
tb.1 <- table(predicted=predict(ct.t,data=dt[s==1,], type="class"),
              observed=dt[s==1, "IWCT"])

tb.1 # Correctly classified (0-0,1-1), missclasified (0-1,1-0)
cc <-(tb.1[1,1]+tb.1[2,2])/(sum(tb.1))
cc # Model proportion correctly classified (correct classification rate)
mc <-(1-cc)*100  # Model misclassificaiton rate 
mc
abs <- (tb.1[1,1])/(tb.1[1,1]+tb.1[2,1])*100 # proportion absent correctly classified
abs
prs <- (tb.1[2,2])/(tb.1[1,2]+tb.1[2,2])*100 # proportion present correctly classified
prs # CT = 1, No CT = 0



####-----------------------------------####
## Durbin -Wattson test for autocorrelation
####-----------------------------------####
library(car)
dw <-glm(formula = IWCT~DEP_m+ACW_m+WID_m+BLD+FIN+LWDP, binomial,data=dt)
dwt(dw$residuals)
dwt(dw)
require(lmtest)
dwtest(dw)


####-------####
## Lag analysis
####-------####
tsy<-ts(ctlogit$residuals)
plot(tsy)
hist(tsy)
summary(tsy)
plot(ctlogit$residuals~ctlogit$y)
lag.plot(tsy)
plot(ctlogit)
summary(ctlogit)

detach(dt)
