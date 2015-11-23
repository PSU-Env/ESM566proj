#### Bern Romey, 23Oct15, ESM566 Term Project, Westslope Cutthroat Trout (WCT)

#### GLOBAL OPTIONS ####
setwd("D:/R/ESM566Proj")
data <- read.csv("wct.csv")
op <- par(mfrow=c(1,1), xpd=NA) # set the default plot grid to 1X1
par(xpd = NA) # set default clip region based on logical value. If TRUE, clipped to figure region, if False, clipped to plot region.
require(dplyr) # require is same as library, loads package
source("cor.matrix.r")

##### ~~oo KNOW YOUR DATA  oo~~ ####
str(data) # look at the structure of the data
summary(data)

## Variables to use in model selection
ct <- na.omit(data[ ,c(6,9,11:14,19:27)])  # variables to use in logit regression
dt <- dplyr::filter(ct, ct$WID_m > 0, ct$DEP > 0, ct$TH > 0) # IWCT P/A=1/0
ct$IWCT <- factor(ct$IWCT, levels = 0:1, labels = c("No CT", "CT")) # turn IWCT interger into a factor
library(dplyr)
wct <- dplyr::filter(ct, ct$WID_m > 0, ct$DEP > 0, ct$TH > 0) # remove habitat with zero water width and depth
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
plot(IWCT,ACW_m, xlab="Presence", ylab="Active Channel Width (m)", col="darkgreen")
plot(IWCT,WID_m, xlab="Presence", ylab="Wet Width (m)", col="green")
plot(IWCT,DEP_m, xlab="Presence", ylab="Channel Depth (m)",col="lightgreen")
plot(IWCT,LWDP, xlab="Presence", ylab="Large Wood (count)",col="brown")
plot(IWCT,FIN, xlab="Presence", ylab="Fine Substrate (%)",col="tan")
plot(IWCT,BLD, xlab="Presence", ylab="Boulder Substrate (%)",col="beige")
par(op)
detach(wct)

#### ~~000ooOO LOGISTIC REGRESSION OOoo000~~ ####
## Normality & Equal variance assumption not required; Predict presence of West Slope Cutthroat Trout;
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

ctlogit <- glm(IWCT ~ .,data=dt, binomial) # Full model, no interaction
ctlogit # short summary
summary(ctlogit)
anova(ctlogit,test="Chisq") # Chi-squared test for differenc in deviance significance. Used as a rough indicator to find good model, sig. predictors

## Variance Invlation Factor (VIF) - Test for multicolinearity: > 4-5 suggest a problem, > 10 highly likely
library(car)
vif(ctlogit)

## Variable Selection - Stepwise
step(ctlogit)

#### ~oo8oo~ REDUCED MODEL ~oo8oo~ ####
#### varaible selection tree & stepwise selection methods
# modR <- glm(formula= IWCT~LENG_m + WID_m + DEP_m + ACW_m + GRA + COB + BLD + BDR + LWDP, binomial,data=dt)
modR <- glm(formula= IWCT~WID_m + WID_m:FIN+DEP_m + ACW_m + BLD + LWDP, binomial,data=dt)
modR
summary(modR)  # Width, depth, and active channel width have a significant negative effect on WCT presence.(large habitat less likely to have CT)
## Also, habitat higher numbers of boulders and LWDP are more likely to have WCT.

pchisq(1162-1164, 2) #  the deviance of the 9 predictor modR is higher than the full model(not good)
pchisq(1162-1141, 2) # the 5 predictor model with interaction is a better model with less residual deviance

## Variance Invlation Factor (VIF) - Test for multicolinearity: > 4-5 suggest a problem, > 10 highly likely
library(car)
vif(modR)

AIC(ctlogit,modR,k=2)

## NULL HYPOTHESIS TEST, the SLOPE of the logit model = 0
mod0 <- glm(formula= IWCT~1, binomial,data=dt) # null model
mod0
anova(modR, mod0, test="Chi")

anova(modR,test="Chisq")# Chi-squared test for differenc in deviance significance. Used as a rough indicator to find good model

## no p-vale if model has interaction term???
anova(ctlogit,modR, test="Chi") # Is the reduced model (modR) model significantly worse than full (ctlogit) model?

## MAXIMUM LOG LIKELIHOOD (ML). it can be thought of as a chi-square value - smallest possible deviance 
## between the observed and predicted values (kind of like finding the best fitting line) 
logLik(modR) # smaller is better
logLik(ctlogit)

## LIKELIHOOD RATIO TEST - G: Is the reduced model an adequate fit?
with(modR, null.deviance - deviance) # difference in deviance between null and modR
with(modR, df.null - df.residual) # df
with(modR, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE)) # p-value
## With a chi-square of 272.0, 6 degrees of freedom, & an associated p-value of less than 0.001 tells us 
## that our model as a whole fits significantly better than a null model. 

## GOODNESS OF FIT: Pseudo R^2
1-(deviance(modR)/deviance(mod0))
1-(deviance(ctlogit)/deviance(mod0))

hist(residuals(modR)) # the residuals should be similar to a chi^2 distribution

## ODDS RATIO COEFFICIENTS & 95% CI
exp(cbind(IWCT = coef(modR), confint(modR)))
## For a one unit increase in LWDP (large wood), the odds of WCT presence (vs not present) increase 
## by a factor of 1.39.  Values < 1, reduce the odds, radio = 1, have no influance, and > 1 increase the odds of Y



## logit Plots
# -------------------------------------
attach(dt)
mod <- glm(IWCT~WID_m, binomial)

par(mfrow=c(1,2))
xv <- seq(0,16,0.01)
yv <- predict(mod, list(WID_m=xv), type="resp")
plot(WID_m,IWCT)
lines(xv,yv,col="red")

### Confusion table - Lect 17 CART & Random Forest
# -------------------------------------
## Slide 12:13,15
dim(dt)
s <- sample(2, nrow(dt), replace = TRUE, prob=c(0.8, 0.2))
table(s)
library(rpart)
ct.t <- rpart(as.factor(IWCT)~IWCT~LENG_m + WID_m + DEP_m + ACW_m + GRA + COB + BLD + BDR + LWDP,method="class",data=dt[s==1,])
plot(ct.t)
text(ct.t, use.n=TRUE)
summary(ct.t)
tb.1 <- table(predicted=predict(ct.t,data=dt[s==1,], type="class"),
              observed=dt[s==1, "IWCT"])
tb.1 # Correctly classified (0-0,1-1), missclasified (0-1,1-0)
cc <-(tb.1[1,1]+tb.1[2,2])/(sum(tb.1))
cc # proportion correctly classified (correct classification rate)
mcc <-1-cc  # misclassificaiton rate 
mcc

## Durbin -Wattson test for autocorrelation
library(car)
dw <-glm(formula = IWCT~DEP_m+ACW_m+WID_m+BLD+FIN+LWDP, binomial,data=dt)
dwt(dw$residuals)
dwt(dw)
require(lmtest)
dwtest(dw)


## Lag analysis 
tsy<-ts(ctlogit$residuals)
plot(tsy)
hist(tsy)
summary(tsy)
plot(ctlogit$residuals~ctlogit$y)
lag.plot(tsy)
plot(ctlogit)
summary(ctlogit)

detach(dt)
