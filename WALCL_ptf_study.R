#Clearing the environment 
rm(list = ls())

#Charging the libraries
library(urca)
library(tseries)
library(xts)
library(xtable)
library(tsDyn)
library(fGarch)
library(FitARMA)
library(fUnitRoots)
library(vars)
library(maps)

library("readxl")
library(tidyverse)
library(ggplot2)
library(plotly)
library(dplyr)
library(quantmod)
library(tseries)

library(forecast)
library(vars)
library(astsa)
library(xts)
install.packages("stargazer")
library(stargazer)

#Loading the dataframe from a xlsx file
my_data <- read_excel("D:\\Master 1 Eco stats\\Projets\\Semestre 2\\Econofi projet supp\\Base de données\\WALCL.xls")
head(my_data)

#1 - Study of each series seperately: transforming into log-series, then into time series + identification of each series' elements + ARIMA specification
#Plotting both series
ggplot(my_data, aes(Date, WALCL)) + geom_line() 
ggplot(my_data, aes(Date, Portefeuille)) + geom_line() 

#To stabilize each series' variance, transform the original series into a log series by applying a log-transformation
my_data$logWALCL <- log(my_data$WALCL)
my_data$logPortefeuille <- log(my_data$Portefeuille)
head(my_data)
#Plotting the log-transformed series
ggplot(my_data, aes(Date, logWALCL)) + geom_line() 
ggplot(my_data, aes(Date, logPortefeuille)) + geom_line() 

#Transforming the series into time series, frequency = 52 since we have weekly data. 
tswalcl <- ts(my_data$WALCL, frequency = 52)
tsPortefeuille <- ts(my_data$Portefeuille,frequency = 52)

#Plotting the multiplicative decomposition 
decWALCL <- decompose(tswalcl, type="multiplicative", filter=NULL)
decPortefeuille <- decompose(tsPortefeuille, type="multiplicative", filter=NULL)
plot(decWALCL)
plot(decPortefeuille)

#Adding x-axis timestamps with the ggplot2 library
n = length(my_data$Date)
df1 = data.frame(date = my_data$Date, name = rep("random", n), data = as.numeric(decWALCL$random))
df2 = data.frame(date = my_data$Date, name = rep("seasonal", n), data = as.numeric(decWALCL$seasonal))
df3 = data.frame(date = my_data$Date, name = rep("trend", n), data = as.numeric(decWALCL$trend))
df4 = data.frame(date = my_data$Date, name = rep("observed", n), data = as.numeric(decWALCL$x))
dfWALCL = rbind(df1, df2, df3, df4)


ggplot(dfWALCL, aes(date, data)) +
  facet_grid(rows=vars(name), scales ="free_y")+
  geom_line() +
  theme_bw() +
  labs(y="Value", x="Date") +
  ggtitle(expression(Decomposition~of~multiplicative~WALCL~series)) +
  theme(plot.title=element_text(hjust=0.5))

n = length(my_data$Date)
df11 = data.frame(date = my_data$Date, name = rep("random", n), data = as.numeric(decPortefeuille$random))
df22 = data.frame(date = my_data$Date, name = rep("seasonal", n), data = as.numeric(decPortefeuille$seasonal))
df33 = data.frame(date = my_data$Date, name = rep("trend", n), data = as.numeric(decPortefeuille$trend))
df44 = data.frame(date = my_data$Date, name = rep("observed", n), data = as.numeric(decPortefeuille$x))
dfPortefeuille = rbind(df11, df22, df33, df44)


ggplot(dfPortefeuille, aes(date, data)) +
  facet_grid(rows=vars(name), scales ="free_y")+
  geom_line() +
  theme_bw() +
  labs(y="Value", x="Date") +
  ggtitle(expression(Decomposition~of~multiplicative~Portefeuille~series)) +
  theme(plot.title=element_text(hjust=0.5))

#remove the intermediary series that are not used anymore
rm(df1, df2, df3, df4, df11, df22, df33, df44)

#Decomposition pt2: making the regression xt = a + bt + St + ut / need to add start + end dates (i.e., start=c(2002,50),end=c(2022,13))? 
#Need to verify: do we need to add the dates as the following format -- xWALCL <- ts(my_data$WALCL, start=c(2002,50),end=c(2022,13), frequency = 52)
#For the WALCL series
xWALCL <- ts(my_data$WALCL, frequency = 52)
WALCL_seas <- decWALCL$seasonal
WALCL_trend <- decWALCL$trend
xWALCL <- window(xWALCL, frequency = 52)
xWALCL_trend <- window(WALCL_trend, frequency = 52)
xWALCL_seas <- window(WALCL_seas, frequency = 52)
xWALCL_rand <- 
reg_WALCL <- lm(xWALCL~xWALCL_trend*(1+xWALCL_seas))
summary1 <- summary(reg_WALCL) #gets the full regression as we want
WALCL_trend <- decWALCL$figure
stargazer(reg_WALCL, type="latex",report=('vc*p'))

#For the Portefeuille series - Didn't do the date on this one, maybe don't need to 
xPtf <- ts(my_data$Portefeuille, frequency = 52)
Ptf_seas <- decPortefeuille$seasonal
Ptf_trend <- decPortefeuille$trend
xPtf <- window(xPtf, frequency=52)
xPtf_trend <- window(Ptf_trend,frequency=52)
xPtf_seas <- window(Ptf_seas, frequency=52)
reg_Ptf <- lm(xPtf~xPtf_trend*(1+xPtf_seas))
summary(reg_Ptf)
Ptf_trend <- decPortefeuille$figure
stargazer(reg_Ptf, type="latex",report=('vc*p'))


#ARMA specification 
#Testing the unit root: both series are containing a unit root, need to 1st differenciate 
print(adf.test((my_data$logWALCL))) 
print(adf.test((my_data$logPortefeuille))) 
print(adf.test(diff(my_data$logWALCL))) 
print(adf.test(diff(my_data$logPortefeuille))) 

#Plotting ACFs + PACFs to determine AR and MAs orders on both log-transformed series 
acf2(diff(my_data$logWALCL)) #Expliquer POURQUOI ON A UN 2???
acf2(diff(my_data$logPortefeuille))

#ARIMA model
#auto.arima finds the best model with or without drift
auto.arima(my_data$logPortefeuille, trace=TRUE) 
auto.arima(my_data$logWALCL, trace=TRUE) 
#Estimation of the ARIMA model based on earlier results + plots whether the residual is a WN or not
sarima(my_data$logPortefeuille,  p = 1, d =1, q = 0)
sarima(my_data$logWALCL,  p = 5, d =1, q = 0)

#Doing the same for differenciated series, becomes respectively AR(1) and AR(5) - not necessary, it is to compare with previous results 
auto.arima(diff(my_data$logPortefeuille), trace=TRUE) 
auto.arima(diff(my_data$logWALCL), trace=TRUE) 
sarima(diff(my_data$logPortefeuille),  p = 1, d =0, q = 0)
sarima(diff(my_data$logWALCL),  p = 5, d =0, q = 0)


#2 - Regards to each series' attributes (stationnary/not, seasonnalities etc), study of a VAR model
#Using the log-transformed data, creating a new dataframe that we will use for our study 
data <- cbind(diff(my_data$logWALCL), diff(my_data$logPortefeuille))
colnames(data) <- c("d(logWALCL)", "d(logPortefeuille)") #Since a VAR model needs to be stationnary, we're using the differenciated series.
head(data)

#Selecting the best VAR order
VARselect(data) #AIC and FPE select 9th order, while HQ & SC select 5th order

#Estimating the VAR model
model <- vars::VAR(y=data, p=9) #Since AIC + FPE select p = 9 
summary(model)


#PROPOSITION DE MODELE - en date du 22/05/2022
#Using data, transform both series into time series with the ts() function
#Starting date: 50th week of 2002 ; Ending point: 13th week of 2022
data_diffts <- ts(data=data, start=c(2002,50), end=c(2022,13), frequency=52, names=c("dlwalcl","dlptf"))
VARselect(data_diffts)$selection #to quickly get to the orders chosen by R: same as before, p=9 for AIC/FPE and p=5 for HQ/SC
varendo <- data.frame(data_diffts) #transform into a dataframe, we don't really care about this detail
#Estimating both models: either with p=9 as chosen by AIC/FPE or with p=5 as chosen by HQ/SC
model_bis1 <- VAR(varendo, type=c("both"),ic=c("AIC")) #type both because we got a constant + trend, as determined previously through the additive decomposition
summary(model_bis1) #p = 9
residvar_1bis <- resid(model_bis1) #residuals model 1 
model_bis2 <- VAR(varendo, type=c("both"),ic=c("HQ"))
summary(model_bis2) #p = 5
residvar_2bis <- resid(model_bis2) #residual model 2

#Cointegration properties: if series are cointegrated - i.e., common trend, in need for a VECM model 
#We need to take the original series, i.e., before the 1st difference: log-transformed series 
#Transformation of logs into time series 
tslwalcl <- ts(my_data$logWALCL, start=c(2002,50), end=c(2022,13), frequency = 52)
tslPortefeuille <- ts(my_data$logPortefeuille, start=c(2002,50), end=c(2022,13), frequency = 52)
data_ts <- cbind(tslwalcl,tslPortefeuille)
colnames(data_ts) <- c("logWALCL", "logPortefeuille")
data_ts #dataframe that contains the log series w/o differenciation

#Computing the Johansen test - eigenvalues/valeurs propres or trace 
#The K stands for the lag order of the series (levels) in the VAR - Here I chose 9 but I'm not sure about my choice, maybe it is K=2 only
jotest <- ca.jo(data_ts, type="trace", K=9, ecdet="none", spec="longrun")
summary(jotest)
eigvpm <- ca.jo(data_ts, type="eigen", K=9, ecdet="none", spec="longrun")
summary(eigvpm) #I can't seem to understand the output - see later on - rank of 1 

#Computing the VECM model - Not sure about the 
vecm <- VECM(data_ts, lag=9, r=1, include = "both", beta= NULL, estim="ML", LRinclude="const") #max likelihood estimator as usual
summary(vecm)
stargazer(vecm$model, type="latex")

sjd.vecm <- ca.jo(data_ts,ecdet="const",type="eigen", K=9, spec="longrun") 
summary(sjd.vecm)
#ENDING



#3 - Forecasting 2 period of time ahead
#Based on our study, we find that the VAR modelling fits better to get a prediction on 2 periods ahead - means 2 weeks normally 
var.pred <- predict(model_bis1, n.ahead = 2) 
plot(var.pred)
fanchart(var.pred) #same plot as before but with shaded colors + add colors 

#-----------------------------------------------------------------------
#COMMENT: 
#Possible to use the VECM model to compute the forecast 2 periods ahead? 
#-----------------------------------------------------------------------


#4 - Orthogonal impulse response with confidence interval analysis
irfdlWALCLrep <- irf(model, n.ahead=10, boot=TRUE, ci=0.95,ortho=TRUE) #OK to comment - do we need to shift the model? 
plot(irfdlWALCLrep)

#Granger causality
grangerDlogWALCL <- causality(model, cause="d.logWALCL.")
grangerDlogPortefeuille <- causality(model, cause="d.logPortefeuille.")

grangerDlogWALCL

grangerDlogPortefeuille





#------------------------------------------------------------------------------
#DID THIS EARLIER

data <- cbind(my_data$logWALCL, my_data$logPortefeuille)
colnames(data) <- c("logWALCL", "logPortefeuille")

# Cointegration properties

# Cointegration rank test (Johansen)
jotest=ca.jo(data, type="trace", K=9, ecdet="none", spec="longrun")
summary(jotest)

sjd.vecm <- ca.jo(data, ecdet="const", type="eigen", K=9, spec="longrun")
summary(sjd.vecm)

# Estimation of the VECM
vecm <- VECM(data, lag=9, r=1, include = "both", beta= NULL, estim="ML", LRinclude="const") 
summary(vecm)

# Plot residuals of the VECM
plotres(sjd.vecm)

#-------------------------------------------------------------------------------

# Granger causality from the VAR in level: ONLY IN THE  CASE WHERE COINTEGRTAION HAS BEEN VALIDATED


causality(model, cause = "d.logWALCL.", vcov.=NULL, boot=FALSE, boot.runs=100)
causality(model, cause = "d.logPortefeuille.", vcov.=NULL, boot=FALSE, boot.runs=100)
