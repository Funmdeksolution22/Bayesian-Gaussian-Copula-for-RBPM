library(sp)
library(spdep)
library(sf)
library(mgcv)
library(tidyverse)
library(GJRM)
library(BayesX)
library(ggplot2)
library(VGAM)
library(spDataLarge)
library(copula)
install.packages("spDataLarge", dependencies = TRUE)
data<-read.csv('C:\\Users\\14483065\\Documents\\Spatial Statistics Document\\main_extract.csv', header=TRUE)
map <- read.bnd('C:\\Users\\14483065\\Documents\\Spatial Statistics Document\\nigeria.bnd')
xt <- list(polys = map)
data$state <- factor(data$state)
Data = na.omit(data)
n=length(data$mobidity)
library(stats)
library(jtools)
library(randomForest)
rf_model <- randomForest(mobidity ~ gm + prim + sec + high + Poorer + Middle + Richer + Richest + Male + never_breast + still_breast + 
                           priv_hospital + govt_hospital + home + hausa + igbo + yoruba + islam + christian +  
                           underweight + normalweight + overweight + sstate, data = data, ntree = n)
oob_error <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
n <- nrow(Data)
p <- length(rf_model$forest$xlevels) - 1
print(aic_like)
######mobidity factor######
data = list(data$mortality, data$mobidity)
modelo_challenger <- glm(formula = mobidity ~ . -state -region -mortality - gm-gm2 -sstate,
                         data = Data,
                         family = binomial(link = "probit"))
modelo_challenger
#Step wise for mobidity for examle 
step_challenger <- step(object = modelo_challenger,
                        k = qchisq(p = 0.05, df = 1, lower.tail = FALSE))

library(xtable)
summ(model = step_challenger, confint = T, digits = 4, ci.width = 0.95) # This should give you the final covariates for each response variable
####################################################################

######mortality ########################
modelo_challenger <- glm(formula = mortality ~ . -state -region - mobidity -regular_antenatal - irregular antenatal - gm-gm2 -sstate,
                         data = Data,
                         family = binomial(link = "probit"))
modelo_challenger
#Step wise for mortality for example 
step_challenger <- step(object = modelo_challenger,
                        k = qchisq(p = 0.05, df = 1, lower.tail = FALSE))
library(xtable)
summ(model = step_challenger, confint = T, digits = 4, ci.width = 0.95) # 
########################################### Model analysis###########################################
eq1 <- mortality  ~ gm  + Male + prim + sec + high + Poorer + Middle + Richer + Richest  + still_breast + irregular_breast + priv_hospital + govt_hospital + home + hausa + igbo + yoruba + islam + christian +  underweight + normalweight + overweight + mobidity +  s(state, bs = "mrf", xt = xt)
eq2 <- mobidity  ~ gm  + Male + prim + sec + high + Poorer + Middle + Richer + Richest  + still_breast + irregular_breast + priv_hospital + govt_hospital + home + hausa + igbo + yoruba + islam + christian +  underweight + normalweight + overweight +s(state, bs = "mrf", xt = xt)


eq3 <- ~  #gm  +  + prim + sec + high + Poorer + Middle + Richer + Richest + Male + never_breast + still_breast + 
  #priv_hospital + govt_hospital + home + hausa + igbo + yoruba + islam + christian + 
  #underweight + normalweight + overweight
    s(state, bs = "mrf", xt = xt)


start_time <- Sys.time()


f.l <- list(eq1 ,eq2)


outDA <- gjrm(f.l,margins = c("probit","probit"), copula="N", model="BPO", data=Data, Chol = TRUE)
end_time <- Sys.time()

conv.check(outDA)
post.check(outDA)
AIC(outDA)
BIC(outDA)
result1 = summary(outDA$gam2)
reslibrary(xtable)
result <- print(result1, digits = 3, prob = c(0.025, 0.975))
xtable(summary(outDA$gam1))
mortality = (summary(outDA$gam1), digit = 4)
summary(outDA$gam3)

par(mfrow = c(1, 3))
plot(outDA,eq1,xlab=" ", ylab="Estimated effect of Child´s Morbidity", ylim = c(-5, 5))
plot(outDA,eq=2,xlab=" ", ylab="Estimated effect of Child´s Morbidity", ylim = c(-5, 5))
plot(outDA,eq=3,xlab=" ", ylab="Estimated effect of Child´s Morbidity", ylim = c(-5, 5))

tidy_summary <- tidy(outDA$gam1)
xtable_summary <- xtable(tidyr::spread(tidy_summary, key = term, value = estimate))


 