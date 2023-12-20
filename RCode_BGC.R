library(rstan)
library(ggplot2)
library(ggpubr)
library(readxl)
library(R2BayesX)
library(tidyverse)
library(Matrix)
require(rgdal)
setwd("C:/Users/14483065/OneDrive/Desktop/Spatial")
secd = read.csv('main_extract.csv', header=TRUE)
##secd$state
y1 <- secd$mortality
y2 <- secd$mobidity
N=length(y1)
W <- read.csv("https://raw.githubusercontent.com/eosafu/Model/master/Nig_adj_Matrix.txt",sep=",",header = F) # Nigeria adjacency matrix
colnames(W) <- NULL
W <- as(as.matrix(W[,-1]),"sparseMatrix")
image(W)
# Compute spatial effect precision matrix
P <- as(diag(rowSums(W)),"sparseMatrix")
q <- nrow(P)
rho <- 0.99
D <- P-rho*W
D_W_inv = solve(D) # enter on stan as data
K=secd$state
y=c(y1, y2)
##X1<-secd$prim
##X2<-secd$urban
##X3=secd$sec
x<-as.matrix(data.frame(X1=secd$Male, X2=secd$prim, X3=secd$sec, X4=secd$gm, X5=secd$high, X6=secd$Poorer, x7=secd$Middle,
                        x8=secd$Richer, x9=secd$Richest, x10=secd$still_breast, x11=secd$irregular_breast, x12=secd$priv_hospital, 
                        x13=secd$govt_hospital, x14=secd$home, x15=secd$hausa, x16=secd$igbo, x17=secd$yoruba, x18=secd$islam, 
                        x19=secd$christian, x20=secd$underweight, x21=secd$normalweight, x22=secd$overweight))
z<-as.matrix(data.frame(x1=secd$Male, X2=secd$prim, X3=secd$sec, X4=secd$gm, X5=secd$high, X6=secd$Poorer, x7=secd$Middle,
                        x8=secd$Richer, x9=secd$Richest, x10=secd$still_breast, x11=secd$irregular_breast, x12=secd$priv_hospital, 
                        x13=secd$govt_hospital, x14=secd$home, x15=secd$hausa, x16=secd$igbo, x17=secd$yoruba, x18=secd$islam, 
                        x19=secd$christian, x20=secd$underweight, x21=secd$normalweight, x22=secd$overweight))
##x<-as.matrix(data.frame(X1=secd$prim, X2=secd$urban, X3=secd$sec))
##z<-as.matrix(data.frame(X1=secd$prim, X2=secd$urban, X3=secd$sec))
data_list <- list(N=N, y=y, x=x, z=z, K=K, mu=rep(0, 37), D_W_inv=as.matrix(D_W_inv))
stan_model <- stan(file = "Main_code.stan", 
                   data = data_list,
                   iter = 5000,
                   warmup = 1000,
                   thin = 4,
                   cores=4,
                   chains = 4,
                   algorithm = "NUTS")
result <- print(stan_model, digits = 4, prob = c(0.025, 0.975))
traceplot(stan_model, pars = c("beta1", "beta10", "alpha"), inc_warmup = TRUE, nrow = 6, xlim = c(1000, 5000))
traceplot(stan_model, pars = c("beta2", "beta20"), inc_warmup = TRUE, nrow = 6, xlim = c(1000, 5000))
traceplot(stan_model, pars = "phi1", inc_warmup = TRUE, nrow = 6, xlim = c(1000, 5000)) ### traceplot of mortality
traceplot(stan_model, pars = "phi2", inc_warmup = TRUE, nrow = 6, xlim = c(1000, 5000)) #### traceplot of morbidity 
traceplot(stan_model, pars = "phi3", inc_warmup = TRUE, nrow = 6, xlim = c(1000, 5000)) #### the traceplot of the joint effect
#### The spatial plot of the child mortatlity #######
phi1 <- stan_model@.MISC[["summary"]][["msd"]][47:83, 1]
nigeria_map <- st_read("C:\\Users\\14483065\\Documents\\Spatial\\nigeria_shp\\gadm36_NGA_1.shp")
require(sf)
require(ggplot2)
require(plotly)
ggplot(nigeria_map) + geom_sf(aes())+
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "yellow", high = "blue") +
  labs(title = " ")
nigeria_map$phi2 <- phi2
mapview(nigeria_map, zcol = 'phi2')

# Plot 02
g <- ggplot(nigeria_map) + geom_sf(aes(fill = phi2))
ggplotly(g)


#### The spatial plot of the child morbidity #######
phi2 <- stan_model@.MISC[["summary"]][["msd"]][84:120, 1]
nigeria_map <- st_read("C:\\Users\\14483065\\Documents\\Spatial\\nigeria_shp\\gadm36_NGA_1.shp")
require(sf)
require(ggplot2)
require(plotly)
ggplot(nigeria_map) + geom_sf(aes())+
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "yellow", high = "blue") +
  labs(title = " ")
nigeria_map$phi2 <- phi2
mapview(nigeria_map, zcol = 'phi2')

# Plot 02
g <- ggplot(nigeria_map) + geom_sf(aes(fill = phi2))
ggplotly(g)


#### The spatial plot of the Joint effect of mortality and morbidity #######
phi3 <- stan_model@.MISC[["summary"]][["msd"]][121:157, 1]
nigeria_map <- st_read("C:\\Users\\14483065\\Documents\\Spatial\\nigeria_shp\\gadm36_NGA_1.shp")
require(sf)
require(ggplot2)
require(plotly)
ggplot(nigeria_map) + geom_sf(aes())+
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "yellow", high = "blue") +
  labs(title = " ")
nigeria_map$phi3 <- phi3
mapview(nigeria_map, zcol = 'phi3')

# Plot 02
g <- ggplot(nigeria_map) + geom_sf(aes(fill = phi3))
ggplotly(g)