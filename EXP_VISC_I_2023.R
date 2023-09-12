# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: VISCOSIDAD CINEMATICA

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Pruebas_Concordancia (Accuracy, Precision, CCC, TDI, CP)
# t-test
# Shapiro_test
# Estadisticas_descriptivas???
# Modelos_potenciales no-lineales 
# Prediccion_(CI + PI)
# Analisis_grafico_ggplot2
# //////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
rm(list = ls())

# Working directory is selected
setwd("/home/shoe/Downloads/LABHYD_Exp_05_Viscosidad_Cinamatica-master")

# CRAN libraries are loaded
# require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(tidyr)
require(reshape)
require(reshape2)
require(visreg)
require(proto)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -3)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Custom Function {as.lm.nls}
# ////////////////////////////////////////////////////////

as.lm.nls <- function(object, ...) {
  require(proto)
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:",
               paste(class(object), collapse = " "))
    warning(w)
  }
  
  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }
  
  response.name <- if (length(formula(object)) == 2) "0" else
    as.character(formula(object)[[2]])
  
  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name
  
  fo <- sprintf("%s ~ %s - 1", response.name,
                paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = as.proto(   as.list(L)   ))
  
  do.call('lm', list(fo, offset = substitute(fitted(object))))
}
# ////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////
# BLOCK: reading input data
# ////////////////////////////////////////////////////////

# Agreement function source is loaded
options(scipen=999)
source("agreement_2020.R")

# Input data is loaded and a data.frame is created
df.base <- read.table("viscosity_SAE_2017_02.txt", header = TRUE)

# Desc {DescTools} function is requested
Desc(df.base, plotit = TRUE)

# names {base} function is requested
names(df.base)

# A ggplot2 object is created
fg01 <- ggplot() +
  geom_point(aes(x = temp_C,y = viscosity_cSt,shape = group,colour = column),data=df.base,
             size = 7.5,alpha = 0.95,position = position_jitter(width = 5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Oil viscosity vs temperature") +
  xlab("Temperature (C)") +
  ylab("Viscosity (cSt)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg01

# A ggplot2 object is created
fg02 <- ggplot() +
  geom_boxplot(aes(y = viscosity_cSt,x = column,colour = column),data=df.base,
               size = 0.85,alpha = 0.95) +
  geom_point(aes(x = column,y = viscosity_cSt,shape = group,colour = column),
             data=df.base,size = 4.5,alpha = 0.95) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Oil viscosity boxplot (by column)") +
  xlab("Column (-)") +
  ylab("Viscosity (cSt)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg02

# A ggplot2 object is created
fg03 <- ggplot() +
  geom_boxplot(aes(y = viscosity_cSt,x = fluid,colour = fluid),data=df.base,
               size = 1.25,alpha = 0.95,outlier.colour = '#ffffff',outlier.size = 0) +
  geom_point(aes(x = fluid,y = viscosity_cSt,shape = group,colour = column),
             data=df.base,size = 5.5,position = position_jitter(width = 0.1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Oil viscosity boxplot") +
  xlab("Column (-)") +
  ylab("Viscosity (cSt)") +
  theme_bw(base_size = 22.0) 

# A ggplot2 object is requested
fg03

# ////////////////////////////////////////////////////////
# BLOCK: Models
# ////////////////////////////////////////////////////////

# A nls {stats} model is fitted (power)
mod.nls.power <- nls(viscosity_cSt ~ a*(temp_C^power),
                     data = df.base,
                     start = list(a=0.5, power = 1),
                     trace = T)

# A summary function is requested
summary(mod.nls.power)

#------------------------------------------
# R2 formula
#------------------------------------------

# Residual sum of squares is calculated 
RSS.p1 <- sum(residuals(mod.nls.power)^2)

# Total sum of squares is calculated 
TSS <- sum((df.base$viscosity_cSt - mean(df.base$viscosity_cSt))^2)

# R-squared is calculated 
R2.nls <- (1 - (RSS.p1/TSS))

# R2 is requested 
print(R2.nls)

# Based on R2 results, nls is selected

# CIs are calculated using custom function {as.lm.nls}
df.predCI <- as.data.frame(predict(as.lm.nls(mod.nls.power),
                                   interval = 'confidence',
                                   level = 0.95))

# df.predCI data.frame is rounded to 3 decimals
df.predCI <- round(df.predCI,3)

# cbind {base} function is used to join data.frames
df.output <- cbind(df.base,df.predCI)

# confint2 {nlstools} function is requested to 
# confidence intervals in nonlinear regression coefficients
confint2(mod.nls.power)

# A ggplot object is created
fg04 <- ggplot() +
  geom_point(aes(x = temp_C,y = viscosity_cSt, shape = group,colour = column),
             data=df.output, size = 7.5, position = position_jitter(width = 5)) +
  geom_line(aes(x = temp_C,y = fit),
            data=df.output,colour = '#0000ff',size = 0.95,alpha = 0.9) +
  geom_line(aes(x = temp_C,y = lwr),
            data=df.output,colour = '#666666',size = 0.95,linetype = 2,alpha = 0.9) +
  geom_line(aes(x = temp_C,y = upr),
            data=df.output,colour = '#666666',size = 0.95,linetype = 2,alpha = 0.9) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Oil viscosity vs temperature") +
  xlab("Temperature (C)") +
  ylab("Viscosity (cSt)") +
  theme_bw(base_size = 22.0)

# A ggplot object is requested
fg04

#-----------------------------------------------
# Agreement comparison, same-type instruments
#-----------------------------------------------

# subsetting data.frame for 300-columns
df.300 <- subset(df.base, column == "300_A" | column == "300_B")

# selecting relevant data.frame variables
df.300 <- df.300[c("viscosity_cSt","column")]

# A rep function is introduced
df.300$ID <- rep(seq(from=1, to=(length(df.300$viscosity_cSt)/2), by=1), 2)

# dcast function is requested to convert from long to wide format
# (unmelt) NAs are introduced by coercion !!!!!!!!!!!!!!!!!
df.300 <- dcast(df.300,
                ID ~ column,
                fun.aggregate = sum,
                value.var="viscosity_cSt")

# data.frame names are changed
names(df.300) <- c("ID","C_300_A","C_300_B")

# A simple-deviation column is created
df.300$SD <- (abs(df.300$C_300_A - df.300$C_300_B))

# the mean of the SD is calculated
mean(df.300$SD)

# An agreement {Agreement} function is executed
test300 <- agreement(x = df.300$C_300_A,
                     y = df.300$C_300_B,
                     error = "const",
                     target = "random",
                     CCC_a = 0.975, # no more of 2.5% should be accepted for instruments comparison
                     TDI_a = 1.00, # an absolute difference when error="const" (units-interval)
                     alpha = 0.05, # 100(1-alpha)
                     CP_a = 0.90, # mirrored of TDI usually taken as 0.90 
                     H_label = "C_300_A",
                     V_label = "C_300_B",
                     dec = 3)

# An agreement-list summary is requested
summary.agreement.Allowance(test300)
summary.agreement.Estimate(test300)
summary.agreement.Conf_Limit(test300)

# subsetting data.frame for 150-columns
df.150 <- subset(df.base, column == "150_A" | column == "150_B")

# selecting relevant data.frame variables
df.150 <- df.150[c("viscosity_cSt","column")]

# A rep function is introduced
df.150$ID <- rep(seq(from=1, to=(length(df.150$viscosity_cSt)/2), by=1), 2)

# dcast function is requested to convert from long to wide format
# (unmelt) NAs are introduced by coercion !!!!!!!!!!!!!!!!!
df.150 <- dcast(df.150,
                ID ~ column,
                fun.aggregate = sum,
                value.var="viscosity_cSt")

# data.frame names are changed
names(df.150) <- c("ID","C_150_A","C_150_B")

# A simple-deviation column is created
df.150$SD <- (abs(df.150$C_150_A - df.150$C_150_B))

# the mean of the SD is calculated
mean(df.150$SD)

# An agreement {Agreement} function is executed
test150 <- agreement(x = df.150$C_150_B,
                     y = df.150$C_150_A,
                     error = "const",
                     target = "random",
                     CCC_a = 0.975, # no more of 2.5% should be accepted for instruments comparison
                     TDI_a = 1.00, # an absolute difference when error="const" (units-interval)
                     alpha = 0.05, # 100(1-alpha)
                     CP_a = 0.90, # mirrored of TDI usually taken as 0.90 
                     H_label = "C_150_A",
                     V_label = "C_150_B",
                     dec = 3)

# An agreement-list summary is requested
summary.agreement.Allowance(test150)
summary.agreement.Estimate(test150)
summary.agreement.Conf_Limit(test150)

# round_df function is applied to relevant data.frames
df.output <- round_df(df=df.output, digits=5)

#//////////////////////////////////////////////////////////////////
# Objects to export:
# Desc(4column, 5group)
# fg02, fg03, fg04
# summary (mod.nls.power)
# R2.nls
# test300
# test150
# df.output
write.csv(df.output, file = "df.output.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
