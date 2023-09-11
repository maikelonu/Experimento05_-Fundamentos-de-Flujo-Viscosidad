# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: PRINCIPIO BERNOULLI + PROPAGACION DE ERROR

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Cuantificacion_incertidumbre
# Propagacion_error
# Modelos_nls_(Nonlinear_Least_Squares)
# An?lisis_gr?fico_ggplot2
# //////////////////////////////////////////////////////////////////////////////////

# Working directory is selected
"/home/shoe/Downloads/LABHYD_Exp_04_Principio_Bernoulli-master"

# CRAN libraries are loaded
# require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(propagate)
require(reshape)
require(visreg)

# Workspace is cleared
#rm(list = ls())

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
# BLOCK: Data Input
# ////////////////////////////////////////////////////////

# Propagate vector containers are reset to NULL values
cont.mean <- NULL
cont.sd   <- NULL
cont.25   <- NULL
cont.975  <- NULL
cont.st.p <- NULL

# Propagate vector containers are defined
cont.mean <- vector()
cont.sd   <- vector()
cont.25   <- vector()
cont.975  <- vector()
cont.st.p <- vector()

# Observed variables are defined as vectors
# User MUST introduced values manually !!!
 dist_m <- c(0, 0.027, 0.056, 0.095, 0.127, 0.170) # distance from the origin (m)
base.CE <- c(122.5, 117.5, 37.5, 82.5, 90.0, 92.5) # observed values of carga estatica (mm)
base.CT <- c(123.0, 123.0, 123.0, 108.0, 103.0, 98.0) # observed values of carga total (mm)
PUERTO  <- c("puerto1","puerto2","puerto3","puerto4","puerto5","puerto6")

# Observed values are converted from mm to m
base.CE <- base.CE / 1000
base.CT <- base.CT / 1000

# Uncertainties are defined based on instruments standard deviation (m)
CE.sd <- 2.5 / 1000 # standard deviation for carga estatica measurement (m)
CT.sd <- 0.5 / 1000 # standard deviation for carga total measurement (m)

# Numeric counters are defined base on data.frame length
countermain <- length(base.CE)

# ////////////////////////////////////////////////////////
# BLOCK: Error propagation on pressure
# ////////////////////////////////////////////////////////

# Propagation loop over carga dinamica (CD) is initialized
for (i in 1 : countermain) {
  
  expr.CD <- expression(CT - CE) # expression for carga dinamica
  CE      <- c(base.CE[i], CE.sd) # vector including observed values and uncertainties
  CT      <- c(base.CT[i], CT.sd) # vector including observed values and uncertainties
  DF      <- cbind(CE, CT) # cbind function is requested to create a class-matrix object
  RES1    <- propagate(expr = expr.CD, # propagate function is requested
                       data = DF,
                       type = "stat",
                       do.sim = TRUE, # MonteCarlo Simulation = 5000 runs!
                       verbose = TRUE)
  
  # summary(RES1)
  
  cont.mean[i] <- as.numeric(RES1$sim[1])
  cont.sd[i]   <- as.numeric(RES1$sim[2])
  cont.25[i]   <- as.numeric(RES1$sim[5])
  cont.975[i]  <- as.numeric(RES1$sim[6])
  
  # A Shapiro-Wilks test is executed on propagation calculation
  ST <- shapiro.test(RES1$datSIM[1:4999])
  cont.st.p[i] <- ST$p.value
  
}

# Vectors are renamed
base.CD <- cont.mean
CD_sd   <- cont.sd
CD_25   <- cont.25
CD_975  <- cont.975
CD_shap_pvalue <- cont.st.p

# Uncertainties are calculated over the observed values according with
# their respective standard deviations
 CT_25 <- base.CT - (2 * CT.sd)
CT_975 <- base.CT + (2 * CT.sd)  
 CE_25 <- base.CE - (2 * CE.sd)
CE_975 <- base.CE + (2 * CE.sd)  

# Uncertainties are calculated over theoritical carga total values
    CT_TEO <- rep((base.CT[1]), countermain)
 CT_TEO_25 <- base.CT[1] - (2 * CT.sd)
CT_TEO_975 <- base.CT[1] + (2 * CT.sd)  
 CT_TEO_25 <- rep(CT_TEO_25, countermain)
CT_TEO_975 <- rep(CT_TEO_975, countermain)

# Headlosses (perdidas) are calculated
    PERD <- CT_TEO - base.CT
 PERD_25 <- PERD - (2 * CT.sd)
PERD_975 <- PERD + (2 * CT.sd)

# A data.frame container is created
df.press <- data.frame(PUERTO,
                       dist_m,
                       CT_TEO,
                       CT_TEO_25,
                       CT_TEO_975,
                       base.CT,
                       CT_25,
                       CT_975,
                       base.CE,
                       CE_25,
                       CE_975,
                       base.CD,
                       CD_25,
                       CD_975,
                       PERD,
                       PERD_25,
                       PERD_975,
                       CD_sd,
                       CD_shap_pvalue)

# df.press data.frame varibles are rounded to 5 decimals
df.press[, 2:19] <- round(df.press[, 2:19], 5)

# A new variable is introduced
df.press$Norm <- c("normal")

# Hypothesis test loop over Shapiro normality test is initialized
for(i in 1:countermain) {
  
  if (df.press$CD_shap_pvalue[i] > 0.05) {
    df.press$Norm[i] = "OK"
  }
  else {
    df.press$Norm[i] = "FAILED"
  }
  
}

# subset {base} function is requested to create a new data.frame
df.press.sub <- subset(df.press, select = c(CT_TEO,
                                            base.CT,
                                            base.CE,
                                            base.CD,
                                            PERD))

# melt {reshape} function is requested to convert data from "wide"
# to "long" format but extracting the first two column
df.press.melt <- melt(df.press.sub)

# First column is repeated 6 times and a data.frame is created
df.temp.presion1 <- as.data.frame(rep(df.press[,1],5))

# Second column is repeated 4 times and a data.frame is created
df.temp.presion2 <- as.data.frame(rep(df.press[,2],5))

# cbind {base} function is used over relevant data.frames
df.press.total <- cbind(df.temp.presion1, df.temp.presion2,df.press.melt)

# df.presion.long.total names are changed
names(df.press.total) <- c("PUERTO", "dist_m", "variable", "value")

# A ggplot2 object is created
fg01 <- ggplot() +
  geom_line(aes(x = dist_m,y = value,colour = variable),data=df.press.total,size = 0.75) +
  geom_point(aes(x = dist_m,y = value,shape = variable,colour = variable),data=df.press.total,size = 4.0) +
  geom_ribbon(aes(x = dist_m,ymin = CT_TEO_25,ymax = CT_TEO_975),data=df.press,fill = '#ff3333',alpha = 0.25) +
  geom_ribbon(aes(x = dist_m,ymin = CT_25,ymax = CT_975),data=df.press,fill = '#999900',alpha = 0.25) +
  geom_ribbon(aes(x = dist_m,ymin = CE_25,ymax = CE_975),data=df.press,fill = '#00cc00',alpha = 0.25) +
  geom_ribbon(aes(x = dist_m,ymin = PERD_25,ymax = PERD_975),data=df.press,fill = '#cc00cc',alpha = 0.25) +
  geom_ribbon(aes(x = dist_m,ymin = CD_25,ymax = CD_975),data=df.press,fill = '#00cccc',alpha = 0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Comportamiento de las componentes energeticas respecto de la distancia") +
  xlab("Distancia (m)") +
  ylab("Energia (m)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg01

# ////////////////////////////////////////////////////////
# BLOCK: Error propagation on velocity
# ////////////////////////////////////////////////////////

# Observed variables are defined as vectors
# User MUST introduced values manually !!!

# Time vector is introduced
TIMEV <- c(89.70, 91.89, 93.95, 94.07, 97.58, 98.04) # time (sec)

# Time-sd (Uncertainty) is calculated (s) due to natural variation
TIMEV.sd <- sd(TIMEV)

# Time-mean is calcculated
TIMEV.mean <- mean(TIMEV)

# Volume vector is calculated (m3)
# bonstant as it reaches a value of 9 Litres
VOLUMEV <- (9/1000)

# Volume-sd (Uncertainty) is calculated (m3)
# based on piezometer uncertainty
VOLUMEV.sd <- 0.0005

# Internal diameter vector is introduced (mm)
# constant for the GUNT Venturi tube
DIAM <- c(28.40, 22.50, 14.00, 17.00, 24.20, 28.40)

# Internal diameter vector is converted to (m)
DIAM <- DIAM/1000

# Area vector is calculated
AREAV <- 0.25*pi*(DIAM^2)

expr.Q.CALC <- expression(volumen / tiempo) # expression for caudal calculado
     tiempo <- c(TIMEV.mean, TIMEV.sd) # vector including observed values and uncertainties
    volumen <- c(VOLUMEV, VOLUMEV.sd) # vector including observed values and uncertainties
      DF.tv <- cbind(volumen, tiempo) # cbind function is requested to create a class-matrix object
       RES2 <- propagate(expr = expr.Q.CALC, # propagate function is requested
                  data = DF.tv,
                  type = "stat",
                  do.sim = TRUE, # MonteCarlo Simulation = 5000 runs!
                  verbose = TRUE)

Q.mean <- as.numeric(RES2$sim[1])
Q.sd   <- as.numeric(RES2$sim[2])
Q.25   <- as.numeric(RES2$sim[5])
Q.975  <- as.numeric(RES2$sim[6])

# Velocidad calculada is calculated (depends on continuity equation)
    V.CALC <- (Q.mean/AREAV)
 V.CALC.25 <- (Q.mean - 2 * Q.sd)/AREAV
V.CALC.975 <- (Q.mean + 2 * Q.sd)/AREAV

# Velocidad medida is calculated (dependes on carga dinamica CD)
    V.MEDIDA <- sqrt(2*9.81*base.CD)
 V.MEDIDA.25 <- sqrt(2*9.81*CD_25)
V.MEDIDA.975 <- sqrt(2*9.81*CD_975)

# A new data.frame is created based on velocity vectors
df.velocity <- data.frame(PUERTO,
                          dist_m,
                          V.CALC,
                          V.MEDIDA,
                          V.CALC.25,
                          V.CALC.975,
                          V.MEDIDA.25,
                          V.MEDIDA.975)

# df.velocity variables are rounded to 5 decimals
df.velocity <- round(df.velocity[, 2:8], 5)

# melt {reshape} function is requested to convert data from "wide"
# to "long" format but extracting the first two column
df.velocity.melt <- melt(df.velocity[ , 2:3])

# First column is repeated 6 times and a data.frame is created
df.temp.velocity1 <- as.data.frame(rep(df.press[,1],2))

# Second column is repeated 4 times and a data.frame is created
df.temp.velocity2 <- as.data.frame(rep(df.press[,2],2))

# cbind {base} function is used over relevant data.frames
df.velocity.total <- cbind(df.temp.velocity1, df.temp.velocity2,df.velocity.melt)

# df.velocity.total names are changed
names(df.velocity.total) <- c("PUERTO", "dist_m", "variable", "value")

# A ggplot2 object is created
fg02 <- ggplot() +
  geom_ribbon(aes(x = dist_m,ymin = V.MEDIDA.25,ymax = V.MEDIDA.975),data=df.velocity,fill = '#009999',alpha = 0.25) +
  geom_ribbon(aes(x = dist_m,ymin = V.CALC.25,ymax = V.CALC.975),data=df.velocity,fill = '#ff6666',alpha = 0.25) +
  geom_point(aes(x = dist_m,y = value,shape = variable,colour = variable),data=df.velocity.total,size = 3.5) +
  geom_line(aes(x = dist_m,y = value,colour = variable),data=df.velocity.total,size = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Comportamiento de la velocidad respecto de la distancia") +
  xlab("Distancia (m)") +
  ylab("Velocidad (m/s)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg02

# Mean absolute deviation (MAD) column is created
df.MAD <- as.data.frame(abs(df.velocity$V.CALC - df.velocity$V.MEDIDA))

# names function is requested
names(df.MAD) <- c("MAD")

# shapiro.test {stats} Normality Test is applied 
# if p-value > 0.05 normality stands true, meaning that
# the comparison is parametric
shapiro.test(df.MAD$MAD)

# Descriptive statistics are requested and round to 5 decimals
df.MAD.desc <- round(stat.desc(df.MAD),5)

# A new "CLASS" character column is created in data.frame df.MAD
df.MAD$CLASS <- c("Deviation")

# A ggplot2 object is created
fg03 <- ggplot(aes(x = CLASS,y = MAD),data=df.MAD) +
  geom_boxplot(size = 0.75,alpha = 0.70,
               outlier.colour = '#ff0033', outlier.size = 4.5) +
  geom_point(colour = '#0000ff', size = 4.5,
             position = position_jitter(width = 0.05)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot of Mean Absolute Deviation (MAD)") +
  xlab("Clase") +
  ylab("simple deviation (SD)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg03

# A diameter vector varible is added to df.velocity data.frame
df.velocity$DIAM <- DIAM

# ////////////////////////////////////////////////////////
# BLOCK: Regression Models
# ////////////////////////////////////////////////////////

# ------------------------------------------------
# nls EXPONENTIAL
# ------------------------------------------------

# A nls {stats} Nonlinear Least Squares model is fitted (exponential)
mod.exp <- nls(V.CALC ~ I(a * exp(b * DIAM)),
               data = df.velocity,
               start = list(a = 1, b = 0),
               trace = TRUE)

# A model summary is requested
summary(mod.exp)

# /////////////////////////////////////////////////
# lm SECTIONS:
# /////////////////////////////////////////////////
# 
# Residuals section:
# it provides a quick summary (min, 1Q, median, 3Q, max) of the distribution.
# 
# Coefficients section: each coefficient is a Gaussian random variable
# Estimate represents the mean distribution of the variable
# Std.Error displays the standard error of the variable
# the t value is Estimate divided by Std.Error
# the p value indicates the probability of getting a value larger than the t value
#
# Residual standard error outputs the standard deviation of residuals
#
# The degree of freedom indicates the differences between the observation 
# in training samples and the number used in the model
#
# Multiple R-squared is obtained by dividing the sum of squares.
#
# Adjusted R-squared uses an unbiased estimate, and will 
# be slightly less than multiple R-squared

# /////////////////////////////////////////////////

# ------------------------------------------------
# nls POTENTIAL
# ------------------------------------------------

# A nls {stats} Nonlinear Least Squares model is fitted (potential)
mod.pot <- nls(V.CALC ~ b * (DIAM ^ a),
               data = df.velocity,
               start = list(a = 1, b = 1),
               trace = TRUE)

# A model summary is requested
summary(mod.pot)

# Total sum of squares is calculated
TSS <- sum((df.velocity$V.CALC - mean(df.velocity$V.CALC))^2)

# ////////////////////////
# R2 for mod.exp
# ////////////////////////

# Residual sum of squares is calculated
RSS.p.exp <- sum(residuals(mod.exp)^2)

# R-squared is calculated
R2.exp <- (1 - (RSS.p.exp/TSS))

# R2 is requested
print(R2.exp)

# ////////////////////////
# R2 for mod.pot
# ////////////////////////

# Residual sum of squares is calculated
RSS.p.pot <- sum(residuals(mod.pot)^2)

# R-squared is calculated
R2.pot <- (1 - (RSS.p.pot/TSS))

# R2 is requested
print(R2.pot)

# A ggplot2 object is created
fg04 <- ggplot(aes(x = DIAM,y = V.CALC),data = df.velocity) +
  geom_point(colour = '#ff0033',size = 5.5,alpha = 0.90) +
  geom_smooth(aes(x = DIAM,y = V.CALC),
              data = df.velocity, method = loess, se = FALSE,
              size = 0.75, alpha = 0.50) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Comportamiento de la velocidad calculada respecto del diametro") +
  xlab("Diametro (mm)") +
  ylab("Velocidad Calculada (m/s)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg04

# round_df function is applied to relevant data.frames
df.press <- round_df(df=df.press, digits=3)
df.velocity <- round_df(df=df.velocity, digits=3)
df.MAD.desc <- round_df(df=df.MAD.desc, digits=3)

# Objects to export
# fg01, fg02, fg03, fg04
# df.press, df.velocity, df.MAD.desc
# shapiro.test(df.MAD$MAD)
# summary(mod.exp), summary(mod.pot), R2.exp, R2.pot
write.csv(df.press, file = "df.press.csv")
write.csv(df.velocity, file = "df.velocity.csv")
write.csv(df.MAD.desc, file = "df.MAD.desc.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
