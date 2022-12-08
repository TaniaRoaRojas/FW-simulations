#### Distancia L2

library(sfsmisc)
library(ggplot2)       #
library(scales)        # working with ggplot2 for label formatting
library(gridExtra)     # working with ggplots for arranging plots
library(ggthemes)      # clean theme for ggplot2
library(viridis)       # color palette
library(DT)
library(statip)
library(plot.matrix)

#######FUNCIONES
##Parámetros
N = 1000
x = seq(0.01, 0.99, length.out = 100)
x0 = seq(0.01, 0.99, length.out = 11)
t0 = c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
t = t0/(N) 
t



#2 beta con parametros adecuados
##beta parametrizada como media y varianza

media = matrix(NA, ncol= 11, nrow = 11)
for (i in 1:11) {
  media[i,1] = x0[1]
  media[i,2] = x0[2]
  media[i,3] = x0[3]
  media[i,4] = x0[4]
  media[i,5] = x0[5]
  media[i,6] = x0[6]
  media[i,7] = x0[7]
  media[i,8] = x0[8]
  media[i,9] = x0[9]
  media[i,10] = x0[10]
  media[i,11] = x0[11]
}

colnames(media) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(media) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
media

variance = matrix(NA, ncol= 11, nrow = 11)
for (i in 1:11) {
  variance[i,1] = x0[1]*(1-x0[1])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,2] = x0[2]*(1-x0[2])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,3] = x0[3]*(1-x0[3])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,4] = x0[4]*(1-x0[4])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,5] = x0[5]*(1-x0[5])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,6] = x0[6]*(1-x0[6])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,7] = x0[7]*(1-x0[7])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,8] = x0[8]*(1-x0[8])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,9] = x0[9]*(1-x0[9])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,10] = x0[10]*(1-x0[10])*(1 - (1 - (1/(2*N)))^t0[i])
  variance[i,11] = x0[11]*(1-x0[11])*(1 - (1 - (1/(2*N)))^t0[i])
}


colnames(variance) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(variance) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
variance

#verificar
variance < media*(1-media)


#####################################################
#####################################################
#####################################################

### Hellinger 3 - Simpson
library('Bolstad2')
## integrate the normal density from -3 to 3
#x=seq(-3,3,length=100)
#fx=dnorm(x)
#estimate=sintegral(x,fx, n.pts = 500)$int
#true.val=diff(pnorm(c(-3,3)))

##función de la raíz cuadrada dentro de la distancia de hellinger
square_hellin = function(x,y,z){
  sh = sqrt(y*z)
  res = sintegral(x, sh)$int
  return(res)
}

###hellinger asymtptical expansion

AEp0 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp0[,i] = expan(x, x0[1], N, t[i]) / (sintegral(x, expan(x, x0[1], N, t[i]))$int)
}

AEp01 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp01[,i] = expan(x, x0[2], N, t[i]) / (sintegral(x, expan(x, x0[2], N, t[i]))$int)
}

AEp02 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp02[,i] = expan(x, x0[3], N, t[i]) / (sintegral(x, expan(x, x0[3], N, t[i]))$int)
}

AEp03 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp03[,i] = expan(x, x0[4], N, t[i]) / (sintegral(x, expan(x, x0[4], N, t[i]))$int)
}

AEp04 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp04[,i] = expan(x, x0[5], N, t[i]) / (sintegral(x, expan(x, x0[5], N, t[i]))$int)
}

AEp05 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp05[,i] = expan(x, x0[6], N, t[i]) / (sintegral(x, expan(x, x0[6], N, t[i]))$int)
}

AEp06 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp06[,i] = expan(x, x0[7], N, t[i]) / (sintegral(x, expan(x, x0[7], N, t[i]))$int)
}

AEp07 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp07[,i] = expan(x, x0[8], N, t[i]) / (sintegral(x, expan(x, x0[8], N, t[i]))$int)
}

AEp08 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp08[,i] = expan(x, x0[9], N, t[i]) / (sintegral(x, expan(x, x0[9], N, t[i]))$int)
}

AEp09 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp09[,i] = expan(x, x0[10], N, t[i]) / (sintegral(x, expan(x, x0[10], N, t[i]))$int)
}

AEp1 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  AEp1[,i] = expan(x, x0[11], N, t[i]) / (sintegral(x, expan(x, x0[11], N, t[i]))$int)
}

hh3_expan = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh3_expan[i,1] =  1 - square_hellin(x, p0_adaptive[,i], AEp0[,i])
  hh3_expan[i,2] =  1 - square_hellin(x, p01_adaptive[,i], AEp01[,i])
  hh3_expan[i,3] =  1 - square_hellin(x, p02_adaptive[,i], AEp02[,i])
  hh3_expan[i,4] =  1 - square_hellin(x, p03_adaptive[,i], AEp03[,i])
  hh3_expan[i,5] =  1 - square_hellin(x, p04_adaptive[,i], AEp04[,i])
  hh3_expan[i,6] =  1 - square_hellin(x, p05_adaptive[,i], AEp05[,i])
  hh3_expan[i,7] =  1 - square_hellin(x, p06_adaptive[,i], AEp06[,i])
  hh3_expan[i,8] =  1 - square_hellin(x, p07_adaptive[,i], AEp07[,i])
  hh3_expan[i,9] =  1 - square_hellin(x, p08_adaptive[,i], AEp08[,i])
  hh3_expan[i,10] =  1 - square_hellin(x, p09_adaptive[,i], AEp09[,i])
  hh3_expan[i,11] =  1 - square_hellin(x, p1_adaptive[,i], AEp1[,i])
}
hh3_expan 

### hellinger gaussian approximation
gaussp0 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp0[,i] =  gauss(x, x0[1], N, t[i]) / (sintegral(x,  gauss(x, x0[1], N, t[i]))$int)
}

gaussp01 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp01[,i] =  gauss(x, x0[2], N, t[i]) / (sintegral(x,  gauss(x, x0[2], N, t[i]))$int)
}

gaussp02 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp02[,i] =  gauss(x, x0[3], N, t[i]) / (sintegral(x,  gauss(x, x0[3], N, t[i]))$int)
}

gaussp03 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp03[,i] =  gauss(x, x0[4], N, t[i]) / (sintegral(x,  gauss(x, x0[4], N, t[i]))$int)
}

gaussp04 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp04[,i] =  gauss(x, x0[5], N, t[i]) / (sintegral(x,  gauss(x, x0[5], N, t[i]))$int)
}

gaussp05 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp05[,i] =  gauss(x, x0[6], N, t[i]) / (sintegral(x,  gauss(x, x0[6], N, t[i]))$int)
}

gaussp06 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp06[,i] =  gauss(x, x0[7], N, t[i]) / (sintegral(x,  gauss(x, x0[7], N, t[i]))$int)
}

gaussp07 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp07[,i] =  gauss(x, x0[8], N, t[i]) / (sintegral(x,  gauss(x, x0[8], N, t[i]))$int)
}

gaussp08 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp08[,i] =  gauss(x, x0[9], N, t[i]) / (sintegral(x,  gauss(x, x0[9], N, t[i]))$int)
}

gaussp09 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp09[,i] =  gauss(x, x0[10], N, t[i]) / (sintegral(x,  gauss(x, x0[10], N, t[i]))$int)
}

gaussp1 = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  gaussp1[,i] =  gauss(x, x0[11], N, t[i]) / (sintegral(x,  gauss(x, x0[11], N, t[i]))$int)
}

hh3_gauss = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh3_gauss[i,1] = 1 - square_hellin(x, p0_adaptive[,i], gaussp0[,i])
  hh3_gauss[i,2] = 1 - square_hellin(x, p01_adaptive[,i], gaussp01[,i])
  hh3_gauss[i,3] = 1 - square_hellin(x, p02_adaptive[,i], gaussp02[,i])
  hh3_gauss[i,4] = 1 - square_hellin(x, p03_adaptive[,i], gaussp03[,i])
  hh3_gauss[i,5] = 1 - square_hellin(x, p04_adaptive[,i], gaussp04[,i])
  hh3_gauss[i,6] = 1 - square_hellin(x, p05_adaptive[,i], gaussp05[,i])
  hh3_gauss[i,7] = 1 - square_hellin(x, p06_adaptive[,i], gaussp06[,i])
  hh3_gauss[i,8] = 1 - square_hellin(x, p07_adaptive[,i], gaussp07[,i])
  hh3_gauss[i,9] = 1 - square_hellin(x, p08_adaptive[,i], gaussp08[,i])
  hh3_gauss[i,10] = 1 - square_hellin(x, p09_adaptive[,i], gaussp09[,i])
  hh3_gauss[i,11] = 1 - square_hellin(x, p1_adaptive[,i], gaussp1[,i])
}
hh3_gauss

###hellinger beta
hh3_beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh3_beta[i,1] = 1 - square_hellin(x, p0_adaptive[,i], p0_dbeta[,i])
  hh3_beta[i,2] = 1 - square_hellin(x, p01_adaptive[,i], p01_dbeta[,i])
  hh3_beta[i,3] = 1 - square_hellin(x, p02_adaptive[,i], p02_dbeta[,i])
  hh3_beta[i,4] = 1 - square_hellin(x, p03_adaptive[,i], p03_dbeta[,i])
  hh3_beta[i,5] = 1 - square_hellin(x, p04_adaptive[,i], p04_dbeta[,i])
  hh3_beta[i,6] = 1 - square_hellin(x, p05_adaptive[,i], p05_dbeta[,i])
  hh3_beta[i,7] = 1 - square_hellin(x, p06_adaptive[,i], p06_dbeta[,i])
  hh3_beta[i,8] = 1 - square_hellin(x, p07_adaptive[,i], p07_dbeta[,i])
  hh3_beta[i,9] = 1 - square_hellin(x, p08_adaptive[,i], p08_dbeta[,i])
  hh3_beta[i,10] = 1 - square_hellin(x, p09_adaptive[,i], p09_dbeta[,i])
  hh3_beta[i,11] = 1 - square_hellin(x, p1_adaptive[,i], p1_dbeta[,i])
}


hh3_beta


####hellinger normal con parámetros 
hh3_norm = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh3_norm[i,1]  = 1 - square_hellin(x, p0_adaptive[,i], p0_norm[,i])
  hh3_norm[i,2]  = 1 - square_hellin(x, p01_adaptive[,i], p01_norm[,i])
  hh3_norm[i,3]  = 1 - square_hellin(x, p02_adaptive[,i], p02_norm[,i])
  hh3_norm[i,4]  = 1 - square_hellin(x, p03_adaptive[,i], p03_norm[,i])
  hh3_norm[i,5]  = 1 - square_hellin(x, p04_adaptive[,i], p04_norm[,i])
  hh3_norm[i,6]  = 1 - square_hellin(x, p05_adaptive[,i], p05_norm[,i])
  hh3_norm[i,7]  = 1 - square_hellin(x, p06_adaptive[,i], p06_norm[,i])
  hh3_norm[i,8]  = 1 - square_hellin(x, p07_adaptive[,i], p07_norm[,i])
  hh3_norm[i,9]  = 1 - square_hellin(x, p08_adaptive[,i], p08_norm[,i])
  hh3_norm[i,10] = 1 - square_hellin(x, p09_adaptive[,i], p09_norm[,i])
  hh3_norm[i,11] = 1 - square_hellin(x, p1_adaptive[,i], p1_norm[,i])
}

hh3_norm


##### GRAFICOS

t_2N = rep(t, 11)
p0_df = c(rep(0,11), rep(0.1, 11), rep(0.2, 11), rep(0.3, 11), rep(0.4, 11), rep(0.5, 11), rep(0.6, 11), rep(0.7, 11), rep(0.8, 11), rep(0.9, 11), rep(1.0, 11))

hh3_expan_df = cbind.data.frame(t_2N, p0_df, rep("AE", 121), as.vector(hh3_expan))
colnames(hh3_expan_df) = c("1","2","3","4")

hh3_gauss_df = cbind.data.frame(t_2N, p0_df, rep("GaussA", 121), as.vector(hh3_gauss))
colnames(hh3_gauss_df) = c("1","2","3","4")

hh3_beta_df = cbind.data.frame(t_2N, p0_df, rep("Beta", 121), as.vector(hh3_beta))
colnames(hh3_beta_df) = c("1","2","3","4")

hh3_norm_df = cbind.data.frame(t_2N, p0_df, rep("Normal", 121), as.vector(hh3_norm))
colnames(hh3_norm_df) = c("1","2","3","4")

#juntando 
hh3_heat_df = rbind(hh3_expan_df, hh3_gauss_df, hh3_beta_df, hh3_norm_df)
colnames(hh3_heat_df) = c("t/2N", "p0", "distr", "Hel")

dim(hh3_heat_df)


###juntando 
hh3_heat_df <- within(hh3_heat_df, distr <- factor(distr, levels=c('AE', 'GaussA',
                                                                 'Beta','Normal')))
gg <- ggplot(hh3_heat_df, aes(x=p0, y=`t/2N`, fill=-log10(Hel)))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y="t", x = expression(x[0]))
gg <- gg + scale_fill_viridis_c(option = "magma", name = expression(-log[10](Hellinger)),  direction = -1)
gg <- gg + facet_wrap(vars(distr), ncol=2, dir = "h")
gg <- gg + theme(strip.background = element_blank())
#gg <- gg + labs(x=NULL, y=NULL, title="Titulo")
#gg <- gg + theme_tufte(base_family="Helvetica")
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(panel.border=element_blank())
#gg <- gg + theme(plot.title=element_text(hjust=0))
#gg <- gg + theme(strip.text=element_text(hjust=0))
gg <- gg + theme(panel.spacing.x=unit(0.5, "cm"))
gg <- gg + theme(panel.spacing.y=unit(0.5, "cm"))
gg <- gg + theme(legend.title=element_text(size=6))
gg <- gg + theme(legend.title.align=1)
gg <- gg + theme(legend.text=element_text(size=6))
gg <- gg + theme(legend.position="bottom")
gg <- gg + theme(legend.key.size=unit(0.6, "cm"))
gg <- gg + theme(legend.key.width=unit(1, "cm"))
gg <- gg + theme(legend.title = element_text(size=10))
gg <- gg + theme(legend.text = element_text(size=8))
gg <- gg + theme(panel.background = element_blank(),
                 plot.background = element_blank())
gg

#####Norma L2

##1. Asymtptotical Expansion
L2_AE = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  L2_AE[i,1] = sqrt(sum((p0_adaptive[,i] - AEp0[,i])^2))
  L2_AE[i,2] = sqrt(sum((p01_adaptive[,i] - AEp01[,i])^2))
  L2_AE[i,3] = sqrt(sum((p02_adaptive[,i] - AEp02[,i])^2))
  L2_AE[i,4] = sqrt(sum((p03_adaptive[,i] - AEp03[,i])^2))
  L2_AE[i,5] = sqrt(sum((p04_adaptive[,i] - AEp04[,i])^2))
  L2_AE[i,6] = sqrt(sum((p05_adaptive[,i] - AEp05[,i])^2))
  L2_AE[i,7] = sqrt(sum((p06_adaptive[,i] - AEp06[,i])^2))
  L2_AE[i,8] = sqrt(sum((p07_adaptive[,i] - AEp07[,i])^2))
  L2_AE[i,9] = sqrt(sum((p08_adaptive[,i] - AEp08[,i])^2))
  L2_AE[i,10] = sqrt(sum((p09_adaptive[,i] - AEp09[,i])^2))
  L2_AE[i,11] = sqrt(sum((p1_adaptive[,i] - AEp1[,i])^2))
}
L2_AE

##Gaussian approximation
L2_Gauss = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  L2_Gauss[i,1] = sqrt(sum((p0_adaptive[,i] - gaussp0[,i])^2))
  L2_Gauss[i,2] = sqrt(sum((p01_adaptive[,i] - gaussp01[,i])^2))
  L2_Gauss[i,3] = sqrt(sum((p02_adaptive[,i] - gaussp02[,i])^2))
  L2_Gauss[i,4] = sqrt(sum((p03_adaptive[,i] - gaussp03[,i])^2))
  L2_Gauss[i,5] = sqrt(sum((p04_adaptive[,i] - gaussp04[,i])^2))
  L2_Gauss[i,6] = sqrt(sum((p05_adaptive[,i] - gaussp05[,i])^2))
  L2_Gauss[i,7] = sqrt(sum((p06_adaptive[,i] - gaussp06[,i])^2))
  L2_Gauss[i,8] = sqrt(sum((p07_adaptive[,i] - gaussp07[,i])^2))
  L2_Gauss[i,9] = sqrt(sum((p08_adaptive[,i] - gaussp08[,i])^2))
  L2_Gauss[i,10]= sqrt(sum((p09_adaptive[,i] - gaussp09[,i])^2))
  L2_Gauss[i,11] = sqrt(sum((p1_adaptive[,i] - gaussp1[,i])^2))
}
L2_Gauss


##Beta con parámetros adecuados
L2_Beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  L2_Beta[i,1]  = sqrt(sum((p0_adaptive[,i] - p0_dbeta[,i])^2))
  L2_Beta[i,2]  = sqrt(sum((p01_adaptive[,i] - p01_dbeta[,i])^2))
  L2_Beta[i,3]  = sqrt(sum((p02_adaptive[,i] - p02_dbeta[,i])^2))
  L2_Beta[i,4]  = sqrt(sum((p03_adaptive[,i] - p03_dbeta[,i])^2))
  L2_Beta[i,5]  = sqrt(sum((p04_adaptive[,i] - p04_dbeta[,i])^2))
  L2_Beta[i,6]  = sqrt(sum((p05_adaptive[,i] - p05_dbeta[,i])^2))
  L2_Beta[i,7]  = sqrt(sum((p06_adaptive[,i] - p06_dbeta[,i])^2))
  L2_Beta[i,8]  = sqrt(sum((p07_adaptive[,i] - p07_dbeta[,i])^2))
  L2_Beta[i,9]  = sqrt(sum((p08_adaptive[,i] - p08_dbeta[,i])^2))
  L2_Beta[i,10] = sqrt(sum((p09_adaptive[,i] - p09_dbeta[,i])^2))
  L2_Beta[i,11] = sqrt(sum((p1_adaptive[,i] - p1_dbeta[,i])^2))
}
L2_Beta


##Normal con parámetros adecuados
L2_Norm = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  L2_Norm[i,1] = sqrt(sum((p0_adaptive[,i] - p0_norm[,i])^2))
  L2_Norm[i,2] = sqrt(sum((p01_adaptive[,i] - p01_norm[,i])^2))
  L2_Norm[i,3] = sqrt(sum((p02_adaptive[,i] - p02_norm[,i])^2))
  L2_Norm[i,4] = sqrt(sum((p03_adaptive[,i] - p03_norm[,i])^2))
  L2_Norm[i,5] = sqrt(sum((p04_adaptive[,i] - p04_norm[,i])^2))
  L2_Norm[i,6] = sqrt(sum((p05_adaptive[,i] - p05_norm[,i])^2))
  L2_Norm[i,7] = sqrt(sum((p06_adaptive[,i] - p06_norm[,i])^2))
  L2_Norm[i,8] = sqrt(sum((p07_adaptive[,i] - p07_norm[,i])^2))
  L2_Norm[i,9] = sqrt(sum((p08_adaptive[,i] - p08_norm[,i])^2))
  L2_Norm[i,10]= sqrt(sum((p09_adaptive[,i] - p09_norm[,i])^2))
  L2_Norm[i,11] = sqrt(sum((p1_adaptive[,i] - p1_norm[,i])^2))
}
L2_Norm


##### GRAFICOS

t_2N = rep(t, 11)
p0_df = c(rep(0, 11), rep(0.1, 11), rep(0.2, 11), rep(0.3, 11), rep(0.4, 11), rep(0.5, 11), rep(0.6, 11), rep(0.7, 11), rep(0.8, 11), rep(0.9, 11), rep(1, 11))

L2_expan_df = cbind.data.frame(t_2N, p0_df, rep("AE", 121), as.vector(L2_AE))
colnames(L2_expan_df) = c("1","2","3","4")

L2_gauss_df = cbind.data.frame(t_2N, p0_df, rep("GaussA", 121), as.vector(L2_Gauss))
colnames(L2_gauss_df) = c("1","2","3","4")

L2_beta_df = cbind.data.frame(t_2N, p0_df, rep("Beta", 121), as.vector(L2_Beta))
colnames(L2_beta_df) = c("1","2","3","4")

L2_norm_df = cbind.data.frame(t_2N, p0_df, rep("Normal", 121), as.vector(L2_Norm))
colnames(L2_norm_df) = c("1","2","3","4")

#juntando 
L2_heat_df = rbind(L2_expan_df, L2_gauss_df, L2_beta_df, L2_norm_df)
colnames(L2_heat_df) = c("t/2N", "p0", "distr", "L2")

dim(L2_heat_df)


###juntando 
##Definiendo el orden 
L2_heat_df <- within(L2_heat_df, distr <- factor(distr, levels=c('AE', 'GaussA',
                                                                 'Beta','Normal')))
##cambio L2
#gg <- ggplot(L2_heat_df, aes(x=p0, y=`t/2N`, fill=L2))

gg <- ggplot(L2_heat_df, aes(x=p0, y=`t/2N`, fill=log10(L2)))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y="t", x = expression(x[0]))
gg <- gg + scale_fill_viridis_c(option = "magma", name = expression(log[10](L^2)))
#gg <- gg + coord_equal()
gg <- gg + facet_wrap(vars(distr), ncol=2, dir = "h")
gg <- gg + theme(strip.background = element_blank())
#gg <- gg + labs(x=NULL, y=NULL, title="Titulo")
#gg <- gg + theme_tufte(base_family="Helvetica")
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(panel.border=element_blank())
#gg <- gg + theme(plot.title=element_text(hjust=0))
#gg <- gg + theme(strip.text=element_text(hjust=0))
gg <- gg + theme(panel.spacing.x=unit(0.5, "cm"))
gg <- gg + theme(panel.spacing.y=unit(0.5, "cm"))
gg <- gg + theme(legend.title=element_text(size=6))
gg <- gg + theme(legend.title.align=1)
gg <- gg + theme(legend.text=element_text(size=6))
gg <- gg + theme(legend.position="bottom")
gg <- gg + theme(legend.key.size=unit(0.6, "cm"))
gg <- gg + theme(legend.key.width=unit(1, "cm"))
gg <- gg + theme(legend.title = element_text(size=10))
gg <- gg + theme(legend.text = element_text(size=8))
gg <- gg + theme(panel.background = element_blank(),
                 plot.background = element_blank())
gg

#ver valores de t más grandes
##porcentajes
sum((hh3_expan < 0.2)*1)/121 
sum((hh3_gauss < 0.2)*1)/121 
sum((hh3_beta < 0.2)*1)/121 
sum((hh3_norm < 0.2)*1)/121 


sum((L2_AE<20)*1)/121
sum((L2_Gauss<20)*1)/121
sum((L2_Beta<20)*1)/121
sum((L2_Norm<20)*1)/121

##### mejora de las funciones AE y GaussA

#t=100
a1 = expan(x1, 0.1, 1000, t[1])
a1[is.nan(a1)] = max(!is.nan(a1))

a2 = expan(x1, 0.5, 1000, t[1])
a2[is.nan(a2)] = max(!is.nan(a2))

a3 = expan(x1, 0.9, 1000, t[1])
a3[is.nan(a3)] = max(!is.nan(a3))

#t=250
a4 = expan(x1, 0.1, 1000, t[6])
a4[is.nan(a4)] = max(!is.nan(a4))

a5 = expan(x1, 0.5, 1000, t[6])
a5[is.nan(a5)] = max(!is.nan(a5))

a6 = expan(x1, 0.9, 1000, t[6])
a6[is.nan(a6)] = max(!is.nan(a6))


#t=450
a7 = expan(x1, 0.1, 1000, t[10])
a7[is.nan(a7)] = max(!is.nan(a7))

a8 = expan(x1, 0.5, 1000, t[10])
a8[is.nan(a8)] = max(!is.nan(a8))

a9 = expan(x1, 0.9, 1000, t[10])
a9[is.nan(a9)] = max(!is.nan(a9))


###Histogramas
#m = matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3, byrow = T)
#layout(mat = m, heights = c(0.8,0.2))

#t=100
par(mfrow=c(1,3))
x1 = seq(0,1, length.out = 100)
hist(p0_inicial01[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, AEp01[,3], lwd = 2, col ="red")
lines(x1, gaussp01[,3], lwd = 2, col ="green")
lines(p01_bt100_adaptive[,2], p01_bt100_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial05[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a2, lwd = 2, col ="red")
lines(x1, gaussp05[,3], lwd = 2, col ="green")
lines(p05_bt100_adaptive[,2], p05_bt100_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial09[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a3, lwd = 2, col ="red")
lines(x1, gaussp09[,3], lwd = 2, col ="green")
lines(p09_bt100_adaptive[,2], p09_bt100_adaptive[,1], lwd = 2, col = 'darkblue')

#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#plot_colors = c('red', 'green', 'darkblue')
#legend(x = "top",inset = 0,
 #      legend = c("AE", "GaussA", "ADE"), 
  #     col=plot_colors, lwd=5, cex=1.2, horiz = TRUE)

#t=250
hist(p0_inicial01[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a4, lwd = 2, col ="red")
lines(x1, gaussp01[,6], lwd = 2, col ="green")
lines(p01_bt250_adaptive[,2], p01_bt250_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial05[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a5, lwd = 2, col ="red")
lines(x1, gaussp05[,6], lwd = 2, col ="green")
lines(p05_bt250_adaptive[,2], p05_bt250_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial09[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a6, lwd = 2, col ="red")
lines(x1, gaussp09[,6], lwd = 2, col ="green")
lines(p09_bt250_adaptive[,2], p09_bt250_adaptive[,1], lwd = 2, col = 'darkblue')


#t=450
hist(p0_inicial01[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a7, lwd = 2, col ="red")
lines(x1, gaussp01[,10], lwd = 2, col ="green")
lines(p01_bt450_adaptive[,2], p01_bt450_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial05[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, a8, lwd = 2, col ="red")
lines(x1, gaussp05[,10], lwd = 2, col ="green")
lines(p05_bt450_adaptive[,2], p05_bt450_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial09[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x1, AEp09[,10], lwd = 2, col ="red")
lines(x1, gaussp09[,10], lwd = 2, col ="green")
lines(p09_bt450_adaptive[,2], p09_bt450_adaptive[,1], lwd = 2, col = 'darkblue')
par(mfrow=c(1,1))


####distancia de Hellinger y norma L2 con log10

#1. distancias de Hellinger

hh3_expan_log = -log10(hh3_expan)
hh3_gauss_log = -log10(hh3_gauss)
hh3_beta_log = -log10(hh3_beta)
hh3_norm_log = -log10(hh3_norm)

##### GRAFICOS
t_2N = rep(t, 11)
p0_df = c(rep(0, 11), rep(0.1, 11), rep(0.2, 11), rep(0.3, 11), rep(0.4, 11), rep(0.5, 11), rep(0.6, 11), rep(0.7, 11), rep(0.8, 11), rep(0.9, 11), rep(1, 11))

log_expan_df = cbind.data.frame(t_2N, p0_df, rep("AE", 121), as.vector(hh3_expan_log))
colnames(log_expan_df) = c("1","2","3","4")

log_gauss_df = cbind.data.frame(t_2N, p0_df, rep("GaussA", 121), as.vector(hh3_gauss_log))
colnames(log_gauss_df) = c("1","2","3","4")

log_beta_df = cbind.data.frame(t_2N, p0_df, rep("Beta", 121), as.vector(hh3_beta_log))
colnames(log_beta_df) = c("1","2","3","4")

log_norm_df = cbind.data.frame(t_2N, p0_df, rep("Normal", 121), as.vector(hh3_norm_log))
colnames(log_norm_df) = c("1","2","3","4")

#juntando 
log_heat_df = rbind(log_expan_df, log_gauss_df, log_beta_df, log_norm_df)
colnames(log_heat_df) = c("t/2N", "p0", "distr", "L2")

dim(log_heat_df)


###juntando 
##Definiendo el orden 
log_heat_df <- within(log_heat_df, distr <- factor(distr, levels=c('AE', 'GaussA',
                                                                 'Beta','Normal')))

gg <- ggplot(log_heat_df, aes(x=p0, y=`t/2N`, fill=L2))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y="t", x = expression(x[0]))
gg <- gg + scale_fill_viridis_c(option = "magma", name = "")
#gg <- gg + coord_equal()
gg <- gg + facet_wrap(vars(distr), ncol=2, dir = "h")
gg <- gg + theme(strip.background = element_blank())
#gg <- gg + labs(x=NULL, y=NULL, title="Titulo")
#gg <- gg + theme_tufte(base_family="Helvetica")
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_text(size=5))
gg <- gg + theme(panel.border=element_blank())
#gg <- gg + theme(plot.title=element_text(hjust=0))
#gg <- gg + theme(strip.text=element_text(hjust=0))
gg <- gg + theme(panel.spacing.x=unit(0.5, "cm"))
gg <- gg + theme(panel.spacing.y=unit(0.5, "cm"))
gg <- gg + theme(legend.title=element_text(size=6))
gg <- gg + theme(legend.title.align=1)
gg <- gg + theme(legend.text=element_text(size=6))
gg <- gg + theme(legend.position="bottom")
gg <- gg + theme(legend.key.size=unit(0.3, "cm"))
gg <- gg + theme(legend.key.width=unit(1, "cm"))
gg <- gg + theme(panel.background = element_blank(),
                 plot.background = element_blank())
gg

##2. norma L2

log_L2AE = log10(L2_AE)
log_L2Gauss = log10(L2_Gauss)
log_L2Beta = log10(L2_Beta)
log_L2Norm = log10(L2_Norm)

##### GRAFICOS
t_2N = rep(t, 11)
p0_df = c(rep(0, 11), rep(0.1, 11), rep(0.2, 11), rep(0.3, 11), rep(0.4, 11), rep(0.5, 11), rep(0.6, 11), rep(0.7, 11), rep(0.8, 11), rep(0.9, 11), rep(1, 11))

logL2_expan_df = cbind.data.frame(t_2N, p0_df, rep("AE", 121), as.vector(log_L2AE))
colnames(logL2_expan_df) = c("1","2","3","4")

logL2_gauss_df = cbind.data.frame(t_2N, p0_df, rep("GaussA", 121), as.vector(log_L2Gauss))
colnames(logL2_gauss_df) = c("1","2","3","4")

logL2_beta_df = cbind.data.frame(t_2N, p0_df, rep("Beta", 121), as.vector(log_L2Beta))
colnames(logL2_beta_df) = c("1","2","3","4")

logL2_norm_df = cbind.data.frame(t_2N, p0_df, rep("Normal", 121), as.vector(log_L2Norm))
colnames(logL2_norm_df) = c("1","2","3","4")

#juntando 
logL2_heat_df = rbind(logL2_expan_df, logL2_gauss_df, logL2_beta_df, logL2_norm_df)
colnames(logL2_heat_df) = c("t/2N", "p0", "distr", "L2")

dim(logL2_heat_df)


###juntando 
##Definiendo el orden 
logL2_heat_df <- within(logL2_heat_df, distr <- factor(distr, levels=c('AE', 'GaussA',
                                                                   'Beta','Normal')))

gg <- ggplot(logL2_heat_df, aes(x=p0, y=`t/2N`, fill=L2))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y="t", x = expression(x[0]))
gg <- gg + scale_fill_viridis_c(option = "magma", name = "")
#gg <- gg + coord_equal()
gg <- gg + facet_wrap(vars(distr), ncol=2, dir = "h")
#gg <- gg + labs(x=NULL, y=NULL, title="Titulo")
#gg <- gg + theme_tufte(base_family="Helvetica")
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_text(size=5))
gg <- gg + theme(panel.border=element_blank())
#gg <- gg + theme(plot.title=element_text(hjust=0))
#gg <- gg + theme(strip.text=element_text(hjust=0))
gg <- gg + theme(panel.spacing.x=unit(0.5, "cm"))
gg <- gg + theme(panel.spacing.y=unit(0.5, "cm"))
gg <- gg + theme(legend.title=element_text(size=6))
gg <- gg + theme(legend.title.align=1)
gg <- gg + theme(legend.text=element_text(size=6))
gg <- gg + theme(legend.position="bottom")
gg <- gg + theme(legend.key.size=unit(0.3, "cm"))
gg <- gg + theme(legend.key.width=unit(1, "cm"))
gg <- gg + theme(panel.background = element_blank(),
                 plot.background = element_blank())
gg
