library(sfsmisc)


#distancia de Hellinger

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
t = t0/(2*N) 
t


## LA formula (asymptotical expansion)
expan = function(x, x0, N, t){
  a_expan = (1/sqrt(2*pi*t)) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1/(2*t)) * (asin(1-2*x0) - asin(1-2*x))^(2))
  #a_expan = a_expan[!is.infinite(a_expan)] #quita el primero y el último que son Inf
  #a_expan[1] = a_expan[N] = 0
  
  ##MC to approximate integral
  samp = runif(1e6, min = 0.01, max = 0.99)
  
  f_sample = (1/sqrt(2*pi*t)) * (x0 * (1-x0))^(1/4) / (samp*(1-samp))^(3/4) * exp(-(1/(2*t)) * (asin(1-2*x0) - asin(1-2*samp))^(2))
  
  #return(mean(f_sample))
  return(a_expan/mean(f_sample))
}



plot(seq(0.01,0.99, length.out = 100), expan(x, 0.01, 100, 100), type = "l")

## LA formula (gaussian approximation)
gauss = function(x, x0, N, t){
  gauss_app = (1/sqrt(2*pi*t)) * (1/ sqrt((x0*(1-x0)))) * exp(-1/(2*t) * ((x-x0)^2) / (x0*(1-x0)) )
  
  ##MC to approximate integral
  samp = runif(1e6, min = 0.01, max = 0.99)
  
  f_sample = (1/sqrt(2*pi*t)) * (1/ sqrt((x0*(1-x0)))) * exp(-1/(2*t) * ((samp-x0)^2) / (x0*(1-x0)) )
  
  #return(mean(f_sample))
  return(gauss_app/mean(f_sample))
}


plot(seq(0.01, 0.99, length.out = 100), gauss(x, 0.2, 100, 100), type = 'l')


### cambio en la función para la distancia de Hellinger
hellinger2 = function (x, y, lower = -Inf, upper = Inf, method = 1,...) 
{
  fx <- function(x) (x)
  fy <- function(y) (y)
  if (method == 1) {
    g <- function(z) (fx(z)^0.5 - fy(z)^0.5)^2
    #increasing tolerance / subdivisions
    h2 <- stats::integrate(g, lower, upper, rel.tol=.Machine$double.eps^.05)$value/2
  }
  else if (method == 2) {
    g <- function(z) (fx(z) * fy(z))^0.5
    h2 <- 1 - stats::integrate(g, lower, upper)$value
  }
  else {
    stop("incorrect 'method' argument", call. = FALSE)
  }
  return(sqrt(h2))
}


##Resultados del estimador adaptativo de B-K
p0_adaptive = cbind(p0_bt1_adaptive[,1] / (sintegral(p0_bt1_adaptive[,2], p0_bt1_adaptive[,1])$int), 
                    p0_bt50_adaptive[,1] / (sintegral(p0_bt50_adaptive[,2], p0_bt50_adaptive[,1])$int),
                    p0_bt100_adaptive[,1] / (sintegral(p0_bt100_adaptive[,2], p0_bt100_adaptive[,1])$int),
                    p0_bt150_adaptive[,1] / (sintegral(p0_bt150_adaptive[,2], p0_bt150_adaptive[,1])$int),
                    p0_bt200_adaptive[,1] / (sintegral(p0_bt200_adaptive[,2], p0_bt200_adaptive[,1])$int),
                    p0_bt250_adaptive[,1] / (sintegral(p0_bt250_adaptive[,2], p0_bt250_adaptive[,1])$int),
                    p0_bt300_adaptive[,1] / (sintegral(p0_bt300_adaptive[,2], p0_bt300_adaptive[,1])$int),
                    p0_bt350_adaptive[,1] / (sintegral(p0_bt350_adaptive[,2], p0_bt350_adaptive[,1])$int),
                    p0_bt400_adaptive[,1] / (sintegral(p0_bt400_adaptive[,2], p0_bt400_adaptive[,1])$int),
                    p0_bt450_adaptive[,1] / (sintegral(p0_bt450_adaptive[,2], p0_bt450_adaptive[,1])$int),
                    p0_bt500_adaptive[,1] / (sintegral(p0_bt500_adaptive[,2], p0_bt500_adaptive[,1])$int))


p01_adaptive = cbind(p01_bt1_adaptive[,1] / (sintegral(p01_bt1_adaptive[,2], p01_bt1_adaptive[,1])$int), 
                     p01_bt50_adaptive[,1] / (sintegral(p01_bt50_adaptive[,2], p01_bt50_adaptive[,1])$int),
                     p01_bt100_adaptive[,1] / (sintegral(p01_bt100_adaptive[,2], p01_bt100_adaptive[,1])$int),
                     p01_bt150_adaptive[,1] / (sintegral(p01_bt150_adaptive[,2], p01_bt150_adaptive[,1])$int),
                     p01_bt200_adaptive[,1] / (sintegral(p01_bt200_adaptive[,2], p01_bt200_adaptive[,1])$int),
                     p01_bt250_adaptive[,1] / (sintegral(p01_bt250_adaptive[,2], p01_bt250_adaptive[,1])$int),
                     p01_bt300_adaptive[,1] / (sintegral(p01_bt300_adaptive[,2], p01_bt300_adaptive[,1])$int),
                     p01_bt350_adaptive[,1] / (sintegral(p01_bt350_adaptive[,2], p01_bt350_adaptive[,1])$int),
                     p01_bt400_adaptive[,1] / (sintegral(p01_bt400_adaptive[,2], p01_bt400_adaptive[,1])$int),
                     p01_bt450_adaptive[,1] / (sintegral(p01_bt450_adaptive[,2], p01_bt450_adaptive[,1])$int),
                     p01_bt500_adaptive[,1] / (sintegral(p01_bt500_adaptive[,2], p01_bt500_adaptive[,1])$int))

p02_adaptive = cbind(p02_bt1_adaptive[,1] / (sintegral(p02_bt1_adaptive[,2], p02_bt1_adaptive[,1])$int), 
                     p02_bt50_adaptive[,1] / (sintegral(p02_bt50_adaptive[,2], p02_bt50_adaptive[,1])$int),
                     p02_bt100_adaptive[,1] / (sintegral(p02_bt100_adaptive[,2], p02_bt100_adaptive[,1])$int),
                     p02_bt150_adaptive[,1] / (sintegral(p02_bt150_adaptive[,2], p02_bt150_adaptive[,1])$int),
                     p02_bt200_adaptive[,1] / (sintegral(p02_bt200_adaptive[,2], p02_bt200_adaptive[,1])$int),
                     p02_bt250_adaptive[,1] / (sintegral(p02_bt250_adaptive[,2], p02_bt250_adaptive[,1])$int),
                     p02_bt300_adaptive[,1] / (sintegral(p02_bt300_adaptive[,2], p02_bt300_adaptive[,1])$int),
                     p02_bt350_adaptive[,1] / (sintegral(p02_bt350_adaptive[,2], p02_bt350_adaptive[,1])$int),
                     p02_bt400_adaptive[,1] / (sintegral(p02_bt400_adaptive[,2], p02_bt400_adaptive[,1])$int),
                     p02_bt450_adaptive[,1] / (sintegral(p02_bt450_adaptive[,2], p02_bt450_adaptive[,1])$int),
                     p02_bt500_adaptive[,1] / (sintegral(p02_bt500_adaptive[,2], p02_bt500_adaptive[,1])$int))

p03_adaptive = cbind(p03_bt1_adaptive[,1] / (sintegral(p03_bt1_adaptive[,2], p03_bt1_adaptive[,1])$int), 
                     p03_bt50_adaptive[,1] / (sintegral(p03_bt50_adaptive[,2], p03_bt50_adaptive[,1])$int),
                     p03_bt100_adaptive[,1] / (sintegral(p03_bt100_adaptive[,2], p03_bt100_adaptive[,1])$int),
                     p03_bt150_adaptive[,1] / (sintegral(p03_bt150_adaptive[,2], p03_bt150_adaptive[,1])$int),
                     p03_bt200_adaptive[,1] / (sintegral(p03_bt200_adaptive[,2], p03_bt200_adaptive[,1])$int),
                     p03_bt250_adaptive[,1] / (sintegral(p03_bt250_adaptive[,2], p03_bt250_adaptive[,1])$int),
                     p03_bt300_adaptive[,1] / (sintegral(p03_bt300_adaptive[,2], p03_bt300_adaptive[,1])$int),
                     p03_bt350_adaptive[,1] / (sintegral(p03_bt350_adaptive[,2], p03_bt350_adaptive[,1])$int),
                     p03_bt400_adaptive[,1] / (sintegral(p03_bt400_adaptive[,2], p03_bt400_adaptive[,1])$int),
                     p03_bt450_adaptive[,1] / (sintegral(p03_bt450_adaptive[,2], p03_bt450_adaptive[,1])$int),
                     p03_bt500_adaptive[,1] / (sintegral(p03_bt500_adaptive[,2], p03_bt500_adaptive[,1])$int))

p04_adaptive = cbind(p04_bt1_adaptive[,1] / (sintegral(p04_bt1_adaptive[,2], p04_bt1_adaptive[,1])$int), 
                     p04_bt50_adaptive[,1] / (sintegral(p04_bt50_adaptive[,2], p04_bt50_adaptive[,1])$int),
                     p04_bt100_adaptive[,1] / (sintegral(p04_bt100_adaptive[,2], p04_bt100_adaptive[,1])$int),
                     p04_bt150_adaptive[,1] / (sintegral(p04_bt150_adaptive[,2], p04_bt150_adaptive[,1])$int),
                     p04_bt200_adaptive[,1] / (sintegral(p04_bt200_adaptive[,2], p04_bt200_adaptive[,1])$int),
                     p04_bt250_adaptive[,1] / (sintegral(p04_bt250_adaptive[,2], p04_bt250_adaptive[,1])$int),
                     p04_bt300_adaptive[,1] / (sintegral(p04_bt300_adaptive[,2], p04_bt300_adaptive[,1])$int),
                     p04_bt350_adaptive[,1] / (sintegral(p04_bt350_adaptive[,2], p04_bt350_adaptive[,1])$int),
                     p04_bt400_adaptive[,1] / (sintegral(p04_bt400_adaptive[,2], p04_bt400_adaptive[,1])$int),
                     p04_bt450_adaptive[,1] / (sintegral(p04_bt450_adaptive[,2], p04_bt450_adaptive[,1])$int),
                     p04_bt500_adaptive[,1] / (sintegral(p04_bt500_adaptive[,2], p04_bt500_adaptive[,1])$int))

p05_adaptive = cbind(p05_bt1_adaptive[,1] / (sintegral(p05_bt1_adaptive[,2], p05_bt1_adaptive[,1])$int), 
                     p05_bt50_adaptive[,1] / (sintegral(p05_bt50_adaptive[,2], p05_bt50_adaptive[,1])$int),
                     p05_bt100_adaptive[,1] / (sintegral(p05_bt100_adaptive[,2], p05_bt100_adaptive[,1])$int),
                     p05_bt150_adaptive[,1] / (sintegral(p05_bt150_adaptive[,2], p05_bt150_adaptive[,1])$int),
                     p05_bt200_adaptive[,1] / (sintegral(p05_bt200_adaptive[,2], p05_bt200_adaptive[,1])$int),
                     p05_bt250_adaptive[,1] / (sintegral(p05_bt250_adaptive[,2], p05_bt250_adaptive[,1])$int),
                     p05_bt300_adaptive[,1] / (sintegral(p05_bt300_adaptive[,2], p05_bt300_adaptive[,1])$int),
                     p05_bt350_adaptive[,1] / (sintegral(p05_bt350_adaptive[,2], p05_bt350_adaptive[,1])$int),
                     p05_bt400_adaptive[,1] / (sintegral(p05_bt400_adaptive[,2], p05_bt400_adaptive[,1])$int),
                     p05_bt450_adaptive[,1] / (sintegral(p05_bt450_adaptive[,2], p05_bt450_adaptive[,1])$int),
                     p05_bt500_adaptive[,1] / (sintegral(p05_bt500_adaptive[,2], p05_bt500_adaptive[,1])$int))

p06_adaptive = cbind(p06_bt1_adaptive[,1] / (sintegral(p06_bt1_adaptive[,2], p06_bt1_adaptive[,1])$int), 
                     p06_bt50_adaptive[,1] / (sintegral(p06_bt50_adaptive[,2], p06_bt50_adaptive[,1])$int),
                     p06_bt100_adaptive[,1] / (sintegral(p06_bt100_adaptive[,2], p06_bt100_adaptive[,1])$int),
                     p06_bt150_adaptive[,1] / (sintegral(p06_bt150_adaptive[,2], p06_bt150_adaptive[,1])$int),
                     p06_bt200_adaptive[,1] / (sintegral(p06_bt200_adaptive[,2], p06_bt200_adaptive[,1])$int),
                     p06_bt250_adaptive[,1] / (sintegral(p06_bt250_adaptive[,2], p06_bt250_adaptive[,1])$int),
                     p06_bt300_adaptive[,1] / (sintegral(p06_bt300_adaptive[,2], p06_bt300_adaptive[,1])$int),
                     p06_bt350_adaptive[,1] / (sintegral(p06_bt350_adaptive[,2], p06_bt350_adaptive[,1])$int),
                     p06_bt400_adaptive[,1] / (sintegral(p06_bt400_adaptive[,2], p06_bt400_adaptive[,1])$int),
                     p06_bt450_adaptive[,1] / (sintegral(p06_bt450_adaptive[,2], p06_bt450_adaptive[,1])$int),
                     p06_bt500_adaptive[,1] / (sintegral(p06_bt500_adaptive[,2], p06_bt500_adaptive[,1])$int))

p07_adaptive = cbind(p07_bt1_adaptive[,1] / (sintegral(p07_bt1_adaptive[,2], p07_bt1_adaptive[,1])$int), 
                     p07_bt50_adaptive[,1] / (sintegral(p07_bt50_adaptive[,2], p07_bt50_adaptive[,1])$int),
                     p07_bt100_adaptive[,1] / (sintegral(p07_bt100_adaptive[,2], p07_bt100_adaptive[,1])$int),
                     p07_bt150_adaptive[,1] / (sintegral(p07_bt150_adaptive[,2], p07_bt150_adaptive[,1])$int),
                     p07_bt200_adaptive[,1] / (sintegral(p07_bt200_adaptive[,2], p07_bt200_adaptive[,1])$int),
                     p07_bt250_adaptive[,1] / (sintegral(p07_bt250_adaptive[,2], p07_bt250_adaptive[,1])$int),
                     p07_bt300_adaptive[,1] / (sintegral(p07_bt300_adaptive[,2], p07_bt300_adaptive[,1])$int),
                     p07_bt350_adaptive[,1] / (sintegral(p07_bt350_adaptive[,2], p07_bt350_adaptive[,1])$int),
                     p07_bt400_adaptive[,1] / (sintegral(p07_bt400_adaptive[,2], p07_bt400_adaptive[,1])$int),
                     p07_bt450_adaptive[,1] / (sintegral(p07_bt450_adaptive[,2], p07_bt450_adaptive[,1])$int),
                     p07_bt500_adaptive[,1] / (sintegral(p07_bt500_adaptive[,2], p07_bt500_adaptive[,1])$int))

p08_adaptive = cbind(p08_bt1_adaptive[,1] / (sintegral(p08_bt1_adaptive[,2], p08_bt1_adaptive[,1])$int), 
                     p08_bt50_adaptive[,1] / (sintegral(p08_bt50_adaptive[,2], p08_bt50_adaptive[,1])$int),
                     p08_bt100_adaptive[,1] / (sintegral(p08_bt100_adaptive[,2], p08_bt100_adaptive[,1])$int),
                     p08_bt150_adaptive[,1] / (sintegral(p08_bt150_adaptive[,2], p08_bt150_adaptive[,1])$int),
                     p08_bt200_adaptive[,1] / (sintegral(p08_bt200_adaptive[,2], p08_bt200_adaptive[,1])$int),
                     p08_bt250_adaptive[,1] / (sintegral(p08_bt250_adaptive[,2], p08_bt250_adaptive[,1])$int),
                     p08_bt300_adaptive[,1] / (sintegral(p08_bt300_adaptive[,2], p08_bt300_adaptive[,1])$int),
                     p08_bt350_adaptive[,1] / (sintegral(p08_bt350_adaptive[,2], p08_bt350_adaptive[,1])$int),
                     p08_bt400_adaptive[,1] / (sintegral(p08_bt400_adaptive[,2], p08_bt400_adaptive[,1])$int),
                     p08_bt450_adaptive[,1] / (sintegral(p08_bt450_adaptive[,2], p08_bt450_adaptive[,1])$int),
                     p08_bt500_adaptive[,1] / (sintegral(p08_bt500_adaptive[,2], p08_bt500_adaptive[,1])$int))

p09_adaptive = cbind(p09_bt1_adaptive[,1] / (sintegral(p09_bt1_adaptive[,2], p09_bt1_adaptive[,1])$int), 
                     p09_bt50_adaptive[,1] / (sintegral(p09_bt50_adaptive[,2], p09_bt50_adaptive[,1])$int),
                     p09_bt100_adaptive[,1] / (sintegral(p09_bt100_adaptive[,2], p09_bt100_adaptive[,1])$int),
                     p09_bt150_adaptive[,1] / (sintegral(p09_bt150_adaptive[,2], p09_bt150_adaptive[,1])$int),
                     p09_bt200_adaptive[,1] / (sintegral(p09_bt200_adaptive[,2], p09_bt200_adaptive[,1])$int),
                     p09_bt250_adaptive[,1] / (sintegral(p09_bt250_adaptive[,2], p09_bt250_adaptive[,1])$int),
                     p09_bt300_adaptive[,1] / (sintegral(p09_bt300_adaptive[,2], p09_bt300_adaptive[,1])$int),
                     p09_bt350_adaptive[,1] / (sintegral(p09_bt350_adaptive[,2], p09_bt350_adaptive[,1])$int),
                     p09_bt400_adaptive[,1] / (sintegral(p09_bt400_adaptive[,2], p09_bt400_adaptive[,1])$int),
                     p09_bt450_adaptive[,1] / (sintegral(p09_bt450_adaptive[,2], p09_bt450_adaptive[,1])$int),
                     p09_bt500_adaptive[,1] / (sintegral(p09_bt500_adaptive[,2], p09_bt500_adaptive[,1])$int))

p1_adaptive = cbind(p1_bt1_adaptive[,1] / (sintegral(p1_bt1_adaptive[,2], p1_bt1_adaptive[,1])$int), 
                    p1_bt50_adaptive[,1] / (sintegral(p1_bt50_adaptive[,2], p1_bt50_adaptive[,1])$int),
                    p1_bt100_adaptive[,1] / (sintegral(p1_bt100_adaptive[,2], p1_bt100_adaptive[,1])$int),
                    p1_bt150_adaptive[,1] / (sintegral(p1_bt150_adaptive[,2], p1_bt150_adaptive[,1])$int),
                    p1_bt200_adaptive[,1] / (sintegral(p1_bt200_adaptive[,2], p1_bt200_adaptive[,1])$int),
                    p1_bt250_adaptive[,1] / (sintegral(p1_bt250_adaptive[,2], p1_bt250_adaptive[,1])$int),
                    p1_bt300_adaptive[,1] / (sintegral(p1_bt300_adaptive[,2], p1_bt300_adaptive[,1])$int),
                    p1_bt350_adaptive[,1] / (sintegral(p1_bt350_adaptive[,2], p1_bt350_adaptive[,1])$int),
                    p1_bt400_adaptive[,1] / (sintegral(p1_bt400_adaptive[,2], p1_bt400_adaptive[,1])$int),
                    p1_bt450_adaptive[,1] / (sintegral(p1_bt450_adaptive[,2], p1_bt450_adaptive[,1])$int),
                    p1_bt500_adaptive[,1] / (sintegral(p1_bt500_adaptive[,2], p1_bt500_adaptive[,1])$int))


##1 asymtotical expansion

hh_expan = matrix(ncol = 11, nrow = 11)
for (i in 1:11) {
  hh_expan[i,1] = hellinger(p0_adaptive[,i], expan(x, x0[1], N, t0[i]), 0, 1)
  hh_expan[i,2] = hellinger(p01_adaptive[,i], expan(x, x0[2], N, t0[i]), 0, 1)
  hh_expan[i,3] = hellinger(p02_adaptive[,i], expan(x, x0[3], N, t0[i]), 0, 1)
  hh_expan[i,4] = hellinger(p03_adaptive[,i], expan(x, x0[4], N, t0[i]), 0, 1)
  hh_expan[i,5] = hellinger(p04_adaptive[,i], expan(x, x0[5], N, t0[i]), 0, 1)
  hh_expan[i,6] = hellinger(p05_adaptive[,i], expan(x, x0[6], N, t0[i]), 0, 1)
  hh_expan[i,7] = hellinger(p06_adaptive[,i], expan(x, x0[7], N, t0[i]), 0, 1)
  hh_expan[i,8] = hellinger(p07_adaptive[,i], expan(x, x0[8], N, t0[i]), 0, 1)
  hh_expan[i,9] = hellinger(p08_adaptive[,i], expan(x, x0[9], N, t0[i]), 0, 1)
  hh_expan[i,10] = hellinger(p09_adaptive[,i], expan(x, x0[10], N, t0[i]), 0, 1)
  hh_expan[i,11] = hellinger(p1_adaptive[,i], expan(x, x0[11], N, t0[i]), 0, 1)
}

colnames(hh_expan) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_expan) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_expan

######2 beta con parametros adecuados
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

##parametros de forma en función a la media y la varianza
sh1 = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  for (j in 1:11){
    sh1[i,j] = (((media[i,j] * (1 - media[i,j])) / variance[i,j]) - 1) * media[i,j]
  }
}

colnames(sh1) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(sh1) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
sh1

sh2 = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  for (j in 1:11){
    sh2[i,j] = (((media[i,j] * (1 - media[i,j])) / variance[i,j]) - 1) * (1 - media[i,j])
  }
}

colnames(sh2) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(sh2) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
sh2

##Para asegurarse que integre 1 los voy a almacenar en una matriz
p0_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p0_dbeta[,i] = dbeta(x, shape1 = sh1[i,1], shape2 = sh2[i,1]) / (sintegral(x, dbeta(x, shape1 = sh1[i,1], shape2 = sh2[i,1]))$int)
}

p01_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p01_dbeta[,i] = dbeta(x, shape1 = sh1[i,2], shape2 = sh2[i,2]) / (sintegral(x, dbeta(x, shape1 = sh1[i,2], shape2 = sh2[i,2]))$int)
}

p02_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p02_dbeta[,i] = dbeta(x, shape1 = sh1[i,3], shape2 = sh2[i,3]) / (sintegral(x, dbeta(x, shape1 = sh1[i,3], shape2 = sh2[i,3]))$int)
}

p03_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p03_dbeta[,i] = dbeta(x, shape1 = sh1[i,4], shape2 = sh2[i,4]) / (sintegral(x, dbeta(x, shape1 = sh1[i,4], shape2 = sh2[i,4]))$int)
}

p04_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p04_dbeta[,i] = dbeta(x, shape1 = sh1[i,5], shape2 = sh2[i,5]) / (sintegral(x, dbeta(x, shape1 = sh1[i,5], shape2 = sh2[i,5]))$int)
}

p05_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p05_dbeta[,i] = dbeta(x, shape1 = sh1[i,6], shape2 = sh2[i,6]) / (sintegral(x, dbeta(x, shape1 = sh1[i,6], shape2 = sh2[i,6]))$int)
}

p06_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p06_dbeta[,i] = dbeta(x, shape1 = sh1[i,7], shape2 = sh2[i,7]) / (sintegral(x, dbeta(x, shape1 = sh1[i,7], shape2 = sh2[i,7]))$int)
}

p07_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p07_dbeta[,i] = dbeta(x, shape1 = sh1[i,8], shape2 = sh2[i,8]) / (sintegral(x, dbeta(x, shape1 = sh1[i,8], shape2 = sh2[i,8]))$int)
}

p08_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p08_dbeta[,i] = dbeta(x, shape1 = sh1[i,9], shape2 = sh2[i,9]) / (sintegral(x, dbeta(x, shape1 = sh1[i,9], shape2 = sh2[i,9]))$int)
}

p09_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p09_dbeta[,i] = dbeta(x, shape1 = sh1[i,10], shape2 = sh2[i,10]) / (sintegral(x, dbeta(x, shape1 = sh1[i,10], shape2 = sh2[i,10]))$int)
}

p1_dbeta = matrix(NA, ncol = 11, nrow = 100)
for (i in 1:11) {
  p1_dbeta[,i] = dbeta(x, shape1 = sh1[i,11], shape2 = sh2[i,11]) / (sintegral(x, dbeta(x, shape1 = sh1[i,11], shape2 = sh2[i,11]))$int)
}


hh_beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh_beta[i,1] = hellinger(p0_adaptive[,i], p0_dbeta[,i], 0, 1)
  hh_beta[i,2] = hellinger(p01_adaptive[,i], p01_dbeta[,i], 0, 1)
  hh_beta[i,3] = hellinger(p02_adaptive[,i], p02_dbeta[,i], 0, 1)
  hh_beta[i,4] = hellinger(p03_adaptive[,i], p03_dbeta[,i], 0, 1)
  hh_beta[i,5] = hellinger(p04_adaptive[,i], p04_dbeta[,i], 0, 1)
  hh_beta[i,6] = hellinger(p05_adaptive[,i], p05_dbeta[,i], 0, 1)
  hh_beta[i,7] = hellinger(p06_adaptive[,i], p06_dbeta[,i], 0, 1)
  hh_beta[i,8] = hellinger(p07_adaptive[,i], p07_dbeta[,i], 0, 1)
  hh_beta[i,9] = hellinger(p08_adaptive[,i], p08_dbeta[,i], 0, 1)
  hh_beta[i,10] = hellinger(p09_adaptive[,i],p09_dbeta[,i], 0, 1)
  hh_beta[i,11] = hellinger(p1_adaptive[,i], p1_dbeta[,i], 0, 1)
}

colnames(hh_beta) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_beta) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_beta



###beta hellinger 2.1 
hh2_beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh2_beta[i,1] = hellinger2(p0_adaptive[,i], dbeta(x, shape1 = sh1[i,1], shape2 = sh2[i,1]), 0, 1)
  hh2_beta[i,2] = hellinger2(p01_adaptive[,i], dbeta(x, shape1 = sh1[i,2], shape2 = sh2[i,2]), 0, 1)
  hh2_beta[i,3] = hellinger2(p02_adaptive[,i], dbeta(x, shape1 = sh1[i,3], shape2 = sh2[i,3]), 0, 1)
  hh2_beta[i,4] = hellinger2(p03_adaptive[,i], dbeta(x, shape1 = sh1[i,4], shape2 = sh2[i,4]), 0, 1)
  hh2_beta[i,5] = hellinger2(p04_adaptive[,i], dbeta(x, shape1 = sh1[i,5], shape2 = sh2[i,5]), 0, 1)
  hh2_beta[i,6] = hellinger2(p05_adaptive[,i], dbeta(x, shape1 = sh1[i,6], shape2 = sh2[i,6]), 0, 1)
  hh2_beta[i,7] = hellinger2(p06_adaptive[,i], dbeta(x, shape1 = sh1[i,7], shape2 = sh2[i,7]), 0, 1)
  hh2_beta[i,8] = hellinger2(p07_adaptive[,i], dbeta(x, shape1 = sh1[i,8], shape2 = sh2[i,8]), 0, 1)
  hh2_beta[i,9] = hellinger2(p08_adaptive[,i], dbeta(x, shape1 = sh1[i,9], shape2 = sh2[i,9]), 0, 1)
  hh2_beta[i,10] = hellinger2(p09_adaptive[,i], dbeta(x, shape1 = sh1[i,10], shape2 = sh2[i,10]), 0, 1)
  hh2_beta[i,11] = hellinger2(p1_adaptive[,i], dbeta(x, shape1 = sh1[i,11], shape2 = sh2[i,11]), 0, 1)
}

colnames(hh2_beta) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh2_beta) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh2_beta


##3 normal con parametros adecuados
p0_norm = cbind(dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[1]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[1]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[2]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[2]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[3]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[3]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[4]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[4]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[5]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[5]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[6]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[6]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[7]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[7]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[8]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[8]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[9]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[9]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[10]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[10]/(2*N)))))$int,
                dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[11]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[1], sd = sqrt(x0[1] * (1-x0[1]) * (t0[11]/(2*N)))))$int)

p01_norm = cbind(dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[2], sd = sqrt(x0[2] * (1-x0[2]) * (t0[11]/(2*N)))))

p02_norm = cbind(dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[3], sd = sqrt(x0[3] * (1-x0[3]) * (t0[11]/(2*N)))))


p03_norm = cbind(dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[4], sd = sqrt(x0[4] * (1-x0[4]) * (t0[11]/(2*N)))))

p04_norm = cbind(dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[5], sd = sqrt(x0[5] * (1-x0[5]) * (t0[11]/(2*N)))))

p05_norm = cbind(dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[6], sd = sqrt(x0[6] * (1-x0[6]) * (t0[11]/(2*N)))))

p06_norm = cbind(dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[7], sd = sqrt(x0[7] * (1-x0[7]) * (t0[11]/(2*N)))))

p07_norm = cbind(dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[8], sd = sqrt(x0[8] * (1-x0[8]) * (t0[11]/(2*N)))))

p08_norm = cbind(dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[9], sd = sqrt(x0[9] * (1-x0[9]) * (t0[11]/(2*N)))))

p09_norm = cbind(dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[1]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[2]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[3]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[4]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[5]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[6]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[7]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[8]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[9]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[10]/(2*N)))),
                 dnorm(x, mean = x0[10], sd = sqrt(x0[10] * (1-x0[10]) * (t0[11]/(2*N)))))

p1_norm = cbind(dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[1]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[1]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[2]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[2]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[3]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[3]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[4]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[4]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[5]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[5]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[6]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[6]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[7]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[7]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[8]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[8]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[9]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[9]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[10]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[10]/(2*N)))))$int,
                dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[11]/(2*N)))) / sintegral(x, dnorm(x, mean = x0[11], sd = sqrt(x0[11] * (1-x0[11]) * (t0[11]/(2*N)))))$int)


hh_nor = matrix(NA, ncol = 11, nrow = 11)

for (i in 1:11) {
  hh_nor[i,1] = hellinger(p0_adaptive[,i],  p0_norm[,i], 0.01, 0.99)
  hh_nor[i,2] = hellinger(p01_adaptive[,i], p01_norm[,i], 0.01, 0.99)
  hh_nor[i,3] = hellinger(p02_adaptive[,i], p02_norm[,i], 0.01, 0.99)
  hh_nor[i,4] = hellinger(p03_adaptive[,i], p03_norm[,i], 0.01, 0.99)
  hh_nor[i,5] = hellinger(p04_adaptive[,i], p04_norm[,i], 0.01, 0.99)
  hh_nor[i,6] = hellinger(p05_adaptive[,i], p05_norm[,i], 0.01, 0.99)
  hh_nor[i,7] = hellinger(p06_adaptive[,i], p06_norm[,i], 0.01, 0.99)
  hh_nor[i,8] = hellinger(p07_adaptive[,i], p07_norm[,i], 0.01, 0.99)
  hh_nor[i,9] = hellinger(p08_adaptive[,i], p08_norm[,i], 0.01, 0.99)
  hh_nor[i,10] = hellinger(p09_adaptive[,i],p09_norm[,i], 0.01, 0.99)
  hh_nor[i,11] = hellinger(p1_adaptive[,i], p1_norm[,i],  0.01, 0.99)
}

colnames(hh_nor) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_nor) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_nor

###normal hellinger 2.1 
hh2_nor = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh2_nor[i,1] = hellinger2(p0_adaptive[,i],  p0_norm[,i], 0.01, 0.99)
  hh2_nor[i,2] = hellinger2(p01_adaptive[,i], p01_norm[,i], 0.01, 0.99)
  hh2_nor[i,3] = hellinger2(p02_adaptive[,i], p02_norm[,i], 0.01, 0.99)
  hh2_nor[i,4] = hellinger2(p03_adaptive[,i], p03_norm[,i], 0.01, 0.99)
  hh2_nor[i,5] = hellinger2(p04_adaptive[,i], p04_norm[,i], 0.01, 0.99)
  hh2_nor[i,6] = hellinger2(p05_adaptive[,i], p05_norm[,i], 0.01, 0.99)
  hh2_nor[i,7] = hellinger2(p06_adaptive[,i], p06_norm[,i], 0.01, 0.99)
  hh2_nor[i,8] = hellinger2(p07_adaptive[,i], p07_norm[,i], 0.01, 0.99)
  hh2_nor[i,9] = hellinger2(p08_adaptive[,i], p08_norm[,i], 0.01, 0.99)
  hh2_nor[i,10] = hellinger2(p09_adaptive[,i],p09_norm[,i], 0.01, 0.99)
  hh2_nor[i,11] = hellinger2(p1_adaptive[,i], p1_norm[,i],  0.01, 0.99)
}

colnames(hh2_nor) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh2_nor) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh2_nor


##4 Gaussian approximation
hh_gauss = matrix(NA, ncol = 11, nrow = 11)

for (i in 1:11) {
  hh_gauss[i,1] = hellinger(p0_adaptive[,i], gauss(x, x0[1], N, t0[i]), 0, 1)
  hh_gauss[i,2] = hellinger(p01_adaptive[,i], gauss(x, x0[2], N, t0[i]), 0, 1)
  hh_gauss[i,3] = hellinger(p02_adaptive[,i], gauss(x, x0[3], N, t0[i]), 0, 1)
  hh_gauss[i,4] = hellinger(p03_adaptive[,i], gauss(x, x0[4], N, t0[i]), 0, 1)
  hh_gauss[i,5] = hellinger(p04_adaptive[,i], gauss(x, x0[5], N, t0[i]), 0, 1)
  hh_gauss[i,6] = hellinger(p05_adaptive[,i], gauss(x, x0[6], N, t0[i]), 0, 1)
  hh_gauss[i,7] = hellinger(p06_adaptive[,i], gauss(x, x0[7], N, t0[i]), 0, 1)
  hh_gauss[i,8] = hellinger(p07_adaptive[,i], gauss(x, x0[8], N, t0[i]), 0, 1)
  hh_gauss[i,9] = hellinger(p08_adaptive[,i], gauss(x, x0[9], N, t0[i]), 0, 1)
  hh_gauss[i,10] = hellinger(p09_adaptive[,i], gauss(x, x0[10], N, t0[i]), 0, 1)
  hh_gauss[i,11] = hellinger(p1_adaptive[,i], gauss(x, x0[11], N, t0[i]), 0, 1)
}
colnames(hh_gauss) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_gauss) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_gauss

##### con la distancia de Hellinger

t_2N = rep(t, 11)
p0_df = c(rep(0,11), rep(0.1, 11), rep(0.2, 11), rep(0.3, 11), rep(0.4, 11), rep(0.5, 11), rep(0.6, 11), rep(0.7, 11), rep(0.8, 11), rep(0.9, 11), rep(1.0, 11))

hell_expan_df = cbind.data.frame(t_2N, p0_df, rep("AE", 121), as.vector(hh_expan))
colnames(hell_expan_df) = c("1","2","3","4")

hell_gauss_df = cbind.data.frame(t_2N, p0_df, rep("GaussA", 121), as.vector(hh_gauss))
colnames(hell_gauss_df) = c("1","2","3","4")

hell_beta_df = cbind.data.frame(t_2N, p0_df, rep("Beta", 121), as.vector(hh_beta))
colnames(hell_beta_df) = c("1","2","3","4")

hell_norm_df = cbind.data.frame(t_2N, p0_df, rep("Gauss", 121), as.vector(hh_nor))
colnames(hell_norm_df) = c("1","2","3","4")

#juntando 
hell_heat_df = rbind(hell_expan_df, hell_norm_df, hell_gauss_df,hell_beta_df)
colnames(hell_heat_df) = c("t", "x0", "dist", "L2")

dim(hell_heat_df)

gg = ggplot(hell_heat_df[1:121,], aes(x=hell_heat_df[1:121,]$x0, y=hell_heat_df[1:121,]$t, fill=hell_heat_df[1:121,]$L2))
gg = gg + geom_tile(color="white", size=0.1)
gg = gg + scale_fill_viridis(name="Hellinger", label=comma)
#gg = gg + coord_equal()
gg = gg + labs(x=NULL, y=NULL, title="Titulo")
gg = gg + theme_tufte(base_family="Arial")
gg = gg + theme(plot.title=element_text(hjust=0))
#gg = gg + theme(axis.ticks = element_blank())
#gg = gg + theme(axis.text = element_text(size=7))

gg = gg + theme(legend.title = element_text(size=8))
gg = gg + theme(legend.text = element_text(size=6))
gg


###juntando 
gg <- ggplot(hell_heat_df, aes(x=x0, y=t, fill=L2))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y="t", x = expression(x[0]))
gg <- gg + scale_fill_viridis(name="Hellinger")
#gg <- gg + coord_equal()
gg <- gg + facet_wrap(~dist, ncol=2)
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
gg



