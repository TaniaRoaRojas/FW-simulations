#direccion de la carpeta
setwd("~/Dropbox/RESEARCH/FISHER WRIGHT DIFFUSION/SIMULACIONES")

###funcion trayectoria

## valores iniciales
#N = 1000 tama??o poblacional
#t = 500 generaciones
#p = 0.5 frecuencia alelica


trayectoria = function(N,p,t){
  frec = rep(0, t)
  frec[1] = p
  
  for (i in 2:t) {
    p = sum( rbinom(n=N, size = 1, prob = p) ) /N
    frec[i] = p
  }
  return(frec)
}

##esta función genera una trayectoria de la variación de la 
##frecuencia alélica. Este proceso se repite 1e^6 veces, 
##se almacenan los resultados en 3 matrices con distintos valores 
##0.1, 0.5 y 0.9 y para diferentes generaciones (11 valores)
##las matrices tienen dimensión 1e^6 x 11.

#p0_inicial01
#p0_inicial05
#p0_inicial09

## Se calcula el asymptotical expansion y gaussian approximation (resultados)
x0 = 0.5
N = 1000
t = 500
x = seq(0,1, length.out = 1000) #para superponer sobre el histograma

## LA formula (asymptotical expansion)
#p0=0.5
expan_1 = (1/sqrt(2*pi*(1/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*1)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_50 = (1/sqrt(2*pi*(50/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*50)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_100 = (1/sqrt(2*pi*(100/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*100)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_150 = (1/sqrt(2*pi*(150/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*150)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_200 = (1/sqrt(2*pi*(200/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*200)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_250 = (1/sqrt(2*pi*(250/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*250)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_300 = (1/sqrt(2*pi*(300/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*300)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_350 = (1/sqrt(2*pi*(350/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*350)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_400 = (1/sqrt(2*pi*(400/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*400)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_450 = (1/sqrt(2*pi*(450/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*450)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_500 = (1/sqrt(2*pi*(500/1000))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*500)) * (asin(1-2*x0) - asin(1-2*x))^(2))

###LA formula (gaussian approximation), p0=0.5
gauss_app_1 = (1/sqrt(2*pi*(1/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*1)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_50 = (1/sqrt(2*pi*(50/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*50)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_100 = (1/sqrt(2*pi*(100/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*100)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_150 = (1/sqrt(2*pi*(150/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*150)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_200 = (1/sqrt(2*pi*(200/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*200)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_250 = (1/sqrt(2*pi*(250/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*250)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_300 = (1/sqrt(2*pi*(300/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*300)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_350 = (1/sqrt(2*pi*(350/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*350)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_400 = (1/sqrt(2*pi*(400/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*400)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_450 = (1/sqrt(2*pi*(450/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*450)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_500 = (1/sqrt(2*pi*(500/1000))) * (1/ sqrt((x0*(1-x0)))) * exp(-(1000/(2*500)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )

#t=1
hist(p0_inicial05[,1], main = "", xlab = "p(t)",breaks = 15, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_1, col = "green", lwd = 2)
lines(x, gauss_app_1, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((1 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=50
hist(p0_inicial05[,2], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F, col = "navy")
lines(x, expan_50, col = "green", lwd = 2)
lines(x, gauss_app_50, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((50 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=100
hist(p0_inicial05[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_100, col = "green", lwd = 2)
lines(x, gauss_app_100, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((100 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#t=150
hist(p0_inicial05[,4], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_150, col = "green", lwd = 2)
lines(x, gauss_app_150, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((100 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=200
hist(p0_inicial05[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_200, col = "green", lwd = 2)
lines(x, gauss_app_200, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((200 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=250
hist(p0_inicial05[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_250, col = "green", lwd = 2)
lines(x, gauss_app_250, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((250 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=300
hist(p0_inicial05[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_300, col = "green", lwd = 2)
lines(x, gauss_app_300, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((300 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=350
hist(p0_inicial05[,8], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_350, col = "green", lwd = 2)
lines(x, gauss_app_350, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((350 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=400
hist(p0_inicial05[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_400, col = "green", lwd = 2)
lines(x, gauss_app_400, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((400 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=450
hist(p0_inicial05[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_450, col = "green", lwd = 2)
lines(x, gauss_app_450, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((450 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#t=500
hist(p0_inicial05[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_500, col = "green", lwd = 2)
lines(x, gauss_app_500, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((500 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

##para estimar, se propone una distribución beta, por lo que se 
##usará el estimador por propuesto por B-K.
##Se crea la funcion beta kernel

beta_kernel = function(x,t,b){
  beta_k = matrix(NA, ncol = length(t), nrow = length(x))
  beta_e = rep(NA, length(t))
  
  for (i in 1:length(x)) {
    for (j in 1:length(t)) {
      
      #kernel
      beta_k[i,j] = (x[i]^(t[j]/b) * (1-x[i])^((1-t[j])/b)) / beta((t[j]/b)+1 , ((1-t[j])/b) +1)
      
    }
  }
  return(beta_k)
}

##y una función para el estimador beta_kernel
beta_est = function(x,b,n){
  y = seq(0,1,length.out = n)
  beta_e = rep(NA, length(y))
  
  for (i in 1:length(y)) {
    beta_e[i] = mean(beta_kernel(x, y[i], b))
  }
  return(cbind(beta_e, y))
}

beta_e_b01_t100 = beta_est(p0_inicial05[,3], 0.1, 100)
beta_e_b03_t100 = beta_est(p0_inicial05[,3], 0.3, 100)
beta_e_b005_t100 = beta_est(p0_inicial05[,3], 0.05, 100)
beta_e_b001_t100 = beta_est(p0_inicial05[,3], 0.01, 100)

beta_e_b01_t450 = beta_est(p0_inicial05[,10], 0.1, 100)
beta_e_b03_t450 = beta_est(p0_inicial05[,10], 0.3, 100)
beta_e_b005_t450 = beta_est(p0_inicial05[,10], 0.05, 100)
beta_e_b001_t450 = beta_est(p0_inicial05[,10], 0.01, 100)


##para la generion t=100 y diferentes valores de b, se obtiene 
par(mfrow=c(1,2))
hist(p0_inicial05[,3], breaks = 35, freq = F, ylim = c(0,4), xlab = expression(p[0]), main = "")
lines(beta_e_b01_t100[,2], beta_e_b01_t100[,1], col = "red", lwd=2)
lines(beta_e_b03_t100[,2], beta_e_b03_t100[,1], col = "green", lwd=2)
lines(beta_e_b005_t100[,2], beta_e_b005_t100[,1], col = "blue", lwd=2)
lines(beta_e_b001_t100[,2], beta_e_b001_t100[,1], col = "orange", lwd=2)
legend("top", legend=c("b=0.1", "b=0.3","b=0.05","b=0.01"),
       col=c("red", "green","blue","orange"), lty=1, cex=0.5, horiz = T)


hist(p0_inicial05[,10], breaks = 35, freq = F, ylim = c(0,4), xlab = expression(p[0]), main = "")
lines(beta_e_b01_t450[,2], beta_e_b01_t450[,1], col = "red", lwd=2)
lines(beta_e_b03_t450[,2], beta_e_b03_t450[,1], col = "green", lwd=2)
lines(beta_e_b005_t450[,2], beta_e_b005_t450[,1], col = "blue", lwd=2)
lines(beta_e_b001_t450[,2], beta_e_b001_t450[,1], col = "orange", lwd=2)
legend("top", legend=c("b=0.1", "b=0.3","b=0.05","b=0.01"),
       col=c("red", "green","blue","orange"), lty=1, cex=0.5, horiz = T)
par(mfrow=c(1,1))
##influye mucho el valor de b, porque lo que a través de un método
##adaptativo, se estimará con el mejor b.
#es necesario definir valores iniciales para la grilla de posibles valores
n = 1e6
varep_1 = 10^(-4) 
gamma_01 = varep_1
gamma_02 = (1/2) * ((1 + varep_1)/(1 - varep_1))

Kn = as.integer(log(n))^2 #169

k_values = seq(0, Kn, 1)

##gamma_k1 (p=1)  
gamma_k1 = c()
for (i in k_values) {
  gamma_k1[i] = gamma_01 + (k_values[i] / Kn) * (2 - gamma_01)
}

## los valores de b tienen que ser menores a 1
gamma_k1 = gamma_k1[gamma_k1 <= 1]
gamma_k1 = gamma_k1[2:85] ##estos son los valores que puede tomar b

#va a 1
plot(n^(-gamma_k1[1:83] / (2*gamma_k1[1:83] + 1)) / n^(-gamma_k1[2:84] / (2*gamma_k1[2:84] + 1)))

##hay 83 valores posibles, para cada uno de esto se calcula 
##el estimador beta_kernel, se almacenan los resultados en 11 (generaciones)
##diferentes matrices de dimensión 68 (elementos de la grilla) y 83 (diferentes
##valor de b).

##INTENTO FALLIDO 1
## se calculó un heatmap donde se ve el MSE entre el estimador beta_kernel
## y gaussian approximation y asymptotical expansion

### comparaciones para t=1 y p0=0.5
beta_est_gammak1_t1 = read.table("est_gamma_k1_t1.txt", sep = " ", header = T)
beta_est_gammak1_t1 = beta_est_gammak1_t1[2:69,2:84]
expan_1 = (1/sqrt(2*pi*(1/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*1)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_1 = expan_1[!is.infinite(expan_1)]
gauss_app_1 = (1/sqrt(2*pi*(1/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*1)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_1 = gauss_app_1[2:69]



### comparaciones para t=50 y p0=0.5
beta_est_gammak1_t50 = read.table("est_gamma_k1_t50.txt", sep = " ", header = T)
beta_est_gammak1_t50 = beta_est_gammak1_t50[2:69,2:84]
expan_50 = (1/sqrt(2*pi*(50/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*50)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_50 = expan_50[!is.infinite(expan_50)]
gauss_app_50 = (1/sqrt(2*pi*(50/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*50)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_50 = gauss_app_50[2:69]

### comparaciones para t=100 y p0=0.5
beta_est_gammak1_t100 = read.table("est_gamma_k1_t100.txt", sep = " ", header = T)
beta_est_gammak1_t100 = beta_est_gammak1_t100[2:69,2:84]
expan_100 = (1/sqrt(2*pi*(100/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*100)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_100 = expan_100[!is.infinite(expan_100)]
gauss_app_100 = (1/sqrt(2*pi*(100/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*100)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_100 = gauss_app_100[2:69]

### comparaciones para t=150 y p0=0.5
beta_est_gammak1_t150 = read.table("est_gamma_k1_t150.txt", sep = " ", header = T)
beta_est_gammak1_t150 = beta_est_gammak1_t150[2:69,2:84]
expan_150 = (1/sqrt(2*pi*(150/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*150)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_150 = expan_150[!is.infinite(expan_150)]
gauss_app_150 = (1/sqrt(2*pi*(150/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*150)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_150 = gauss_app_150[2:69]

### comparaciones para t=200 y p0=0.5
beta_est_gammak1_t200 = read.table("est_gamma_k1_t200.txt", sep = " ", header = T)
beta_est_gammak1_t200 = beta_est_gammak1_t200[2:69,2:84]
expan_200 = (1/sqrt(2*pi*(200/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*200)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_200 = expan_200[!is.infinite(expan_200)]
gauss_app_200 = (1/sqrt(2*pi*(200/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*200)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_200 = gauss_app_200[2:69]

### comparaciones para t=250 y p0=0.5
beta_est_gammak1_t250 = read.table("est_gamma_k1_t250.txt", sep = " ", header = T)
beta_est_gammak1_t250 = beta_est_gammak1_t250[2:69,2:84]
expan_250 = (1/sqrt(2*pi*(250/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*250)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_250 = expan_250[!is.infinite(expan_250)]
gauss_app_250 = (1/sqrt(2*pi*(250/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*250)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_250 = gauss_app_250[2:69]

### comparaciones para t=300 y p0=0.5
beta_est_gammak1_t300 = read.table("est_gamma_k1_t300.txt", sep = " ", header = T)
beta_est_gammak1_t300 = beta_est_gammak1_t300[2:69,2:84]
expan_300 = (1/sqrt(2*pi*(300/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*300)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_300 = expan_300[!is.infinite(expan_300)]
gauss_app_300 = (1/sqrt(2*pi*(300/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*300)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_300 = gauss_app_300[2:69]

### comparaciones para t=350 y p0=0.5
beta_est_gammak1_t350 = read.table("est_gamma_k1_t350.txt", sep = " ", header = T)
beta_est_gammak1_t350 = beta_est_gammak1_t350[2:69,2:84]
expan_350 = (1/sqrt(2*pi*(350/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*350)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_350 = expan_350[!is.infinite(expan_350)]
gauss_app_350 = (1/sqrt(2*pi*(350/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*350)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_350 = gauss_app_350[2:69]

### comparaciones para t=400 y p0=0.5
beta_est_gammak1_t400 = read.table("est_gamma_k1_t400.txt", sep = " ", header = T)
beta_est_gammak1_t400 = beta_est_gammak1_t400[2:69,2:84]
expan_400 = (1/sqrt(2*pi*(400/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*400)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_400 = expan_400[!is.infinite(expan_400)]
gauss_app_400 = (1/sqrt(2*pi*(400/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*400)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_400 = gauss_app_400[2:69]

### comparaciones para t=450 y p0=0.5
beta_est_gammak1_t450 = read.table("est_gamma_k1_t450.txt", sep = " ", header = T)
beta_est_gammak1_t450 = beta_est_gammak1_t450[2:69,2:84]
expan_450 = (1/sqrt(2*pi*(450/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*450)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_450 = expan_450[!is.infinite(expan_450)]
gauss_app_450 = (1/sqrt(2*pi*(450/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*450)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_450 = gauss_app_450[2:69]

### comparaciones para t=500 y p0=0.5
beta_est_gammak1_t500 = read.table("est_gamma_k1_t500.txt", sep = " ", header = T)
beta_est_gammak1_t500 = beta_est_gammak1_t500[2:69,2:84]
expan_500 = (1/sqrt(2*pi*(500/70))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-(70/(2*500)) * (asin(1-2*x0) - asin(1-2*x))^(2))
expan_500 = expan_500[!is.infinite(expan_500)]
gauss_app_500 = (1/sqrt(2*pi*(500/70))) * (1/ sqrt((x0*(1-x0)))) * exp(-(70/(2*500)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
gauss_app_500 = gauss_app_500[2:69]

## matriz con asymptotical expansion
expan = cbind(expan_1, expan_50, expan_100, expan_150, expan_200, expan_250, 
              expan_300, expan_350, expan_400, expan_450, expan_500)

gauss_app = cbind(gauss_app_1, gauss_app_50, gauss_app_100, gauss_app_150,
                  gauss_app_200, gauss_app_250, gauss_app_300, 
                  gauss_app_350, gauss_app_400,
                  gauss_app_450, gauss_app_500)

## hay que formar matrices para los heatmap, cada elemento de la matriz 
## es un MSE, habrán 168 heatmaps (por cada valor gamma_k1 (168))

# matriz de MSE's, para los valores de gamma_k1(fila) y t = 1, 50,.., 500
# comparación con asypmtotical expansion
mse_expan = matrix(0, ncol = 83, nrow = 11)
for (i in 1:11) {
  for (j in 1:83){
    mse_expan[1,j] = mean((beta_est_gammak1_t1[,j] - gauss_app[,1])^2) 
    mse_expan[2,j] = mean((beta_est_gammak1_t50[,j] - gauss_app[,2])^2)  
    mse_expan[3,j] = mean((beta_est_gammak1_t100[,j] - gauss_app[,3])^2) 
    mse_expan[4,j] = mean((beta_est_gammak1_t150[,j] - gauss_app[,4])^2) 
    mse_expan[5,j] = mean((beta_est_gammak1_t200[,j] - gauss_app[,5])^2) 
    mse_expan[6,j] = mean((beta_est_gammak1_t250[,j] - gauss_app[,6])^2)  
    mse_expan[7,j] = mean((beta_est_gammak1_t300[,j] - gauss_app[,7])^2)  
    mse_expan[8,j] = mean((beta_est_gammak1_t350[,j] - gauss_app[,8])^2)  
    mse_expan[9,j] = mean((beta_est_gammak1_t400[,j] - gauss_app[,9])^2) 
    mse_expan[10,j] = mean((beta_est_gammak1_t450[,j] - gauss_app[,10])^2)  
    mse_expan[11,j] = mean((beta_est_gammak1_t500[,j] - gauss_app[,11])^2) 
  }
}
mse_expan


# matriz de MSE's, para los valores de gamma_k1(fila) y t = 1, 50,.., 500
# comparación con gaussian approximation
mse_gauss = matrix(0, ncol = 83, nrow = 11)
for (i in 1:11) {
  for (j in 1:83){
    mse_gauss[1,j] = mean((beta_est_gammak1_t1[,j] - expan[,1])^2)
    mse_gauss[2,j] = mean((beta_est_gammak1_t50[,j] - expan[,2])^2)
    mse_gauss[3,j] = mean((beta_est_gammak1_t100[,j] - expan[,3])^2)
    mse_gauss[4,j] = mean((beta_est_gammak1_t150[,j] - expan[,4])^2) 
    mse_gauss[5,j] = mean((beta_est_gammak1_t200[,j] - expan[,5])^2) 
    mse_gauss[6,j] = mean((beta_est_gammak1_t250[,j] - expan[,6])^2) 
    mse_gauss[7,j] = mean((beta_est_gammak1_t300[,j] - expan[,7])^2) 
    mse_gauss[8,j] = mean((beta_est_gammak1_t350[,j] - expan[,8])^2) 
    mse_gauss[9,j] = mean((beta_est_gammak1_t400[,j] - expan[,9])^2) 
    mse_gauss[10,j] = mean((beta_est_gammak1_t450[,j] - expan[,10])^2) 
    mse_gauss[11,j] = mean((beta_est_gammak1_t500[,j] - expan[,11])^2) 
  }
}
mse_gauss

rownames(mse_expan) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(mse_expan) = paste(expression(gamma[k1]), 2:84, sep="")
mse_expan

rownames(mse_gauss) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(mse_gauss) = paste(expression(gamma[k1]), 2:84, sep="")
mse_gauss

x <- list(
  title = expression(gamma[k1])
)
y <- list(
  title = "t"
)


plot_ly(z = mse_expan, type = "heatmap", 
        x = colnames(mse_expan), y = rownames(mse_expan),
        colors = colorRamp(c("cyan2", "blue4")), 
        colorbar = list(title = "MSE")) %>% layout(xaxis = x, yaxis = y) 

plot_ly(z = mse_gauss, type = "heatmap", 
        x = colnames(mse_gauss), y = rownames(mse_gauss),
        colors = colorRamp(c("cyan2", "blue4"), 1848), 
        colorbar = list(title = "MSE")) %>% layout(xaxis = x, yaxis = y) 

#### encontrar el minimo por fila (por valor de t)

# MSE gaussian approximation
#t=1
gauss_min_b_t1 = sort(mse_gauss[1,])[1]
gauss_min_b_t1 

#t=50
gauss_min_b_t50 = sort(mse_gauss[2,])[1]
gauss_min_b_t50

#t=100
gauss_min_b_t100 = sort(mse_gauss[3,])[1]
gauss_min_b_t100

#t=150
gauss_min_b_t150 = sort(mse_gauss[4,])[1]
gauss_min_b_t150

#t=200
gauss_min_b_t200 = sort(mse_gauss[5,])[1]
gauss_min_b_t200

#t=250
gauss_min_b_t250 = sort(mse_gauss[6,])[1]
gauss_min_b_t250

#t=300
gauss_min_b_t300 = sort(mse_gauss[7,])[1]
gauss_min_b_t300

#t=350
gauss_min_b_t350 = sort(mse_gauss[8,])[1]
gauss_min_b_t350

#t=400
gauss_min_b_t400 = sort(mse_gauss[9,])[1]
gauss_min_b_t400

#t=450
gauss_min_b_t450 = sort(mse_gauss[10,])[1]
gauss_min_b_t450

#t=500
gauss_min_b_t500 = sort(mse_gauss[11,])[1]
gauss_min_b_t500

# MSE asymptotical expansion
#t=1
expan_min_b_t1 = sort(mse_expan[1,])[1]
expan_min_b_t1 


#t=50
expan_min_b_t50 = sort(mse_expan[2,])[1]
expan_min_b_t50

#t=100
expan_min_b_t100 = sort(mse_expan[3,])[1]
expan_min_b_t100

#t=150
expan_min_b_t150 = sort(mse_expan[4,])[1]
expan_min_b_t150

#t=200
expan_min_b_t200 = sort(mse_expan[5,])[1]
expan_min_b_t200

#t=250
expan_min_b_t250 = sort(mse_expan[6,])[1]
expan_min_b_t250

#t=300
expan_min_b_t300 = sort(mse_expan[7,])[1]
expan_min_b_t300

#t=350
expan_min_b_t350 = sort(mse_expan[8,])[1]
expan_min_b_t350

#t=400
expan_min_b_t400 = sort(mse_expan[9,])[1]
expan_min_b_t400

#t=450
expan_min_b_t450 = sort(mse_expan[10,])[1]
expan_min_b_t450

#t=500
expan_min_b_t500 = sort(mse_expan[11,])[1]
expan_min_b_t500

# t = 200
bt200_expan1 = beta_est(p0_inicial05[,5], gamma_k1[84], 50)
bt200_expan2 = beta_est(p0_inicial05[,5], gamma_k1[1], 50)
hist(p0_inicial05[,5], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt200_expan1[,2], bt200_expan1[,1], col = "red", lwd=2)
lines(bt200_expan2[,2], bt200_expan2[,1], col = "blue", lwd=2)

##el valor resultanto ajusta mal

## INTENTO 2. 
## encontrar los índices admisibles de acuerdo al criterio de B-K.

#los valores de gamma_k son
gamma_k1 = c()
for (i in k_values) {
  gamma_k1[i] = gamma_01 + (k_values[i] / Kn) * (2 - gamma_01)
}

## los valores de b tienen que ser menores a 1
gamma_k1 = gamma_k1[gamma_k1 <= 1]
gamma_k1 = gamma_k1[2:85] ##estos son los valores que puede tomar b

#psi(beta) pag3
psi_gamma = n^(-gamma_k1/(2*gamma_k1 + 1))
2 * psi_gamma

#b(beta)
b_beta = n^(-2/(2*gamma_k1 + 1))
b_beta


## ||gamma_m - gamma_l|| < 2*C1*psi(gamma_m), con l < m

#t=50 
#est_gamma_k1_t50 = as.matrix(est_gamma_k1_t50[,2:85])

norm2_t50 = matrix(0, ncol = 84, nrow = 84)

for (l in 1:84) {
  for (m in 1:84) {
    if ((m <= l) & sum(abs(est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])) <= 2*n^(-gamma_k1[m]/(2*gamma_k1[m] + 1))) {
        norm2_t50[i,j] = k
    } else {
        norm2_t50[i,j] = 2
      }
    }
  }


##1 ver que calcule los ECM

#para la desigualdad para cada gamma_k
psi_gamma = n^(-gamma_k1/(2*gamma_k1 + 1))
2 * psi_beta

#t=50
norm2_t50 = matrix(NA, ncol = 84, nrow = 84)
for (m in 1:84) {
  for (l in 1:84) {
    norm2_t50[m,l] = sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]
  }
}
norm2_t50

##para k = 3 fijo
norm2_t50 = matrix(NA, nco = 3, nrow = 3)
for (l in 1:3) {
  for (m in 1:3) {
    if(sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_t50[m,l] = sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2)
    }else{
      norm2_t50[m,l] = NA
    }
  }
}
 
norm2_t50

##para k = 4 fijo
k=4
norm2_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_t50[m,l] = sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2)
    }else{
      norm2_t50[m,l] = NA
    }
    diag(norm2_t50) = NA
  }
}

norm2_t50


##para k = 4 fijo
k=60
norm2_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_t50[m,l] = sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2)
    }else{
      norm2_t50[m,l] = NA
    }
    diag(norm2_t50) = NA
  }
}

norm2_t50

## que me indique la fila y la columna donde se encuentra el máximo?
which(norm2_t50 == max(norm2_t50, na.rm = T), arr.ind=T)

##programa para todos los valores de k
#t=1
sup_k_t1 = rep(NA, 83)
for (k in 2:83) {
  norm2_t1 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t1[,m] - est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
        norm2_t1[m,l] = sum((est_gamma_k1_t1[,m] - est_gamma_k1_t1[,l])^2)
      }else{
        norm2_t1[m,l] = NA
      }
      diag(norm2_t1) = NA
    }
  }
  
  sup_k_t1[k] = which(norm2_t1 == max(norm2_t1, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t1, ylim = c(0,10), type = 'l')

max_k_t1 = max(sup_k_t1, na.rm = T)
max_k_t1

bt1_adaptive = beta_est(p0_inicial05[,1], gamma_k1[max_k_t1], 100)
hist(p0_inicial05[,1], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt1_adaptive[,2], bt1_adaptive[,1], col = "blue", lwd=2)

#t=50
sup_k_t50 = rep(NA, 83)
for (k in 2:83) {
  norm2_t50 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
        norm2_t50[m,l] = sum((est_gamma_k1_t50[,m] - est_gamma_k1_t50[,l])^2)
      }else{
        norm2_t50[m,l] = NA
      }
      diag(norm2_t50) = NA
    }
  }
  
  sup_k_t50[k] = which(norm2_t50 == max(norm2_t50, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t50, ylim = c(0,10), type = 'l')

max_k_t50 = max(sup_k_t50, na.rm = T)
max_k_t50

bt50_adaptive = beta_est(p0_inicial05[,2], gamma_k1[max_k_t50], 100)
hist(p0_inicial05[,3], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt100_adaptive[,2], bt100_adaptive[,1], col = "blue", lwd=2)


#t=100
sup_k_t100 = rep(NA, 83)
for (k in 2:83) {
  norm2_t100 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t100[,m] - est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
        norm2_t100[m,l] = sum((est_gamma_k1_t100[,m] - est_gamma_k1_t100[,l])^2)
      }else{
        norm2_t100[m,l] = NA
      }
      diag(norm2_t100) = NA
    }
  }
  
  sup_k_t100[k] = which(norm2_t100 == max(norm2_t100, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t100, ylim = c(0,10), type = 'l')

max_k_t100 = max(sup_k_t100, na.rm = T)
max_k_t100

bt100_adaptive = beta_est(p0_inicial05[,3], gamma_k1[max_k_t100], 100)
hist(p0_inicial05[,3], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt100_adaptive[,2], bt100_adaptive[,1], col = "blue", lwd=2)

#t=150
sup_k_t150 = rep(NA, 83)
for (k in 2:83) {
  norm2_t150 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t150[,m] - est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
        norm2_t150[m,l] = sum((est_gamma_k1_t150[,m] - est_gamma_k1_t150[,l])^2)
      }else{
        norm2_t150[m,l] = NA
      }
      diag(norm2_t150) = NA
    }
  }
  
  sup_k_t150[k] = which(norm2_t150 == max(norm2_t150, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t150, ylim = c(0,10), type = 'l')

max_k_t150 = max(sup_k_t150, na.rm = T)
max_k_t150

bt150_adaptive = beta_est(p0_inicial05[,4], gamma_k1[max_k_t150], 100)
hist(p0_inicial05[,4], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt150_adaptive[,2], bt150_adaptive[,1], col = "blue", lwd=2)


#t=200
sup_k_t200 = rep(NA, 83)
for (k in 2:83) {
  norm2_t200 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t200[,m] - est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
        norm2_t200[m,l] = sum((est_gamma_k1_t200[,m] - est_gamma_k1_t200[,l])^2)
      }else{
        norm2_t200[m,l] = NA
      }
      diag(norm2_t200) = NA
    }
  }
  
  sup_k_t200[k] = which(norm2_t200 == max(norm2_t200, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t200, ylim = c(0,10), type = 'l')

max_k_t200 = max(sup_k_t200, na.rm = T)
max_k_t200

bt200_adaptive = beta_est(p0_inicial05[,5], gamma_k1[max_k_t200], 100)
hist(p0_inicial05[,5], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt200_adaptive[,2], bt200_adaptive[,1], col = "blue", lwd=2)


#t=250
sup_k_t250 = rep(NA, 83)
for (k in 2:83) {
  norm2_t250 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t250[,m] - est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
        norm2_t250[m,l] = sum((est_gamma_k1_t250[,m] - est_gamma_k1_t250[,l])^2)
      }else{
        norm2_t250[m,l] = NA
      }
      diag(norm2_t250) = NA
    }
  }
  
  sup_k_t250[k] = which(norm2_t250 == max(norm2_t250, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t250, ylim = c(0,10), type = 'l')

max_k_t250 = max(sup_k_t250, na.rm = T)
max_k_t250

bt250_adaptive = beta_est(p0_inicial05[,6], gamma_k1[max_k_t250], 100)
hist(p0_inicial05[,6], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt250_adaptive[,2], bt250_adaptive[,1], col = "blue", lwd=2)


#t=300
sup_k_t300 = rep(NA, 83)
for (k in 2:83) {
  norm2_t300 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t300[,m] - est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
        norm2_t300[m,l] = sum((est_gamma_k1_t300[,m] - est_gamma_k1_t300[,l])^2)
      }else{
        norm2_t300[m,l] = NA
      }
      diag(norm2_t300) = NA
    }
  }
  
  sup_k_t300[k] = which(norm2_t300 == max(norm2_t300, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t300, ylim = c(0,10), type = 'l')

max_k_t300 = max(sup_k_t300, na.rm = T)
max_k_t300

bt300_adaptive = beta_est(p0_inicial05[,7], gamma_k1[max_k_t300], 100)
hist(p0_inicial05[,7], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt300_adaptive[,2], bt300_adaptive[,1], col = "blue", lwd=2)


#t=350
sup_k_t350 = matrix(NA, nrow = 83, ncol = 2)
for (k in 2:83) {
  #k=84
  norm2_t350 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t350[,m] - est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
        norm2_t350[m,l] = sum((est_gamma_k1_t350[,m] - est_gamma_k1_t350[,l])^2)
      }else{
        norm2_t350[m,l] = NA
      }
      diag(norm2_t350) = NA
    }
  }
  
  sup_k_t350[k,] = which(norm2_t350 == max(norm2_t350, na.rm = T), arr.ind=T)
  
}

plot(sup_k_t350, ylim = c(0,10), type = 'l')

max_k_t300 = max(sup_k_t300, na.rm = T)
max_k_t300


k=84
  norm2_t350 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t350[,m] - est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
        norm2_t350[m,l] = sum((est_gamma_k1_t350[,m] - est_gamma_k1_t350[,l])^2)
      }else{
        norm2_t350[m,l] = NA
      }
      diag(norm2_t350) = NA
    }
  }

max_k_t350 = which(norm2_t350 == max(norm2_t350, na.rm = T), arr.ind=T)[1,1]

bt350_adaptive = beta_est(p0_inicial05[,8], gamma_k1[max_k_t350], 100)
hist(p0_inicial05[,8], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt350_adaptive[,2], bt350_adaptive[,1], col = "blue", lwd=2)

#t=400
sup_k_t400 = rep(NA, 83)
for (k in 2:83) {
  norm2_t400 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t400[,m] - est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
        norm2_t400[m,l] = sum((est_gamma_k1_t400[,m] - est_gamma_k1_t400[,l])^2)
      }else{
        norm2_t400[m,l] = NA
      }
      diag(norm2_t400) = NA
    }
  }
  
  sup_k_t400[k] = which(norm2_t400 == max(norm2_t400, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t400, ylim = c(0,10), type = 'l')

max_k_t400 = max(sup_k_t400, na.rm = T)
max_k_t400
k=84
norm2_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((est_gamma_k1_t400[,m] - est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_t400[m,l] = sum((est_gamma_k1_t400[,m] - est_gamma_k1_t400[,l])^2)
    }else{
      norm2_t400[m,l] = NA
    }
    diag(norm2_t400) = NA
  }
}

which(norm2_t400 == max(norm2_t400, na.rm = T), arr.ind=T)
max_k_t400 = which(norm2_t400 == max(norm2_t400, na.rm = T), arr.ind=T)[1,1]

bt400_adaptive = beta_est(p0_inicial05[,9], gamma_k1[max_k_t400], 100)
hist(p0_inicial05[,9], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt400_adaptive[,2], bt400_adaptive[,1], col = "blue", lwd=2)

#t=450
sup_k_t450 = rep(NA, 83)
for (k in 2:83) {
  norm2_t450 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t450[,m] - est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
        norm2_t450[m,l] = sum((est_gamma_k1_t450[,m] - est_gamma_k1_t450[,l])^2)
      }else{
        norm2_t450[m,l] = NA
      }
      diag(norm2_t450) = NA
    }
  }
  
  sup_k_t450[k] = which(norm2_t450 == max(norm2_t450, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t450, ylim = c(0,10), type = 'l')

max_k_t450 = max(sup_k_t450, na.rm = T)
max_k_t450

k=84
norm2_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((est_gamma_k1_t450[,m] - est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_t450[m,l] = sum((est_gamma_k1_t450[,m] - est_gamma_k1_t450[,l])^2)
    }else{
      norm2_t450[m,l] = NA
    }
    diag(norm2_t450) = NA
  }
}

which(norm2_t450 == max(norm2_t450, na.rm = T), arr.ind=T)
max_k_t450 = which(norm2_t450 == max(norm2_t450, na.rm = T), arr.ind=T)[1,1]

bt450_adaptive = beta_est(p0_inicial05[,10], gamma_k1[max_k_t450], 100)
hist(p0_inicial05[,10], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt450_adaptive[,2], bt450_adaptive[,1], col = "blue", lwd=2)


#t=500
sup_k_t500 = rep(NA, 83)
for (k in 2:83) {
  norm2_t500 = matrix(NA, nco = k, nrow = k)
  for (l in 1:k) {
    for (m in 1:k) {
      if((m <= l) & sum((est_gamma_k1_t500[,m] - est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
        norm2_t500[m,l] = sum((est_gamma_k1_t500[,m] - est_gamma_k1_t500[,l])^2)
      }else{
        norm2_t500[m,l] = NA
      }
      diag(norm2_t500) = NA
    }
  }
  
  sup_k_t500[k] = which(norm2_t500 == max(norm2_t500, na.rm = T), arr.ind=T)[1,1]
  
}

plot(sup_k_t500, ylim = c(0,10), type = 'l')

k=84
norm2_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((est_gamma_k1_t500[,m] - est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_t500[m,l] = sum((est_gamma_k1_t500[,m] - est_gamma_k1_t500[,l])^2)
    }else{
      norm2_t500[m,l] = NA
    }
    diag(norm2_t500) = NA
  }
}

which(norm2_t500 == max(norm2_t500, na.rm = T), arr.ind=T)
max_k_t500 = which(norm2_t500 == max(norm2_t500, na.rm = T), arr.ind=T)[1,1]

bt500_adaptive = beta_est(p0_inicial05[,11], gamma_k1[max_k_t500], 100)
hist(p0_inicial05[,10], breaks = 45, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(bt500_adaptive[,2], bt500_adaptive[,1], col = "blue", lwd=2)

