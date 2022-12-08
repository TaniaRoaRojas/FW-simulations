setwd("~/Dropbox/RESEARCH/FISHER WRIGHT DIFFUSION/SIMULACIONES")
#direccion de la carpeta

##valores iniciales de p0
p0_inicial0 = read.table('rep_0.txt', header = T, sep = ' ')
p0_inicial01 = read.table('rep_01.txt', header = T, sep = ' ')
p0_inicial02 = read.table('rep_02.txt', header = T, sep = ' ')
p0_inicial03 = read.table('rep_03.txt', header = T, sep = ' ')
p0_inicial04 = read.table('rep_04.txt', header = T, sep = ' ')
p0_inicial05 = read.table('rep_05.txt', header = T, sep = ' ')
p0_inicial06 = read.table('rep_06.txt', header = T, sep = ' ')
p0_inicial07 = read.table('rep_07.txt', header = T, sep = ' ')
p0_inicial08 = read.table('rep_08.txt', header = T, sep = ' ')
p0_inicial09 = read.table('rep_09.txt', header = T, sep = ' ')
p0_inicial1 = read.table('rep_1.txt', header = T, sep = ' ')

colnames(p0_inicial0) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial01) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial02) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial03) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial04) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial05) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial06) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial07) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial08) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial09) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
colnames(p0_inicial1) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")

## Se calcula el asymptotical expansion y gaussian approximation (resultados)
#x0 = 0.5
N = 70
t = c(1,50,100,150,200,250,300,350,400,450,500)
x0 = seq(0,1, 1/10) #valores iniciales


## LA formula (asymptotical expansion)
expan = function(x0, N, t){
  x = seq(0,1, length.out = N) #para superponer sobre el histograma
  a_expan = (1/sqrt(2*pi*(t/(2*N)))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-((2*N)/(2*t)) * (asin(1-2*x0) - asin(1-2*x))^(2))
  #a_expan = a_expan[!is.infinite(a_expan)] #quita el primero y el último que son Inf
  a_expan[1] = a_expan[N] = 0
  return(a_expan)
}


expan(0.2, 100, 100)

## LA formula (gaussian approximation)
gauss = function(x0, N, t){
  x = seq(0,1, length.out = N) #para superponer sobre el histograma
  gauss_app = (1/sqrt(2*pi*(t/(2*N)))) * (1/ sqrt((x0*(1-x0)))) * exp(-((2*N)/(2*t)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
  return(gauss_app)
}

gauss(0.2, 100, 100)

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


##  encontrar los valores de b óptimos
# hay que encontrar el supremo del conjunto de los valores óptimos

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

#para la desigualdad para cada gamma_k
psi_gamma = n^(-gamma_k1/(2*gamma_k1 + 1))

#b(beta)
b_beta = n^(-2/(2*gamma_k1 + 1))

######################################################################################
######### p0 = 0 ######### 
######### t = 1 ######### 
############################# 
p0_est_gamma_k1_t1 = read.table("p0_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p0_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t1[,m] - p0_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t1[m,l] = sum((p0_est_gamma_k1_t1[,m] - p0_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p0_t1[m,l] = NA
    }
    diag(norm2_p0_t1) = NA
  }
}


max_k_p0_t1 = which(norm2_p0_t1 == max(norm2_p0_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t1

p0_bt1_adaptive = beta_est(p0_inicial0[,1], gamma_k1[max_k_p0_t1], 100)
hist(p0_inicial0[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt1_adaptive[,2], p0_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0 ######### 
######### t = 50 ######### 
############################# 
p0_est_gamma_k1_t50 = read.table("p0_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p0_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t50[,m] - p0_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t50[m,l] = sum((p0_est_gamma_k1_t50[,m] - p0_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p0_t50[m,l] = NA
    }
    diag(norm2_p0_t50) = NA
  }
}


max_k_p0_t50 = which(norm2_p0_t50 == max(norm2_p0_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t50

p0_bt50_adaptive = beta_est(p0_inicial0[,2], gamma_k1[max_k_p0_t50], 100)
hist(p0_inicial0[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt50_adaptive[,2], p0_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0 ######### 
######### t = 100 ######### 
############################# 
p0_est_gamma_k1_t100 = read.table("p0_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p0_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t100[,m] - p0_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t100[m,l] = sum((p0_est_gamma_k1_t100[,m] - p0_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p0_t100[m,l] = NA
    }
    diag(norm2_p0_t100) = NA
  }
}


max_k_p0_t100 = which(norm2_p0_t100 == max(norm2_p0_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t100

p0_bt100_adaptive = beta_est(p0_inicial0[,3], gamma_k1[max_k_p0_t100], 100)
hist(p0_inicial0[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt100_adaptive[,2], p0_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0 ######### 
######### t = 150 ######### 
############################# 
p0_est_gamma_k1_t150 = read.table("p0_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p0_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t150[,m] - p0_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t150[m,l] = sum((p0_est_gamma_k1_t150[,m] - p0_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p0_t150[m,l] = NA
    }
    diag(norm2_p0_t150) = NA
  }
}


max_k_p0_t150 = which(norm2_p0_t150 == max(norm2_p0_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t150

p0_bt150_adaptive = beta_est(p0_inicial0[,4], gamma_k1[max_k_p0_t150], 100)
hist(p0_inicial0[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt150_adaptive[,2], p0_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0 ######### 
######### t = 200 ######### 
############################# 
p0_est_gamma_k1_t200 = read.table("p0_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p0_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t200[,m] - p0_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t200[m,l] = sum((p0_est_gamma_k1_t200[,m] - p0_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p0_t200[m,l] = NA
    }
    diag(norm2_p0_t200) = NA
  }
}


max_k_p0_t200 = which(norm2_p0_t200 == max(norm2_p0_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t200

p0_bt200_adaptive = beta_est(p0_inicial0[,5], gamma_k1[max_k_p0_t200], 100)
hist(p0_inicial0[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt200_adaptive[,2], p0_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0 ######### 
######### t = 250 ######### 
############################# 
p0_est_gamma_k1_t250 = read.table("p0_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p0_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t250[,m] - p0_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t250[m,l] = sum((p0_est_gamma_k1_t250[,m] - p0_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p0_t250[m,l] = NA
    }
    diag(norm2_p0_t250) = NA
  }
}


max_k_p0_t250 = which(norm2_p0_t250 == max(norm2_p0_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t250

p0_bt250_adaptive = beta_est(p0_inicial0[,6], gamma_k1[max_k_p0_t250], 100)
hist(p0_inicial0[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt250_adaptive[,2], p0_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0 ######### 
######### t = 300 ######### 
############################# 
p0_est_gamma_k1_t300 = read.table("p0_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p0_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t300[,m] - p0_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t300[m,l] = sum((p0_est_gamma_k1_t300[,m] - p0_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p0_t300[m,l] = NA
    }
    diag(norm2_p0_t300) = NA
  }
}


max_k_p0_t300 = which(norm2_p0_t300 == max(norm2_p0_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t300

p0_bt300_adaptive = beta_est(p0_inicial0[,7], gamma_k1[max_k_p0_t300], 100)
hist(p0_inicial0[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt300_adaptive[,2], p0_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0 ######### 
######### t = 350 ######### 
############################# 
p0_est_gamma_k1_t350 = read.table("p0_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p0_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t350[,m] - p0_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t350[m,l] = sum((p0_est_gamma_k1_t350[,m] - p0_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p0_t350[m,l] = NA
    }
    diag(norm2_p0_t350) = NA
  }
}


max_k_p0_t350 = which(norm2_p0_t350 == max(norm2_p0_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t350

p0_bt350_adaptive = beta_est(p0_inicial0[,8], gamma_k1[max_k_p0_t350], 100)
hist(p0_inicial0[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt350_adaptive[,2], p0_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0 ######### 
######### t = 400 ######### 
############################# 
p0_est_gamma_k1_t400 = read.table("p0_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p0_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t400[,m] - p0_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t400[m,l] = sum((p0_est_gamma_k1_t400[,m] - p0_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p0_t400[m,l] = NA
    }
    diag(norm2_p0_t400) = NA
  }
}


max_k_p0_t400 = which(norm2_p0_t400 == max(norm2_p0_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t400

p0_bt400_adaptive = beta_est(p0_inicial0[,9], gamma_k1[max_k_p0_t400], 100)
hist(p0_inicial0[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt400_adaptive[,2], p0_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0 ######### 
######### t = 450 ######### 
############################# 
p0_est_gamma_k1_t450 = read.table("p0_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p0_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t450[,m] - p0_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t450[m,l] = sum((p0_est_gamma_k1_t450[,m] - p0_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p0_t450[m,l] = NA
    }
    diag(norm2_p0_t450) = NA
  }
}


max_k_p0_t450 = which(norm2_p0_t450 == max(norm2_p0_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t450

p0_bt450_adaptive = beta_est(p0_inicial0[,10], gamma_k1[max_k_p0_t450], 100)
hist(p0_inicial0[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt450_adaptive[,2], p0_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0 ######### 
######### t = 500 ######### 
############################# 
p0_est_gamma_k1_t500 = read.table("p0_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p0_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p0_est_gamma_k1_t500[,m] - p0_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p0_t500[m,l] = sum((p0_est_gamma_k1_t500[,m] - p0_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p0_t500[m,l] = NA
    }
    diag(norm2_p0_t500) = NA
  }
}


max_k_p0_t500 = which(norm2_p0_t500 == max(norm2_p0_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p0_t500

p0_bt500_adaptive = beta_est(p0_inicial0[,11], gamma_k1[max_k_p0_t500], 100)
hist(p0_inicial0[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p0_bt500_adaptive[,2], p0_bt500_adaptive[,1], col = "blue", lwd=2)


######################################################################################
######### p0 = 0.1 ######### 
######### t = 1 ######### 
############################# 
p01_est_gamma_k1_t1 = read.table("p01_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p01_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t1[,m] - p01_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t1[m,l] = sum((p01_est_gamma_k1_t1[,m] - p01_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p01_t1[m,l] = NA
    }
    diag(norm2_p01_t1) = NA
  }
}


max_k_p01_t1 = which(norm2_p01_t1 == max(norm2_p01_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t1

p01_bt1_adaptive = beta_est(p0_inicial01[,1], gamma_k1[max_k_p01_t1], 100)
hist(p0_inicial01[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt1_adaptive[,2], p01_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.1 ######### 
######### t = 50 ######### 
############################# 
p01_est_gamma_k1_t50 = read.table("p01_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p01_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t50[,m] - p01_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t50[m,l] = sum((p01_est_gamma_k1_t50[,m] - p01_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p01_t50[m,l] = NA
    }
    diag(norm2_p01_t50) = NA
  }
}


max_k_p01_t50 = which(norm2_p01_t50 == max(norm2_p01_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t50

p01_bt50_adaptive = beta_est(p0_inicial01[,2], gamma_k1[max_k_p01_t50], 100)
hist(p0_inicial01[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt50_adaptive[,2], p01_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.1 ######### 
######### t = 100 ######### 
############################# 
p01_est_gamma_k1_t100 = read.table("p01_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p01_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t100[,m] - p01_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t100[m,l] = sum((p01_est_gamma_k1_t100[,m] - p01_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p01_t100[m,l] = NA
    }
    diag(norm2_p01_t100) = NA
  }
}


max_k_p01_t100 = which(norm2_p01_t100 == max(norm2_p01_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t100

p01_bt100_adaptive = beta_est(p0_inicial01[,3], gamma_k1[max_k_p01_t100], 100)
hist(p0_inicial01[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt100_adaptive[,2], p01_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.1 ######### 
######### t = 150 ######### 
############################# 
p01_est_gamma_k1_t150 = read.table("p01_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p01_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t150[,m] - p01_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t150[m,l] = sum((p01_est_gamma_k1_t150[,m] - p01_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p01_t150[m,l] = NA
    }
    diag(norm2_p01_t150) = NA
  }
}


max_k_p01_t150 = which(norm2_p01_t150 == max(norm2_p01_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t150

p01_bt150_adaptive = beta_est(p0_inicial01[,4], gamma_k1[max_k_p01_t150], 100)
hist(p0_inicial01[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt150_adaptive[,2], p01_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.1 ######### 
######### t = 200 ######### 
############################# 
p01_est_gamma_k1_t200 = read.table("p01_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p01_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t200[,m] - p01_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t200[m,l] = sum((p01_est_gamma_k1_t200[,m] - p01_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p01_t200[m,l] = NA
    }
    diag(norm2_p01_t200) = NA
  }
}


max_k_p01_t200 = which(norm2_p01_t200 == max(norm2_p01_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t200

p01_bt200_adaptive = beta_est(p0_inicial01[,5], gamma_k1[max_k_p01_t200], 100)
hist(p0_inicial01[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt200_adaptive[,2], p01_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.1 ######### 
######### t = 250 ######### 
############################# 
p01_est_gamma_k1_t250 = read.table("p01_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p01_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t250[,m] - p01_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t250[m,l] = sum((p01_est_gamma_k1_t250[,m] - p01_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p01_t250[m,l] = NA
    }
    diag(norm2_p01_t250) = NA
  }
}


max_k_p01_t250 = which(norm2_p01_t250 == max(norm2_p01_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t250

p01_bt250_adaptive = beta_est(p0_inicial01[,6], gamma_k1[max_k_p01_t250], 100)
hist(p0_inicial01[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt250_adaptive[,2], p01_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.1 ######### 
######### t = 300 ######### 
############################# 
p01_est_gamma_k1_t300 = read.table("p01_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p01_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t300[,m] - p01_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t300[m,l] = sum((p01_est_gamma_k1_t300[,m] - p01_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p01_t300[m,l] = NA
    }
    diag(norm2_p01_t300) = NA
  }
}


max_k_p01_t300 = which(norm2_p01_t300 == max(norm2_p01_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t300

p01_bt300_adaptive = beta_est(p0_inicial01[,7], gamma_k1[max_k_p01_t300], 100)
hist(p0_inicial01[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt300_adaptive[,2], p01_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.1 ######### 
######### t = 350 ######### 
############################# 
p01_est_gamma_k1_t350 = read.table("p01_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p01_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t350[,m] - p01_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t350[m,l] = sum((p01_est_gamma_k1_t350[,m] - p01_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p01_t350[m,l] = NA
    }
    diag(norm2_p01_t350) = NA
  }
}


max_k_p01_t350 = which(norm2_p01_t350 == max(norm2_p01_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t350

p01_bt350_adaptive = beta_est(p0_inicial01[,8], gamma_k1[max_k_p01_t350], 100)
hist(p0_inicial01[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt350_adaptive[,2], p01_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.1 ######### 
######### t = 400 ######### 
############################# 
p01_est_gamma_k1_t400 = read.table("p01_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p01_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t400[,m] - p01_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t400[m,l] = sum((p01_est_gamma_k1_t400[,m] - p01_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p01_t400[m,l] = NA
    }
    diag(norm2_p01_t400) = NA
  }
}


max_k_p01_t400 = which(norm2_p01_t400 == max(norm2_p01_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t400

p01_bt400_adaptive = beta_est(p0_inicial01[,9], gamma_k1[max_k_p01_t400], 100)
hist(p0_inicial01[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt400_adaptive[,2], p01_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.1 ######### 
######### t = 450 ######### 
############################# 
p01_est_gamma_k1_t450 = read.table("p01_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p01_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t450[,m] - p01_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t450[m,l] = sum((p01_est_gamma_k1_t450[,m] - p01_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p01_t450[m,l] = NA
    }
    diag(norm2_p01_t450) = NA
  }
}


max_k_p01_t450 = which(norm2_p01_t450 == max(norm2_p01_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t450

p01_bt450_adaptive = beta_est(p0_inicial01[,10], gamma_k1[max_k_p01_t450], 100)
hist(p0_inicial01[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt450_adaptive[,2], p01_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.1 ######### 
######### t = 500 ######### 
############################# 
p01_est_gamma_k1_t500 = read.table("p01_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p01_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p01_est_gamma_k1_t500[,m] - p01_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p01_t500[m,l] = sum((p01_est_gamma_k1_t500[,m] - p01_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p01_t500[m,l] = NA
    }
    diag(norm2_p01_t500) = NA
  }
}


max_k_p01_t500 = which(norm2_p01_t500 == max(norm2_p01_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p01_t500

p01_bt500_adaptive = beta_est(p0_inicial01[,11], gamma_k1[max_k_p01_t500], 100)
hist(p0_inicial01[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p01_bt500_adaptive[,2], p01_bt500_adaptive[,1], col = "blue", lwd=2)

###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.2 ######### 
######### t = 1 ######### 
############################# 
p02_est_gamma_k1_t1 = read.table("p02_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p02_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t1[,m] - p02_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t1[m,l] = sum((p02_est_gamma_k1_t1[,m] - p02_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p02_t1[m,l] = NA
    }
    diag(norm2_p02_t1) = NA
  }
}


max_k_p02_t1 = which(norm2_p02_t1 == max(norm2_p02_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t1

p02_bt1_adaptive = beta_est(p0_inicial02[,1], gamma_k1[max_k_p02_t1], 100)
hist(p0_inicial02[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt1_adaptive[,2], p02_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.2 ######### 
######### t = 50 ######### 
############################# 
p02_est_gamma_k1_t50 = read.table("p02_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p02_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t50[,m] - p02_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t50[m,l] = sum((p02_est_gamma_k1_t50[,m] - p02_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p02_t50[m,l] = NA
    }
    diag(norm2_p02_t50) = NA
  }
}


max_k_p02_t50 = which(norm2_p02_t50 == max(norm2_p02_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t50

p02_bt50_adaptive = beta_est(p0_inicial02[,2], gamma_k1[max_k_p02_t50], 100)
hist(p0_inicial02[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt50_adaptive[,2], p02_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.2 ######### 
######### t = 100 ######### 
############################# 
p02_est_gamma_k1_t100 = read.table("p02_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p02_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t100[,m] - p02_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t100[m,l] = sum((p02_est_gamma_k1_t100[,m] - p02_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p02_t100[m,l] = NA
    }
    diag(norm2_p02_t100) = NA
  }
}


max_k_p02_t100 = which(norm2_p02_t100 == max(norm2_p02_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t100

p02_bt100_adaptive = beta_est(p0_inicial02[,3], gamma_k1[max_k_p02_t100], 100)
hist(p0_inicial02[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt100_adaptive[,2], p02_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.2 ######### 
######### t = 150 ######### 
############################# 
p02_est_gamma_k1_t150 = read.table("p02_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p02_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t150[,m] - p02_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t150[m,l] = sum((p02_est_gamma_k1_t150[,m] - p02_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p02_t150[m,l] = NA
    }
    diag(norm2_p02_t150) = NA
  }
}


max_k_p02_t150 = which(norm2_p02_t150 == max(norm2_p02_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t150

p02_bt150_adaptive = beta_est(p0_inicial02[,4], gamma_k1[max_k_p02_t150], 100)
hist(p0_inicial02[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt150_adaptive[,2], p02_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.2 ######### 
######### t = 200 ######### 
############################# 
p02_est_gamma_k1_t200 = read.table("p02_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p02_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t200[,m] - p02_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t200[m,l] = sum((p02_est_gamma_k1_t200[,m] - p02_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p02_t200[m,l] = NA
    }
    diag(norm2_p02_t200) = NA
  }
}


max_k_p02_t200 = which(norm2_p02_t200 == max(norm2_p02_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t200

p02_bt200_adaptive = beta_est(p0_inicial02[,5], gamma_k1[max_k_p02_t200], 100)
hist(p0_inicial02[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt200_adaptive[,2], p02_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.2 ######### 
######### t = 250 ######### 
############################# 
p02_est_gamma_k1_t250 = read.table("p02_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p02_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t250[,m] - p02_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t250[m,l] = sum((p02_est_gamma_k1_t250[,m] - p02_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p02_t250[m,l] = NA
    }
    diag(norm2_p02_t250) = NA
  }
}


max_k_p02_t250 = which(norm2_p02_t250 == max(norm2_p02_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t250

p02_bt250_adaptive = beta_est(p0_inicial02[,6], gamma_k1[max_k_p02_t250], 100)
hist(p0_inicial02[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt250_adaptive[,2], p02_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.2 ######### 
######### t = 300 ######### 
############################# 
p02_est_gamma_k1_t300 = read.table("p02_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p02_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t300[,m] - p02_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t300[m,l] = sum((p02_est_gamma_k1_t300[,m] - p02_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p02_t300[m,l] = NA
    }
    diag(norm2_p02_t300) = NA
  }
}


max_k_p02_t300 = which(norm2_p02_t300 == max(norm2_p02_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t300

p02_bt300_adaptive = beta_est(p0_inicial02[,7], gamma_k1[max_k_p02_t300], 100)
hist(p0_inicial02[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt300_adaptive[,2], p02_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.2 ######### 
######### t = 350 ######### 
############################# 
p02_est_gamma_k1_t350 = read.table("p02_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p02_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t350[,m] - p02_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t350[m,l] = sum((p02_est_gamma_k1_t350[,m] - p02_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p02_t350[m,l] = NA
    }
    diag(norm2_p02_t350) = NA
  }
}


max_k_p02_t350 = which(norm2_p02_t350 == max(norm2_p02_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t350

p02_bt350_adaptive = beta_est(p0_inicial02[,8], gamma_k1[max_k_p02_t350], 100)
hist(p0_inicial02[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt350_adaptive[,2], p02_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.2 ######### 
######### t = 400 ######### 
############################# 
p02_est_gamma_k1_t400 = read.table("p02_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p02_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t400[,m] - p02_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t400[m,l] = sum((p02_est_gamma_k1_t400[,m] - p02_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p02_t400[m,l] = NA
    }
    diag(norm2_p02_t400) = NA
  }
}


max_k_p02_t400 = which(norm2_p02_t400 == max(norm2_p02_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t400

p02_bt400_adaptive = beta_est(p0_inicial02[,9], gamma_k1[max_k_p02_t400], 100)
hist(p0_inicial02[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt400_adaptive[,2], p02_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.2 ######### 
######### t = 450 ######### 
############################# 
p02_est_gamma_k1_t450 = read.table("p02_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p02_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t450[,m] - p02_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t450[m,l] = sum((p02_est_gamma_k1_t450[,m] - p02_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p02_t450[m,l] = NA
    }
    diag(norm2_p02_t450) = NA
  }
}


max_k_p02_t450 = which(norm2_p02_t450 == max(norm2_p02_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t450

p02_bt450_adaptive = beta_est(p0_inicial02[,10], gamma_k1[max_k_p02_t450], 100)
hist(p0_inicial02[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt450_adaptive[,2], p02_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.2 ######### 
######### t = 500 ######### 
############################# 
p02_est_gamma_k1_t500 = read.table("p02_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p02_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p02_est_gamma_k1_t500[,m] - p02_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p02_t500[m,l] = sum((p02_est_gamma_k1_t500[,m] - p02_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p02_t500[m,l] = NA
    }
    diag(norm2_p02_t500) = NA
  }
}


max_k_p02_t500 = which(norm2_p02_t500 == max(norm2_p02_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p02_t500

p02_bt500_adaptive = beta_est(p0_inicial02[,11], gamma_k1[max_k_p02_t500], 100)
hist(p0_inicial02[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p02_bt500_adaptive[,2], p02_bt500_adaptive[,1], col = "blue", lwd=2)


###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.3 ######### 
######### t = 1 ######### 
############################# 
p03_est_gamma_k1_t1 = read.table("p03_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p03_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t1[,m] - p03_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t1[m,l] = sum((p03_est_gamma_k1_t1[,m] - p03_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p03_t1[m,l] = NA
    }
    diag(norm2_p03_t1) = NA
  }
}


max_k_p03_t1 = which(norm2_p03_t1 == max(norm2_p03_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t1

p03_bt1_adaptive = beta_est(p0_inicial03[,1], gamma_k1[max_k_p03_t1], 100)
hist(p0_inicial03[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt1_adaptive[,2], p03_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.3 ######### 
######### t = 50 ######### 
############################# 
p03_est_gamma_k1_t50 = read.table("p03_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p03_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t50[,m] - p03_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t50[m,l] = sum((p03_est_gamma_k1_t50[,m] - p03_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p03_t50[m,l] = NA
    }
    diag(norm2_p03_t50) = NA
  }
}


max_k_p03_t50 = which(norm2_p03_t50 == max(norm2_p03_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t50

p03_bt50_adaptive = beta_est(p0_inicial03[,2], gamma_k1[max_k_p03_t50], 100)
hist(p0_inicial03[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt50_adaptive[,2], p03_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.3 ######### 
######### t = 100 ######### 
############################# 
p03_est_gamma_k1_t100 = read.table("p03_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p03_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t100[,m] - p03_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t100[m,l] = sum((p03_est_gamma_k1_t100[,m] - p03_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p03_t100[m,l] = NA
    }
    diag(norm2_p03_t100) = NA
  }
}


max_k_p03_t100 = which(norm2_p03_t100 == max(norm2_p03_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t100

p03_bt100_adaptive = beta_est(p0_inicial03[,3], gamma_k1[max_k_p03_t100], 100)
hist(p0_inicial03[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt100_adaptive[,2], p03_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.3 ######### 
######### t = 150 ######### 
############################# 
p03_est_gamma_k1_t150 = read.table("p03_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p03_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t150[,m] - p03_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t150[m,l] = sum((p03_est_gamma_k1_t150[,m] - p03_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p03_t150[m,l] = NA
    }
    diag(norm2_p03_t150) = NA
  }
}


max_k_p03_t150 = which(norm2_p03_t150 == max(norm2_p03_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t150

p03_bt150_adaptive = beta_est(p0_inicial03[,4], gamma_k1[max_k_p03_t150], 100)
hist(p0_inicial03[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt150_adaptive[,2], p03_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.3 ######### 
######### t = 200 ######### 
############################# 
p03_est_gamma_k1_t200 = read.table("p03_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p03_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t200[,m] - p03_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t200[m,l] = sum((p03_est_gamma_k1_t200[,m] - p03_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p03_t200[m,l] = NA
    }
    diag(norm2_p03_t200) = NA
  }
}


max_k_p03_t200 = which(norm2_p03_t200 == max(norm2_p03_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t200

p03_bt200_adaptive = beta_est(p0_inicial03[,5], gamma_k1[max_k_p03_t200], 100)
hist(p0_inicial03[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt200_adaptive[,2], p03_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.3 ######### 
######### t = 250 ######### 
############################# 
p03_est_gamma_k1_t250 = read.table("p03_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p03_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t250[,m] - p03_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t250[m,l] = sum((p03_est_gamma_k1_t250[,m] - p03_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p03_t250[m,l] = NA
    }
    diag(norm2_p03_t250) = NA
  }
}


max_k_p03_t250 = which(norm2_p03_t250 == max(norm2_p03_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t250

p03_bt250_adaptive = beta_est(p0_inicial03[,6], gamma_k1[max_k_p03_t250], 100)
hist(p0_inicial03[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt250_adaptive[,2], p03_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.3 ######### 
######### t = 300 ######### 
############################# 
p03_est_gamma_k1_t300 = read.table("p03_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p03_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t300[,m] - p03_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t300[m,l] = sum((p03_est_gamma_k1_t300[,m] - p03_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p03_t300[m,l] = NA
    }
    diag(norm2_p03_t300) = NA
  }
}


max_k_p03_t300 = which(norm2_p03_t300 == max(norm2_p03_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t300

p03_bt300_adaptive = beta_est(p0_inicial03[,7], gamma_k1[max_k_p03_t300], 100)
hist(p0_inicial03[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt300_adaptive[,2], p03_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.3 ######### 
######### t = 350 ######### 
############################# 
p03_est_gamma_k1_t350 = read.table("p03_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p03_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t350[,m] - p03_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t350[m,l] = sum((p03_est_gamma_k1_t350[,m] - p03_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p03_t350[m,l] = NA
    }
    diag(norm2_p03_t350) = NA
  }
}


max_k_p03_t350 = which(norm2_p03_t350 == max(norm2_p03_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t350

p03_bt350_adaptive = beta_est(p0_inicial03[,8], gamma_k1[max_k_p03_t350], 100)
hist(p0_inicial03[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt350_adaptive[,2], p03_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.3 ######### 
######### t = 400 ######### 
############################# 
p03_est_gamma_k1_t400 = read.table("p03_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p03_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t400[,m] - p03_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t400[m,l] = sum((p03_est_gamma_k1_t400[,m] - p03_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p03_t400[m,l] = NA
    }
    diag(norm2_p03_t400) = NA
  }
}


max_k_p03_t400 = which(norm2_p03_t400 == max(norm2_p03_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t400

p03_bt400_adaptive = beta_est(p0_inicial03[,9], gamma_k1[max_k_p03_t400], 100)
hist(p0_inicial03[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt400_adaptive[,2], p03_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.3 ######### 
######### t = 450 ######### 
############################# 
p03_est_gamma_k1_t450 = read.table("p03_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p03_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t450[,m] - p03_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t450[m,l] = sum((p03_est_gamma_k1_t450[,m] - p03_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p03_t450[m,l] = NA
    }
    diag(norm2_p03_t450) = NA
  }
}


max_k_p03_t450 = which(norm2_p03_t450 == max(norm2_p03_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t450

p03_bt450_adaptive = beta_est(p0_inicial03[,10], gamma_k1[max_k_p03_t450], 100)
hist(p0_inicial03[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt450_adaptive[,2], p03_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.3 ######### 
######### t = 500 ######### 
############################# 
p03_est_gamma_k1_t500 = read.table("p03_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p03_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p03_est_gamma_k1_t500[,m] - p03_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p03_t500[m,l] = sum((p03_est_gamma_k1_t500[,m] - p03_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p03_t500[m,l] = NA
    }
    diag(norm2_p03_t500) = NA
  }
}


max_k_p03_t500 = which(norm2_p03_t500 == max(norm2_p03_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p03_t500

p03_bt500_adaptive = beta_est(p0_inicial03[,11], gamma_k1[max_k_p03_t500], 100)
hist(p0_inicial03[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p03_bt500_adaptive[,2], p03_bt500_adaptive[,1], col = "blue", lwd=2)


###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.4 ######### 
######### t = 1 ######### 
############################# 
p04_est_gamma_k1_t1 = read.table("p04_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p04_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t1[,m] - p04_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t1[m,l] = sum((p04_est_gamma_k1_t1[,m] - p04_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p04_t1[m,l] = NA
    }
    diag(norm2_p04_t1) = NA
  }
}


max_k_p04_t1 = which(norm2_p04_t1 == max(norm2_p04_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t1

p04_bt1_adaptive = beta_est(p0_inicial04[,1], gamma_k1[max_k_p04_t1], 100)
hist(p0_inicial04[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt1_adaptive[,2], p04_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.4 ######### 
######### t = 50 ######### 
############################# 
p04_est_gamma_k1_t50 = read.table("p04_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p04_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t50[,m] - p04_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t50[m,l] = sum((p04_est_gamma_k1_t50[,m] - p04_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p04_t50[m,l] = NA
    }
    diag(norm2_p04_t50) = NA
  }
}


max_k_p04_t50 = which(norm2_p04_t50 == max(norm2_p04_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t50

p04_bt50_adaptive = beta_est(p0_inicial04[,2], gamma_k1[max_k_p04_t50], 100)
hist(p0_inicial04[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt50_adaptive[,2], p04_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.4 ######### 
######### t = 100 ######### 
############################# 
p04_est_gamma_k1_t100 = read.table("p04_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p04_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t100[,m] - p04_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t100[m,l] = sum((p04_est_gamma_k1_t100[,m] - p04_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p04_t100[m,l] = NA
    }
    diag(norm2_p04_t100) = NA
  }
}


max_k_p04_t100 = which(norm2_p04_t100 == max(norm2_p04_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t100

p04_bt100_adaptive = beta_est(p0_inicial04[,3], gamma_k1[max_k_p04_t100], 100)
hist(p0_inicial04[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt100_adaptive[,2], p04_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.4 ######### 
######### t = 150 ######### 
############################# 
p04_est_gamma_k1_t150 = read.table("p04_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p04_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t150[,m] - p04_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t150[m,l] = sum((p04_est_gamma_k1_t150[,m] - p04_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p04_t150[m,l] = NA
    }
    diag(norm2_p04_t150) = NA
  }
}


max_k_p04_t150 = which(norm2_p04_t150 == max(norm2_p04_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t150

p04_bt150_adaptive = beta_est(p0_inicial04[,4], gamma_k1[max_k_p04_t150], 100)
hist(p0_inicial04[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt150_adaptive[,2], p04_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.4 ######### 
######### t = 200 ######### 
############################# 
p04_est_gamma_k1_t200 = read.table("p04_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p04_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t200[,m] - p04_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t200[m,l] = sum((p04_est_gamma_k1_t200[,m] - p04_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p04_t200[m,l] = NA
    }
    diag(norm2_p04_t200) = NA
  }
}


max_k_p04_t200 = which(norm2_p04_t200 == max(norm2_p04_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t200

p04_bt200_adaptive = beta_est(p0_inicial04[,5], gamma_k1[max_k_p04_t200], 100)
hist(p0_inicial04[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt200_adaptive[,2], p04_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.4 ######### 
######### t = 250 ######### 
############################# 
p04_est_gamma_k1_t250 = read.table("p04_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p04_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t250[,m] - p04_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t250[m,l] = sum((p04_est_gamma_k1_t250[,m] - p04_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p04_t250[m,l] = NA
    }
    diag(norm2_p04_t250) = NA
  }
}


max_k_p04_t250 = which(norm2_p04_t250 == max(norm2_p04_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t250

p04_bt250_adaptive = beta_est(p0_inicial04[,6], gamma_k1[max_k_p04_t250], 100)
hist(p0_inicial04[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt250_adaptive[,2], p04_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.4 ######### 
######### t = 300 ######### 
############################# 
p04_est_gamma_k1_t300 = read.table("p04_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p04_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t300[,m] - p04_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t300[m,l] = sum((p04_est_gamma_k1_t300[,m] - p04_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p04_t300[m,l] = NA
    }
    diag(norm2_p04_t300) = NA
  }
}


max_k_p04_t300 = which(norm2_p04_t300 == max(norm2_p04_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t300

p04_bt300_adaptive = beta_est(p0_inicial04[,7], gamma_k1[max_k_p04_t300], 100)
hist(p0_inicial04[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt300_adaptive[,2], p04_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.4 ######### 
######### t = 350 ######### 
############################# 
p04_est_gamma_k1_t350 = read.table("p04_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p04_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t350[,m] - p04_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t350[m,l] = sum((p04_est_gamma_k1_t350[,m] - p04_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p04_t350[m,l] = NA
    }
    diag(norm2_p04_t350) = NA
  }
}


max_k_p04_t350 = which(norm2_p04_t350 == max(norm2_p04_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t350

p04_bt350_adaptive = beta_est(p0_inicial04[,8], gamma_k1[max_k_p04_t350], 100)
hist(p0_inicial04[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt350_adaptive[,2], p04_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.4 ######### 
######### t = 400 ######### 
############################# 
p04_est_gamma_k1_t400 = read.table("p04_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p04_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t400[,m] - p04_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t400[m,l] = sum((p04_est_gamma_k1_t400[,m] - p04_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p04_t400[m,l] = NA
    }
    diag(norm2_p04_t400) = NA
  }
}


max_k_p04_t400 = which(norm2_p04_t400 == max(norm2_p04_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t400

p04_bt400_adaptive = beta_est(p0_inicial04[,9], gamma_k1[max_k_p04_t400], 100)
hist(p0_inicial04[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt400_adaptive[,2], p04_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.4 ######### 
######### t = 450 ######### 
############################# 
p04_est_gamma_k1_t450 = read.table("p04_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p04_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t450[,m] - p04_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t450[m,l] = sum((p04_est_gamma_k1_t450[,m] - p04_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p04_t450[m,l] = NA
    }
    diag(norm2_p04_t450) = NA
  }
}


max_k_p04_t450 = which(norm2_p04_t450 == max(norm2_p04_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t450

p04_bt450_adaptive = beta_est(p0_inicial04[,10], gamma_k1[max_k_p04_t450], 100)
hist(p0_inicial04[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt450_adaptive[,2], p04_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.4 ######### 
######### t = 500 ######### 
############################# 
p04_est_gamma_k1_t500 = read.table("p04_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p04_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p04_est_gamma_k1_t500[,m] - p04_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p04_t500[m,l] = sum((p04_est_gamma_k1_t500[,m] - p04_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p04_t500[m,l] = NA
    }
    diag(norm2_p04_t500) = NA
  }
}


max_k_p04_t500 = which(norm2_p04_t500 == max(norm2_p04_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p04_t500

p04_bt500_adaptive = beta_est(p0_inicial04[,11], gamma_k1[max_k_p04_t500], 100)
hist(p0_inicial04[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p04_bt500_adaptive[,2], p04_bt500_adaptive[,1], col = "blue", lwd=2)

###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.5 ######### 
######### t = 1 ######### 
############################# 
p05_est_gamma_k1_t1 = read.table("p05_est_gamma_k1_t1.txt", header = T, sep = ' ')
p05_est_gamma_k1_t1 = p05_est_gamma_k1_t1[,2:85]

k=84
norm2_p05_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t1[,m] - p05_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t1[m,l] = sum((p05_est_gamma_k1_t1[,m] - p05_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p05_t1[m,l] = NA
    }
    diag(norm2_p05_t1) = NA
  }
}


max_k_p05_t1 = which(norm2_p05_t1 == max(norm2_p05_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t1

p05_bt1_adaptive = beta_est(p0_inicial05[,1], gamma_k1[max_k_p05_t1], 100)
hist(p0_inicial05[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt1_adaptive[,2], p05_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.5 ######### 
######### t = 50 ######### 
############################# 
p05_est_gamma_k1_t50 = read.table("p05_est_gamma_k1_t50.txt", header = T, sep = ' ')
p05_est_gamma_k1_t50 = p05_est_gamma_k1_t50[,2:85]

k=84
norm2_p05_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t50[,m] - p05_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t50[m,l] = sum((p05_est_gamma_k1_t50[,m] - p05_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p05_t50[m,l] = NA
    }
    diag(norm2_p05_t50) = NA
  }
}


max_k_p05_t50 = which(norm2_p05_t50 == max(norm2_p05_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t50

p05_bt50_adaptive = beta_est(p0_inicial05[,2], gamma_k1[max_k_p05_t50], 100)
hist(p0_inicial05[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt50_adaptive[,2], p05_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.5 ######### 
######### t = 100 ######### 
############################# 
p05_est_gamma_k1_t100 = read.table("p05_est_gamma_k1_t100.txt", header = T, sep = ' ')
p05_est_gamma_k1_t100 = p05_est_gamma_k1_t100[,2:85]

k=84
norm2_p05_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t100[,m] - p05_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t100[m,l] = sum((p05_est_gamma_k1_t100[,m] - p05_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p05_t100[m,l] = NA
    }
    diag(norm2_p05_t100) = NA
  }
}


max_k_p05_t100 = which(norm2_p05_t100 == max(norm2_p05_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t100

p05_bt100_adaptive = beta_est(p0_inicial05[,3], gamma_k1[max_k_p05_t100], 100)
hist(p0_inicial05[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt100_adaptive[,2], p05_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.5 ######### 
######### t = 150 ######### 
############################# 
p05_est_gamma_k1_t150 = read.table("p05_est_gamma_k1_t150.txt", header = T, sep = ' ')
p05_est_gamma_k1_t150 = p05_est_gamma_k1_t150[,2:85]

k=84
norm2_p05_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t150[,m] - p05_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t150[m,l] = sum((p05_est_gamma_k1_t150[,m] - p05_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p05_t150[m,l] = NA
    }
    diag(norm2_p05_t150) = NA
  }
}


max_k_p05_t150 = which(norm2_p05_t150 == max(norm2_p05_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t150

p05_bt150_adaptive = beta_est(p0_inicial05[,4], gamma_k1[max_k_p05_t150], 100)
hist(p0_inicial05[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt150_adaptive[,2], p05_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.5 ######### 
######### t = 200 ######### 
############################# 
p05_est_gamma_k1_t200 = read.table("p05_est_gamma_k1_t200.txt", header = T, sep = ' ')
p05_est_gamma_k1_t200 = p05_est_gamma_k1_t200[,2:85]

k=84
norm2_p05_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t200[,m] - p05_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t200[m,l] = sum((p05_est_gamma_k1_t200[,m] - p05_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p05_t200[m,l] = NA
    }
    diag(norm2_p05_t200) = NA
  }
}


max_k_p05_t200 = which(norm2_p05_t200 == max(norm2_p05_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t200

p05_bt200_adaptive = beta_est(p0_inicial05[,5], gamma_k1[max_k_p05_t200], 100)
hist(p0_inicial05[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt200_adaptive[,2], p05_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.5 ######### 
######### t = 250 ######### 
############################# 
p05_est_gamma_k1_t250 = read.table("p05_est_gamma_k1_t250.txt", header = T, sep = ' ')
p05_est_gamma_k1_t250 = p05_est_gamma_k1_t250[,2:85]

k=84
norm2_p05_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t250[,m] - p05_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t250[m,l] = sum((p05_est_gamma_k1_t250[,m] - p05_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p05_t250[m,l] = NA
    }
    diag(norm2_p05_t250) = NA
  }
}


max_k_p05_t250 = which(norm2_p05_t250 == max(norm2_p05_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t250

p05_bt250_adaptive = beta_est(p0_inicial05[,6], gamma_k1[max_k_p05_t250], 100)
hist(p0_inicial05[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt250_adaptive[,2], p05_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.5 ######### 
######### t = 300 ######### 
############################# 
p05_est_gamma_k1_t300 = read.table("p05_est_gamma_k1_t300.txt", header = T, sep = ' ')
p05_est_gamma_k1_t300 = p05_est_gamma_k1_t300[,2:85]

k=84
norm2_p05_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t300[,m] - p05_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t300[m,l] = sum((p05_est_gamma_k1_t300[,m] - p05_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p05_t300[m,l] = NA
    }
    diag(norm2_p05_t300) = NA
  }
}


max_k_p05_t300 = which(norm2_p05_t300 == max(norm2_p05_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t300

p05_bt300_adaptive = beta_est(p0_inicial05[,7], gamma_k1[max_k_p05_t300], 100)
hist(p0_inicial05[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt300_adaptive[,2], p05_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.5 ######### 
######### t = 350 ######### 
############################# 
p05_est_gamma_k1_t350 = read.table("p05_est_gamma_k1_t350.txt", header = T, sep = ' ')
p05_est_gamma_k1_t350 = p05_est_gamma_k1_t350[,2:85]

k=84
norm2_p05_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t350[,m] - p05_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t350[m,l] = sum((p05_est_gamma_k1_t350[,m] - p05_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p05_t350[m,l] = NA
    }
    diag(norm2_p05_t350) = NA
  }
}


max_k_p05_t350 = which(norm2_p05_t350 == max(norm2_p05_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t350

p05_bt350_adaptive = beta_est(p0_inicial05[,8], gamma_k1[max_k_p05_t350], 100)
hist(p0_inicial05[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt350_adaptive[,2], p05_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.5 ######### 
######### t = 400 ######### 
############################# 
p05_est_gamma_k1_t400 = read.table("p05_est_gamma_k1_t400.txt", header = T, sep = ' ')
p05_est_gamma_k1_t400 = p05_est_gamma_k1_t400[,2:85]

k=84
norm2_p05_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t400[,m] - p05_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t400[m,l] = sum((p05_est_gamma_k1_t400[,m] - p05_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p05_t400[m,l] = NA
    }
    diag(norm2_p05_t400) = NA
  }
}


max_k_p05_t400 = which(norm2_p05_t400 == max(norm2_p05_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t400

p05_bt400_adaptive = beta_est(p0_inicial05[,9], gamma_k1[max_k_p05_t400], 100)
hist(p0_inicial05[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt400_adaptive[,2], p05_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.5 ######### 
######### t = 450 ######### 
############################# 
p05_est_gamma_k1_t450 = read.table("p05_est_gamma_k1_t450.txt", header = T, sep = ' ')
p05_est_gamma_k1_t450 = p05_est_gamma_k1_t450[,2:85]

k=84
norm2_p05_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t450[,m] - p05_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t450[m,l] = sum((p05_est_gamma_k1_t450[,m] - p05_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p05_t450[m,l] = NA
    }
    diag(norm2_p05_t450) = NA
  }
}


max_k_p05_t450 = which(norm2_p05_t450 == max(norm2_p05_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t450

p05_bt450_adaptive = beta_est(p0_inicial05[,10], gamma_k1[max_k_p05_t450], 100)
hist(p0_inicial05[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt450_adaptive[,2], p05_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.5 ######### 
######### t = 500 ######### 
############################# 
p05_est_gamma_k1_t500 = read.table("p05_est_gamma_k1_t500.txt", header = T, sep = ' ')
p05_est_gamma_k1_t500 = p05_est_gamma_k1_t500[,2:85]

k=84
norm2_p05_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p05_est_gamma_k1_t500[,m] - p05_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p05_t500[m,l] = sum((p05_est_gamma_k1_t500[,m] - p05_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p05_t500[m,l] = NA
    }
    diag(norm2_p05_t500) = NA
  }
}


max_k_p05_t500 = which(norm2_p05_t500 == max(norm2_p05_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p05_t500

p05_bt500_adaptive = beta_est(p0_inicial05[,11], gamma_k1[max_k_p05_t500], 100)
hist(p0_inicial05[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p05_bt500_adaptive[,2], p05_bt500_adaptive[,1], col = "blue", lwd=2)

###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.6 ######### 
######### t = 1 ######### 
############################# 
p06_est_gamma_k1_t1 = read.table("p06_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p06_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t1[,m] - p06_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t1[m,l] = sum((p06_est_gamma_k1_t1[,m] - p06_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p06_t1[m,l] = NA
    }
    diag(norm2_p06_t1) = NA
  }
}


max_k_p06_t1 = which(norm2_p06_t1 == max(norm2_p06_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t1

p06_bt1_adaptive = beta_est(p0_inicial06[,1], gamma_k1[max_k_p06_t1], 100)
hist(p0_inicial06[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt1_adaptive[,2], p06_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.6 ######### 
######### t = 50 ######### 
############################# 
p06_est_gamma_k1_t50 = read.table("p06_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p06_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t50[,m] - p06_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t50[m,l] = sum((p06_est_gamma_k1_t50[,m] - p06_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p06_t50[m,l] = NA
    }
    diag(norm2_p06_t50) = NA
  }
}


max_k_p06_t50 = which(norm2_p06_t50 == max(norm2_p06_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t50

p06_bt50_adaptive = beta_est(p0_inicial06[,2], gamma_k1[max_k_p06_t50], 100)
hist(p0_inicial06[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt50_adaptive[,2], p06_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.6 ######### 
######### t = 100 ######### 
############################# 
p06_est_gamma_k1_t100 = read.table("p06_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p06_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t100[,m] - p06_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t100[m,l] = sum((p06_est_gamma_k1_t100[,m] - p06_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p06_t100[m,l] = NA
    }
    diag(norm2_p06_t100) = NA
  }
}


max_k_p06_t100 = which(norm2_p06_t100 == max(norm2_p06_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t100

p06_bt100_adaptive = beta_est(p0_inicial06[,3], gamma_k1[max_k_p06_t100], 100)
hist(p0_inicial06[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt100_adaptive[,2], p06_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.6 ######### 
######### t = 150 ######### 
############################# 
p06_est_gamma_k1_t150 = read.table("p06_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p06_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t150[,m] - p06_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t150[m,l] = sum((p06_est_gamma_k1_t150[,m] - p06_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p06_t150[m,l] = NA
    }
    diag(norm2_p06_t150) = NA
  }
}


max_k_p06_t150 = which(norm2_p06_t150 == max(norm2_p06_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t150

p06_bt150_adaptive = beta_est(p0_inicial06[,4], gamma_k1[max_k_p06_t150], 100)
hist(p0_inicial06[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt150_adaptive[,2], p06_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.6 ######### 
######### t = 200 ######### 
############################# 
p06_est_gamma_k1_t200 = read.table("p06_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p06_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t200[,m] - p06_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t200[m,l] = sum((p06_est_gamma_k1_t200[,m] - p06_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p06_t200[m,l] = NA
    }
    diag(norm2_p06_t200) = NA
  }
}


max_k_p06_t200 = which(norm2_p06_t200 == max(norm2_p06_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t200

p06_bt200_adaptive = beta_est(p0_inicial06[,5], gamma_k1[max_k_p06_t200], 100)
hist(p0_inicial06[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt200_adaptive[,2], p06_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.6 ######### 
######### t = 250 ######### 
############################# 
p06_est_gamma_k1_t250 = read.table("p06_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p06_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t250[,m] - p06_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t250[m,l] = sum((p06_est_gamma_k1_t250[,m] - p06_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p06_t250[m,l] = NA
    }
    diag(norm2_p06_t250) = NA
  }
}


max_k_p06_t250 = which(norm2_p06_t250 == max(norm2_p06_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t250

p06_bt250_adaptive = beta_est(p0_inicial06[,6], gamma_k1[max_k_p06_t250], 100)
hist(p0_inicial06[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt250_adaptive[,2], p06_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.6 ######### 
######### t = 300 ######### 
############################# 
p06_est_gamma_k1_t300 = read.table("p06_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p06_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t300[,m] - p06_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t300[m,l] = sum((p06_est_gamma_k1_t300[,m] - p06_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p06_t300[m,l] = NA
    }
    diag(norm2_p06_t300) = NA
  }
}


max_k_p06_t300 = which(norm2_p06_t300 == max(norm2_p06_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t300

p06_bt300_adaptive = beta_est(p0_inicial06[,7], gamma_k1[max_k_p06_t300], 100)
hist(p0_inicial06[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt300_adaptive[,2], p06_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.6 ######### 
######### t = 350 ######### 
############################# 
p06_est_gamma_k1_t350 = read.table("p06_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p06_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t350[,m] - p06_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t350[m,l] = sum((p06_est_gamma_k1_t350[,m] - p06_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p06_t350[m,l] = NA
    }
    diag(norm2_p06_t350) = NA
  }
}


max_k_p06_t350 = which(norm2_p06_t350 == max(norm2_p06_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t350

p06_bt350_adaptive = beta_est(p0_inicial06[,8], gamma_k1[max_k_p06_t350], 100)
hist(p0_inicial06[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt350_adaptive[,2], p06_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.6 ######### 
######### t = 400 ######### 
############################# 
p06_est_gamma_k1_t400 = read.table("p06_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p06_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t400[,m] - p06_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t400[m,l] = sum((p06_est_gamma_k1_t400[,m] - p06_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p06_t400[m,l] = NA
    }
    diag(norm2_p06_t400) = NA
  }
}


max_k_p06_t400 = which(norm2_p06_t400 == max(norm2_p06_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t400

p06_bt400_adaptive = beta_est(p0_inicial06[,9], gamma_k1[max_k_p06_t400], 100)
hist(p0_inicial06[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt400_adaptive[,2], p06_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.6 ######### 
######### t = 450 ######### 
############################# 
p06_est_gamma_k1_t450 = read.table("p06_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p06_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t450[,m] - p06_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t450[m,l] = sum((p06_est_gamma_k1_t450[,m] - p06_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p06_t450[m,l] = NA
    }
    diag(norm2_p06_t450) = NA
  }
}


max_k_p06_t450 = which(norm2_p06_t450 == max(norm2_p06_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t450

p06_bt450_adaptive = beta_est(p0_inicial06[,10], gamma_k1[max_k_p06_t450], 100)
hist(p0_inicial06[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt450_adaptive[,2], p06_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.6 ######### 
######### t = 500 ######### 
############################# 
p06_est_gamma_k1_t500 = read.table("p06_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p06_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p06_est_gamma_k1_t500[,m] - p06_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p06_t500[m,l] = sum((p06_est_gamma_k1_t500[,m] - p06_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p06_t500[m,l] = NA
    }
    diag(norm2_p06_t500) = NA
  }
}


max_k_p06_t500 = which(norm2_p06_t500 == max(norm2_p06_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p06_t500

p06_bt500_adaptive = beta_est(p0_inicial06[,11], gamma_k1[max_k_p06_t500], 100)
hist(p0_inicial06[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p06_bt500_adaptive[,2], p06_bt500_adaptive[,1], col = "blue", lwd=2)


###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.7 ######### 
######### t = 1 ######### 
############################# 
p07_est_gamma_k1_t1 = read.table("p07_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p07_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t1[,m] - p07_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t1[m,l] = sum((p07_est_gamma_k1_t1[,m] - p07_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p07_t1[m,l] = NA
    }
    diag(norm2_p07_t1) = NA
  }
}


max_k_p07_t1 = which(norm2_p07_t1 == max(norm2_p07_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t1

p07_bt1_adaptive = beta_est(p0_inicial07[,1], gamma_k1[max_k_p07_t1], 100)
hist(p0_inicial07[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt1_adaptive[,2], p07_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.7 ######### 
######### t = 50 ######### 
############################# 
p07_est_gamma_k1_t50 = read.table("p07_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p07_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t50[,m] - p07_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t50[m,l] = sum((p07_est_gamma_k1_t50[,m] - p07_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p07_t50[m,l] = NA
    }
    diag(norm2_p07_t50) = NA
  }
}


max_k_p07_t50 = which(norm2_p07_t50 == max(norm2_p07_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t50

p07_bt50_adaptive = beta_est(p0_inicial07[,2], gamma_k1[max_k_p07_t50], 100)
hist(p0_inicial07[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt50_adaptive[,2], p07_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.7 ######### 
######### t = 100 ######### 
############################# 
p07_est_gamma_k1_t100 = read.table("p07_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p07_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t100[,m] - p07_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t100[m,l] = sum((p07_est_gamma_k1_t100[,m] - p07_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p07_t100[m,l] = NA
    }
    diag(norm2_p07_t100) = NA
  }
}


max_k_p07_t100 = which(norm2_p07_t100 == max(norm2_p07_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t100

p07_bt100_adaptive = beta_est(p0_inicial07[,3], gamma_k1[max_k_p07_t100], 100)
hist(p0_inicial07[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt100_adaptive[,2], p07_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.7 ######### 
######### t = 150 ######### 
############################# 
p07_est_gamma_k1_t150 = read.table("p07_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p07_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t150[,m] - p07_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t150[m,l] = sum((p07_est_gamma_k1_t150[,m] - p07_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p07_t150[m,l] = NA
    }
    diag(norm2_p07_t150) = NA
  }
}


max_k_p07_t150 = which(norm2_p07_t150 == max(norm2_p07_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t150

p07_bt150_adaptive = beta_est(p0_inicial07[,4], gamma_k1[max_k_p07_t150], 100)
hist(p0_inicial07[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt150_adaptive[,2], p07_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.7 ######### 
######### t = 200 ######### 
############################# 
p07_est_gamma_k1_t200 = read.table("p07_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p07_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t200[,m] - p07_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t200[m,l] = sum((p07_est_gamma_k1_t200[,m] - p07_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p07_t200[m,l] = NA
    }
    diag(norm2_p07_t200) = NA
  }
}


max_k_p07_t200 = which(norm2_p07_t200 == max(norm2_p07_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t200

p07_bt200_adaptive = beta_est(p0_inicial07[,5], gamma_k1[max_k_p07_t200], 100)
hist(p0_inicial07[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt200_adaptive[,2], p07_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.7 ######### 
######### t = 250 ######### 
############################# 
p07_est_gamma_k1_t250 = read.table("p07_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p07_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t250[,m] - p07_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t250[m,l] = sum((p07_est_gamma_k1_t250[,m] - p07_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p07_t250[m,l] = NA
    }
    diag(norm2_p07_t250) = NA
  }
}


max_k_p07_t250 = which(norm2_p07_t250 == max(norm2_p07_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t250

p07_bt250_adaptive = beta_est(p0_inicial07[,6], gamma_k1[max_k_p07_t250], 100)
hist(p0_inicial07[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt250_adaptive[,2], p07_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.7 ######### 
######### t = 300 ######### 
############################# 
p07_est_gamma_k1_t300 = read.table("p07_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p07_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t300[,m] - p07_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t300[m,l] = sum((p07_est_gamma_k1_t300[,m] - p07_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p07_t300[m,l] = NA
    }
    diag(norm2_p07_t300) = NA
  }
}


max_k_p07_t300 = which(norm2_p07_t300 == max(norm2_p07_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t300

p07_bt300_adaptive = beta_est(p0_inicial07[,7], gamma_k1[max_k_p07_t300], 100)
hist(p0_inicial07[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt300_adaptive[,2], p07_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.7 ######### 
######### t = 350 ######### 
############################# 
p07_est_gamma_k1_t350 = read.table("p07_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p07_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t350[,m] - p07_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t350[m,l] = sum((p07_est_gamma_k1_t350[,m] - p07_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p07_t350[m,l] = NA
    }
    diag(norm2_p07_t350) = NA
  }
}


max_k_p07_t350 = which(norm2_p07_t350 == max(norm2_p07_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t350

p07_bt350_adaptive = beta_est(p0_inicial07[,8], gamma_k1[max_k_p07_t350], 100)
hist(p0_inicial07[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt350_adaptive[,2], p07_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.7 ######### 
######### t = 400 ######### 
############################# 
p07_est_gamma_k1_t400 = read.table("p07_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p07_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t400[,m] - p07_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t400[m,l] = sum((p07_est_gamma_k1_t400[,m] - p07_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p07_t400[m,l] = NA
    }
    diag(norm2_p07_t400) = NA
  }
}


max_k_p07_t400 = which(norm2_p07_t400 == max(norm2_p07_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t400

p07_bt400_adaptive = beta_est(p0_inicial07[,9], gamma_k1[max_k_p07_t400], 100)
hist(p0_inicial07[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt400_adaptive[,2], p07_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.7 ######### 
######### t = 450 ######### 
############################# 
p07_est_gamma_k1_t450 = read.table("p07_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p07_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t450[,m] - p07_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t450[m,l] = sum((p07_est_gamma_k1_t450[,m] - p07_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p07_t450[m,l] = NA
    }
    diag(norm2_p07_t450) = NA
  }
}


max_k_p07_t450 = which(norm2_p07_t450 == max(norm2_p07_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t450

p07_bt450_adaptive = beta_est(p0_inicial07[,10], gamma_k1[max_k_p07_t450], 100)
hist(p0_inicial07[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt450_adaptive[,2], p07_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.7 ######### 
######### t = 500 ######### 
############################# 
p07_est_gamma_k1_t500 = read.table("p07_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p07_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p07_est_gamma_k1_t500[,m] - p07_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p07_t500[m,l] = sum((p07_est_gamma_k1_t500[,m] - p07_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p07_t500[m,l] = NA
    }
    diag(norm2_p07_t500) = NA
  }
}


max_k_p07_t500 = which(norm2_p07_t500 == max(norm2_p07_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p07_t500

p07_bt500_adaptive = beta_est(p0_inicial07[,11], gamma_k1[max_k_p07_t500], 100)
hist(p0_inicial07[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p07_bt500_adaptive[,2], p07_bt500_adaptive[,1], col = "blue", lwd=2)

###################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.8 ######### 
######### t = 1 ######### 
############################# 
p08_est_gamma_k1_t1 = read.table("p08_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p08_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t1[,m] - p08_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t1[m,l] = sum((p08_est_gamma_k1_t1[,m] - p08_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p08_t1[m,l] = NA
    }
    diag(norm2_p08_t1) = NA
  }
}


max_k_p08_t1 = which(norm2_p08_t1 == max(norm2_p08_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t1

p08_bt1_adaptive = beta_est(p0_inicial08[,1], gamma_k1[max_k_p08_t1], 100)
hist(p0_inicial08[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt1_adaptive[,2], p08_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.8 ######### 
######### t = 50 ######### 
############################# 
p08_est_gamma_k1_t50 = read.table("p08_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p08_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t50[,m] - p08_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t50[m,l] = sum((p08_est_gamma_k1_t50[,m] - p08_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p08_t50[m,l] = NA
    }
    diag(norm2_p08_t50) = NA
  }
}


max_k_p08_t50 = which(norm2_p08_t50 == max(norm2_p08_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t50

p08_bt50_adaptive = beta_est(p0_inicial08[,2], gamma_k1[max_k_p08_t50], 100)
hist(p0_inicial08[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt50_adaptive[,2], p08_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.8 ######### 
######### t = 100 ######### 
############################# 
p08_est_gamma_k1_t100 = read.table("p08_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p08_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t100[,m] - p08_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t100[m,l] = sum((p08_est_gamma_k1_t100[,m] - p08_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p08_t100[m,l] = NA
    }
    diag(norm2_p08_t100) = NA
  }
}


max_k_p08_t100 = which(norm2_p08_t100 == max(norm2_p08_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t100

p08_bt100_adaptive = beta_est(p0_inicial08[,3], gamma_k1[max_k_p08_t100], 100)
hist(p0_inicial08[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt100_adaptive[,2], p08_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.8 ######### 
######### t = 150 ######### 
############################# 
p08_est_gamma_k1_t150 = read.table("p08_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p08_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t150[,m] - p08_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t150[m,l] = sum((p08_est_gamma_k1_t150[,m] - p08_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p08_t150[m,l] = NA
    }
    diag(norm2_p08_t150) = NA
  }
}


max_k_p08_t150 = which(norm2_p08_t150 == max(norm2_p08_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t150

p08_bt150_adaptive = beta_est(p0_inicial08[,4], gamma_k1[max_k_p08_t150], 100)
hist(p0_inicial08[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt150_adaptive[,2], p08_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.8 ######### 
######### t = 200 ######### 
############################# 
p08_est_gamma_k1_t200 = read.table("p08_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p08_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t200[,m] - p08_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t200[m,l] = sum((p08_est_gamma_k1_t200[,m] - p08_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p08_t200[m,l] = NA
    }
    diag(norm2_p08_t200) = NA
  }
}


max_k_p08_t200 = which(norm2_p08_t200 == max(norm2_p08_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t200

p08_bt200_adaptive = beta_est(p0_inicial08[,5], gamma_k1[max_k_p08_t200], 100)
hist(p0_inicial08[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt200_adaptive[,2], p08_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.8 ######### 
######### t = 250 ######### 
############################# 
p08_est_gamma_k1_t250 = read.table("p08_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p08_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t250[,m] - p08_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t250[m,l] = sum((p08_est_gamma_k1_t250[,m] - p08_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p08_t250[m,l] = NA
    }
    diag(norm2_p08_t250) = NA
  }
}


max_k_p08_t250 = which(norm2_p08_t250 == max(norm2_p08_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t250

p08_bt250_adaptive = beta_est(p0_inicial08[,6], gamma_k1[max_k_p08_t250], 100)
hist(p0_inicial08[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt250_adaptive[,2], p08_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.8 ######### 
######### t = 300 ######### 
############################# 
p08_est_gamma_k1_t300 = read.table("p08_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p08_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t300[,m] - p08_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t300[m,l] = sum((p08_est_gamma_k1_t300[,m] - p08_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p08_t300[m,l] = NA
    }
    diag(norm2_p08_t300) = NA
  }
}


max_k_p08_t300 = which(norm2_p08_t300 == max(norm2_p08_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t300

p08_bt300_adaptive = beta_est(p0_inicial08[,7], gamma_k1[max_k_p08_t300], 100)
hist(p0_inicial08[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt300_adaptive[,2], p08_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.8 ######### 
######### t = 350 ######### 
############################# 
p08_est_gamma_k1_t350 = read.table("p08_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p08_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t350[,m] - p08_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t350[m,l] = sum((p08_est_gamma_k1_t350[,m] - p08_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p08_t350[m,l] = NA
    }
    diag(norm2_p08_t350) = NA
  }
}


max_k_p08_t350 = which(norm2_p08_t350 == max(norm2_p08_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t350

p08_bt350_adaptive = beta_est(p0_inicial08[,8], gamma_k1[max_k_p08_t350], 100)
hist(p0_inicial08[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt350_adaptive[,2], p08_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.8 ######### 
######### t = 400 ######### 
############################# 
p08_est_gamma_k1_t400 = read.table("p08_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p08_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t400[,m] - p08_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t400[m,l] = sum((p08_est_gamma_k1_t400[,m] - p08_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p08_t400[m,l] = NA
    }
    diag(norm2_p08_t400) = NA
  }
}


max_k_p08_t400 = which(norm2_p08_t400 == max(norm2_p08_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t400

p08_bt400_adaptive = beta_est(p0_inicial08[,9], gamma_k1[max_k_p08_t400], 100)
hist(p0_inicial08[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt400_adaptive[,2], p08_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.8 ######### 
######### t = 450 ######### 
############################# 
p08_est_gamma_k1_t450 = read.table("p08_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p08_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t450[,m] - p08_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t450[m,l] = sum((p08_est_gamma_k1_t450[,m] - p08_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p08_t450[m,l] = NA
    }
    diag(norm2_p08_t450) = NA
  }
}


max_k_p08_t450 = which(norm2_p08_t450 == max(norm2_p08_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t450

p08_bt450_adaptive = beta_est(p0_inicial08[,10], gamma_k1[max_k_p08_t450], 100)
hist(p0_inicial08[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt450_adaptive[,2], p08_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.8 ######### 
######### t = 500 ######### 
############################# 
p08_est_gamma_k1_t500 = read.table("p08_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p08_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p08_est_gamma_k1_t500[,m] - p08_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p08_t500[m,l] = sum((p08_est_gamma_k1_t500[,m] - p08_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p08_t500[m,l] = NA
    }
    diag(norm2_p08_t500) = NA
  }
}


max_k_p08_t500 = which(norm2_p08_t500 == max(norm2_p08_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p08_t500

p08_bt500_adaptive = beta_est(p0_inicial08[,11], gamma_k1[max_k_p08_t500], 100)
hist(p0_inicial08[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p08_bt500_adaptive[,2], p08_bt500_adaptive[,1], col = "blue", lwd=2)

##################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 0.9 ######### 
######### t = 1 ######### 
############################# 
p09_est_gamma_k1_t1 = read.table("p09_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p09_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t1[,m] - p09_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t1[m,l] = sum((p09_est_gamma_k1_t1[,m] - p09_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p09_t1[m,l] = NA
    }
    diag(norm2_p09_t1) = NA
  }
}


max_k_p09_t1 = which(norm2_p09_t1 == max(norm2_p09_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t1

p09_bt1_adaptive = beta_est(p0_inicial09[,1], gamma_k1[max_k_p09_t1], 100)
hist(p0_inicial09[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt1_adaptive[,2], p09_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.9 ######### 
######### t = 50 ######### 
############################# 
p09_est_gamma_k1_t50 = read.table("p09_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p09_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t50[,m] - p09_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t50[m,l] = sum((p09_est_gamma_k1_t50[,m] - p09_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p09_t50[m,l] = NA
    }
    diag(norm2_p09_t50) = NA
  }
}


max_k_p09_t50 = which(norm2_p09_t50 == max(norm2_p09_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t50

p09_bt50_adaptive = beta_est(p0_inicial09[,2], gamma_k1[max_k_p09_t50], 100)
hist(p0_inicial09[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt50_adaptive[,2], p09_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.9 ######### 
######### t = 100 ######### 
############################# 
p09_est_gamma_k1_t100 = read.table("p09_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p09_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t100[,m] - p09_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t100[m,l] = sum((p09_est_gamma_k1_t100[,m] - p09_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p09_t100[m,l] = NA
    }
    diag(norm2_p09_t100) = NA
  }
}


max_k_p09_t100 = which(norm2_p09_t100 == max(norm2_p09_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t100

p09_bt100_adaptive = beta_est(p0_inicial09[,3], gamma_k1[max_k_p09_t100], 100)
hist(p0_inicial09[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt100_adaptive[,2], p09_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.9 ######### 
######### t = 150 ######### 
############################# 
p09_est_gamma_k1_t150 = read.table("p09_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p09_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t150[,m] - p09_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t150[m,l] = sum((p09_est_gamma_k1_t150[,m] - p09_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p09_t150[m,l] = NA
    }
    diag(norm2_p09_t150) = NA
  }
}


max_k_p09_t150 = which(norm2_p09_t150 == max(norm2_p09_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t150

p09_bt150_adaptive = beta_est(p0_inicial09[,4], gamma_k1[max_k_p09_t150], 100)
hist(p0_inicial09[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt150_adaptive[,2], p09_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.9 ######### 
######### t = 200 ######### 
############################# 
p09_est_gamma_k1_t200 = read.table("p09_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p09_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t200[,m] - p09_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t200[m,l] = sum((p09_est_gamma_k1_t200[,m] - p09_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p09_t200[m,l] = NA
    }
    diag(norm2_p09_t200) = NA
  }
}


max_k_p09_t200 = which(norm2_p09_t200 == max(norm2_p09_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t200

p09_bt200_adaptive = beta_est(p0_inicial09[,5], gamma_k1[max_k_p09_t200], 100)
hist(p0_inicial09[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt200_adaptive[,2], p09_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.9 ######### 
######### t = 250 ######### 
############################# 
p09_est_gamma_k1_t250 = read.table("p09_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p09_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t250[,m] - p09_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t250[m,l] = sum((p09_est_gamma_k1_t250[,m] - p09_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p09_t250[m,l] = NA
    }
    diag(norm2_p09_t250) = NA
  }
}


max_k_p09_t250 = which(norm2_p09_t250 == max(norm2_p09_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t250

p09_bt250_adaptive = beta_est(p0_inicial09[,6], gamma_k1[max_k_p09_t250], 100)
hist(p0_inicial09[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt250_adaptive[,2], p09_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.9 ######### 
######### t = 300 ######### 
############################# 
p09_est_gamma_k1_t300 = read.table("p09_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p09_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t300[,m] - p09_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t300[m,l] = sum((p09_est_gamma_k1_t300[,m] - p09_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p09_t300[m,l] = NA
    }
    diag(norm2_p09_t300) = NA
  }
}


max_k_p09_t300 = which(norm2_p09_t300 == max(norm2_p09_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t300

p09_bt300_adaptive = beta_est(p0_inicial09[,7], gamma_k1[max_k_p09_t300], 100)
hist(p0_inicial09[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt300_adaptive[,2], p09_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.9 ######### 
######### t = 350 ######### 
############################# 
p09_est_gamma_k1_t350 = read.table("p09_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p09_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t350[,m] - p09_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t350[m,l] = sum((p09_est_gamma_k1_t350[,m] - p09_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p09_t350[m,l] = NA
    }
    diag(norm2_p09_t350) = NA
  }
}


max_k_p09_t350 = which(norm2_p09_t350 == max(norm2_p09_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t350

p09_bt350_adaptive = beta_est(p0_inicial09[,8], gamma_k1[max_k_p09_t350], 100)
hist(p0_inicial09[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt350_adaptive[,2], p09_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.9 ######### 
######### t = 400 ######### 
############################# 
p09_est_gamma_k1_t400 = read.table("p09_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p09_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t400[,m] - p09_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t400[m,l] = sum((p09_est_gamma_k1_t400[,m] - p09_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p09_t400[m,l] = NA
    }
    diag(norm2_p09_t400) = NA
  }
}


max_k_p09_t400 = which(norm2_p09_t400 == max(norm2_p09_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t400

p09_bt400_adaptive = beta_est(p0_inicial09[,9], gamma_k1[max_k_p09_t400], 100)
hist(p0_inicial09[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt400_adaptive[,2], p09_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 0.9 ######### 
######### t = 450 ######### 
############################# 
p09_est_gamma_k1_t450 = read.table("p09_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p09_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t450[,m] - p09_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t450[m,l] = sum((p09_est_gamma_k1_t450[,m] - p09_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p09_t450[m,l] = NA
    }
    diag(norm2_p09_t450) = NA
  }
}


max_k_p09_t450 = which(norm2_p09_t450 == max(norm2_p09_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t450

p09_bt450_adaptive = beta_est(p0_inicial09[,10], gamma_k1[max_k_p09_t450], 100)
hist(p0_inicial09[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt450_adaptive[,2], p09_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 0.9 ######### 
######### t = 500 ######### 
############################# 
p09_est_gamma_k1_t500 = read.table("p09_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p09_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p09_est_gamma_k1_t500[,m] - p09_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p09_t500[m,l] = sum((p09_est_gamma_k1_t500[,m] - p09_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p09_t500[m,l] = NA
    }
    diag(norm2_p09_t500) = NA
  }
}


max_k_p09_t500 = which(norm2_p09_t500 == max(norm2_p09_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p09_t500

p09_bt500_adaptive = beta_est(p0_inicial09[,11], gamma_k1[max_k_p09_t500], 100)
hist(p0_inicial09[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p09_bt500_adaptive[,2], p09_bt500_adaptive[,1], col = "blue", lwd=2)

##################################################################################### 
###################################################################################### 
###################################################################################### 
############################# 
######### p0 = 1 ######### 
######### t = 1 ######### 
############################# 
p1_est_gamma_k1_t1 = read.table("p1_est_gamma_k1_t1.txt", header = T, sep = ' ')

k=84
norm2_p1_t1 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t1[,m] - p1_est_gamma_k1_t1[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t1[m,l] = sum((p1_est_gamma_k1_t1[,m] - p1_est_gamma_k1_t1[,l])^2)
    }else{
      norm2_p1_t1[m,l] = NA
    }
    diag(norm2_p1_t1) = NA
  }
}


max_k_p1_t1 = which(norm2_p1_t1 == max(norm2_p1_t1, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t1

p1_bt1_adaptive = beta_est(p0_inicial1[,1], gamma_k1[max_k_p1_t1], 100)
hist(p0_inicial1[,1], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt1_adaptive[,2], p1_bt1_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 1 ######### 
######### t = 50 ######### 
############################# 
p1_est_gamma_k1_t50 = read.table("p1_est_gamma_k1_t50.txt", header = T, sep = ' ')

k=84
norm2_p1_t50 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t50[,m] - p1_est_gamma_k1_t50[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t50[m,l] = sum((p1_est_gamma_k1_t50[,m] - p1_est_gamma_k1_t50[,l])^2)
    }else{
      norm2_p1_t50[m,l] = NA
    }
    diag(norm2_p1_t50) = NA
  }
}


max_k_p1_t50 = which(norm2_p1_t50 == max(norm2_p1_t50, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t50

p1_bt50_adaptive = beta_est(p0_inicial1[,2], gamma_k1[max_k_p1_t50], 100)
hist(p0_inicial1[,2], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt50_adaptive[,2], p1_bt50_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 1 ######### 
######### t = 100 ######### 
############################# 
p1_est_gamma_k1_t100 = read.table("p1_est_gamma_k1_t100.txt", header = T, sep = ' ')

k=84
norm2_p1_t100 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t100[,m] - p1_est_gamma_k1_t100[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t100[m,l] = sum((p1_est_gamma_k1_t100[,m] - p1_est_gamma_k1_t100[,l])^2)
    }else{
      norm2_p1_t100[m,l] = NA
    }
    diag(norm2_p1_t100) = NA
  }
}


max_k_p1_t100 = which(norm2_p1_t100 == max(norm2_p1_t100, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t100

p1_bt100_adaptive = beta_est(p0_inicial1[,3], gamma_k1[max_k_p1_t100], 100)
hist(p0_inicial1[,3], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt100_adaptive[,2], p1_bt100_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 1 ######### 
######### t = 150 ######### 
############################# 
p1_est_gamma_k1_t150 = read.table("p1_est_gamma_k1_t150.txt", header = T, sep = ' ')

k=84
norm2_p1_t150 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t150[,m] - p1_est_gamma_k1_t150[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t150[m,l] = sum((p1_est_gamma_k1_t150[,m] - p1_est_gamma_k1_t150[,l])^2)
    }else{
      norm2_p1_t150[m,l] = NA
    }
    diag(norm2_p1_t150) = NA
  }
}


max_k_p1_t150 = which(norm2_p1_t150 == max(norm2_p1_t150, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t150

p1_bt150_adaptive = beta_est(p0_inicial1[,4], gamma_k1[max_k_p1_t150], 100)
hist(p0_inicial1[,4], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt150_adaptive[,2], p1_bt150_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 1 ######### 
######### t = 200 ######### 
############################# 
p1_est_gamma_k1_t200 = read.table("p1_est_gamma_k1_t200.txt", header = T, sep = ' ')

k=84
norm2_p1_t200 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t200[,m] - p1_est_gamma_k1_t200[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t200[m,l] = sum((p1_est_gamma_k1_t200[,m] - p1_est_gamma_k1_t200[,l])^2)
    }else{
      norm2_p1_t200[m,l] = NA
    }
    diag(norm2_p1_t200) = NA
  }
}


max_k_p1_t200 = which(norm2_p1_t200 == max(norm2_p1_t200, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t200

p1_bt200_adaptive = beta_est(p0_inicial1[,5], gamma_k1[max_k_p1_t200], 100)
hist(p0_inicial1[,5], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt200_adaptive[,2], p1_bt200_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 1 ######### 
######### t = 250 ######### 
############################# 
p1_est_gamma_k1_t250 = read.table("p1_est_gamma_k1_t250.txt", header = T, sep = ' ')

k=84
norm2_p1_t250 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t250[,m] - p1_est_gamma_k1_t250[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t250[m,l] = sum((p1_est_gamma_k1_t250[,m] - p1_est_gamma_k1_t250[,l])^2)
    }else{
      norm2_p1_t250[m,l] = NA
    }
    diag(norm2_p1_t250) = NA
  }
}


max_k_p1_t250 = which(norm2_p1_t250 == max(norm2_p1_t250, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t250

p1_bt250_adaptive = beta_est(p0_inicial1[,6], gamma_k1[max_k_p1_t250], 100)
hist(p0_inicial1[,6], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt250_adaptive[,2], p1_bt250_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 1 ######### 
######### t = 300 ######### 
############################# 
p1_est_gamma_k1_t300 = read.table("p1_est_gamma_k1_t300.txt", header = T, sep = ' ')

k=84
norm2_p1_t300 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t300[,m] - p1_est_gamma_k1_t300[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t300[m,l] = sum((p1_est_gamma_k1_t300[,m] - p1_est_gamma_k1_t300[,l])^2)
    }else{
      norm2_p1_t300[m,l] = NA
    }
    diag(norm2_p1_t300) = NA
  }
}


max_k_p1_t300 = which(norm2_p1_t300 == max(norm2_p1_t300, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t300

p1_bt300_adaptive = beta_est(p0_inicial1[,7], gamma_k1[max_k_p1_t300], 100)
hist(p0_inicial1[,7], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt300_adaptive[,2], p1_bt300_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 1 ######### 
######### t = 350 ######### 
############################# 
p1_est_gamma_k1_t350 = read.table("p1_est_gamma_k1_t350.txt", header = T, sep = ' ')

k=84
norm2_p1_t350 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t350[,m] - p1_est_gamma_k1_t350[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t350[m,l] = sum((p1_est_gamma_k1_t350[,m] - p1_est_gamma_k1_t350[,l])^2)
    }else{
      norm2_p1_t350[m,l] = NA
    }
    diag(norm2_p1_t350) = NA
  }
}


max_k_p1_t350 = which(norm2_p1_t350 == max(norm2_p1_t350, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t350

p1_bt350_adaptive = beta_est(p0_inicial1[,8], gamma_k1[max_k_p1_t350], 100)
hist(p0_inicial1[,8], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt350_adaptive[,2], p1_bt350_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 1 ######### 
######### t = 400 ######### 
############################# 
p1_est_gamma_k1_t400 = read.table("p1_est_gamma_k1_t400.txt", header = T, sep = ' ')

k=84
norm2_p1_t400 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t400[,m] - p1_est_gamma_k1_t400[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t400[m,l] = sum((p1_est_gamma_k1_t400[,m] - p1_est_gamma_k1_t400[,l])^2)
    }else{
      norm2_p1_t400[m,l] = NA
    }
    diag(norm2_p1_t400) = NA
  }
}


max_k_p1_t400 = which(norm2_p1_t400 == max(norm2_p1_t400, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t400

p1_bt400_adaptive = beta_est(p0_inicial1[,9], gamma_k1[max_k_p1_t400], 100)
hist(p0_inicial1[,9], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt400_adaptive[,2], p1_bt400_adaptive[,1], col = "blue", lwd=2)

############################# 
######### p0 = 1 ######### 
######### t = 450 ######### 
############################# 
p1_est_gamma_k1_t450 = read.table("p1_est_gamma_k1_t450.txt", header = T, sep = ' ')

k=84
norm2_p1_t450 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t450[,m] - p1_est_gamma_k1_t450[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t450[m,l] = sum((p1_est_gamma_k1_t450[,m] - p1_est_gamma_k1_t450[,l])^2)
    }else{
      norm2_p1_t450[m,l] = NA
    }
    diag(norm2_p1_t450) = NA
  }
}


max_k_p1_t450 = which(norm2_p1_t450 == max(norm2_p1_t450, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t450

p1_bt450_adaptive = beta_est(p0_inicial1[,10], gamma_k1[max_k_p1_t450], 100)
hist(p0_inicial1[,10], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt450_adaptive[,2], p1_bt450_adaptive[,1], col = "blue", lwd=2)


############################# 
######### p0 = 1 ######### 
######### t = 500 ######### 
############################# 
p1_est_gamma_k1_t500 = read.table("p1_est_gamma_k1_t500.txt", header = T, sep = ' ')

k=84
norm2_p1_t500 = matrix(NA, nco = k, nrow = k)
for (l in 1:k) {
  for (m in 1:k) {
    if((m <= l) & sum((p1_est_gamma_k1_t500[,m] - p1_est_gamma_k1_t500[,l])^2) <= 2*psi_gamma[m]){
      norm2_p1_t500[m,l] = sum((p1_est_gamma_k1_t500[,m] - p1_est_gamma_k1_t500[,l])^2)
    }else{
      norm2_p1_t500[m,l] = NA
    }
    diag(norm2_p1_t500) = NA
  }
}


max_k_p1_t500 = which(norm2_p1_t500 == max(norm2_p1_t500, na.rm = T), arr.ind=T)[1,2]
max_k_p1_t500

p1_bt500_adaptive = beta_est(p0_inicial1[,11], gamma_k1[max_k_p1_t500], 100)
hist(p0_inicial1[,11], breaks = 45, freq = F, ylim = c(0,7), xlim = c(0,1) ,xlab = expression(p[0]), main = "")
lines(p1_bt500_adaptive[,2], p1_bt500_adaptive[,1], col = "blue", lwd=2)

#############
## hay que construir los heatmaps para ver las distancias entre las 
## filas son las generaciones
## columnas con los valores de p0 inicial 
## cada celda es la distancia L2

#1. asymptotical expansion (11x9)
t = c(1,50,100,150,200,250,300,350,400,450,500)
x0 = seq(0,1, 1/10) #valores iniciales

p0_beta = cbind(p0_bt1_adaptive[,1], p0_bt50_adaptive[,1], p0_bt100_adaptive[,1],
                 p0_bt150_adaptive[,1], p0_bt200_adaptive[,1], p0_bt250_adaptive[,1],
                 p0_bt300_adaptive[,1], p0_bt350_adaptive[,1], p0_bt400_adaptive[,1],
                 p0_bt450_adaptive[,1], p0_bt500_adaptive[,1])

p01_beta = cbind(p01_bt1_adaptive[,1], p01_bt50_adaptive[,1], p01_bt100_adaptive[,1],
                 p01_bt150_adaptive[,1], p01_bt200_adaptive[,1], p01_bt250_adaptive[,1],
                 p01_bt300_adaptive[,1], p01_bt350_adaptive[,1], p01_bt400_adaptive[,1],
                 p01_bt450_adaptive[,1], p01_bt500_adaptive[,1])

p02_beta = cbind(p02_bt1_adaptive[,1], p02_bt50_adaptive[,1], p02_bt100_adaptive[,1],
                 p02_bt150_adaptive[,1], p02_bt200_adaptive[,1], p02_bt250_adaptive[,1],
                 p02_bt300_adaptive[,1], p02_bt350_adaptive[,1], p02_bt400_adaptive[,1],
                 p02_bt450_adaptive[,1], p02_bt500_adaptive[,1])

p03_beta = cbind(p03_bt1_adaptive[,1], p03_bt50_adaptive[,1], p03_bt100_adaptive[,1],
                 p03_bt150_adaptive[,1], p03_bt200_adaptive[,1], p03_bt250_adaptive[,1],
                 p03_bt300_adaptive[,1], p03_bt350_adaptive[,1], p03_bt400_adaptive[,1],
                 p03_bt450_adaptive[,1], p03_bt500_adaptive[,1])

p04_beta = cbind(p04_bt1_adaptive[,1], p04_bt50_adaptive[,1], p04_bt100_adaptive[,1],
                 p04_bt150_adaptive[,1], p04_bt200_adaptive[,1], p04_bt250_adaptive[,1],
                 p04_bt300_adaptive[,1], p04_bt350_adaptive[,1], p04_bt400_adaptive[,1],
                 p04_bt450_adaptive[,1], p04_bt500_adaptive[,1])

p05_beta = cbind(p05_bt1_adaptive[,1], p05_bt50_adaptive[,1], p05_bt100_adaptive[,1],
                 p05_bt150_adaptive[,1], p05_bt200_adaptive[,1], p05_bt250_adaptive[,1],
                 p05_bt300_adaptive[,1], p05_bt350_adaptive[,1], p05_bt400_adaptive[,1],
                 p05_bt450_adaptive[,1], p05_bt500_adaptive[,1])

p06_beta = cbind(p06_bt1_adaptive[,1], p06_bt50_adaptive[,1], p06_bt100_adaptive[,1],
                 p06_bt150_adaptive[,1], p06_bt200_adaptive[,1], p06_bt250_adaptive[,1],
                 p06_bt300_adaptive[,1], p06_bt350_adaptive[,1], p06_bt400_adaptive[,1],
                 p06_bt450_adaptive[,1], p06_bt500_adaptive[,1])

p07_beta = cbind(p07_bt1_adaptive[,1], p07_bt50_adaptive[,1], p07_bt100_adaptive[,1],
                 p07_bt150_adaptive[,1], p07_bt200_adaptive[,1], p07_bt250_adaptive[,1],
                 p07_bt300_adaptive[,1], p07_bt350_adaptive[,1], p07_bt400_adaptive[,1],
                 p07_bt450_adaptive[,1], p07_bt500_adaptive[,1])

p08_beta = cbind(p08_bt1_adaptive[,1], p08_bt50_adaptive[,1], p08_bt100_adaptive[,1],
                 p08_bt150_adaptive[,1], p08_bt200_adaptive[,1], p08_bt250_adaptive[,1],
                 p08_bt300_adaptive[,1], p08_bt350_adaptive[,1], p08_bt400_adaptive[,1],
                 p08_bt450_adaptive[,1], p08_bt500_adaptive[,1])


p09_beta = cbind(p09_bt1_adaptive[,1], p09_bt50_adaptive[,1], p09_bt100_adaptive[,1],
                 p09_bt150_adaptive[,1], p09_bt200_adaptive[,1], p09_bt250_adaptive[,1],
                 p09_bt300_adaptive[,1], p09_bt350_adaptive[,1], p09_bt400_adaptive[,1],
                 p09_bt450_adaptive[,1], p09_bt500_adaptive[,1])

p1_beta = cbind(p1_bt1_adaptive[,1], p1_bt50_adaptive[,1], p1_bt100_adaptive[,1],
                 p1_bt150_adaptive[,1], p1_bt200_adaptive[,1], p1_bt250_adaptive[,1],
                 p1_bt300_adaptive[,1], p1_bt350_adaptive[,1], p1_bt400_adaptive[,1],
                 p1_bt450_adaptive[,1], p1_bt500_adaptive[,1])

heat_expan = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
 heat_expan[i,1] = sqrt(sum((p0_adaptive[,i] - rep(0, 100))^2))
 heat_expan[i,2] = sqrt(sum((p01_adaptive[,i] - expan(x0[2], 100, t[i]))^2))
 heat_expan[i,3] = sqrt(sum((p02_adaptive[,i] - expan(x0[3], 100, t[i]))^2))
 heat_expan[i,4] = sqrt(sum((p03_adaptive[,i] - expan(x0[4], 100, t[i]))^2))
 heat_expan[i,5] = sqrt(sum((p04_adaptive[,i] - expan(x0[5], 100, t[i]))^2))
 heat_expan[i,6] = sqrt(sum((p05_adaptive[,i] - expan(x0[6], 100, t[i]))^2))
 heat_expan[i,7] = sqrt(sum((p06_adaptive[,i] - expan(x0[7], 100, t[i]))^2))
 heat_expan[i,8] = sqrt(sum((p07_adaptive[,i] - expan(x0[8], 100, t[i]))^2))
 heat_expan[i,9] = sqrt(sum((p08_adaptive[,i] - expan(x0[9], 100, t[i]))^2))
 heat_expan[i,10] = sqrt(sum((p09_adaptive[,i] - expan(x0[10], 100, t[i]))^2))
 heat_expan[i,11] = sqrt(sum((p1_adaptive[,i] - rep(0, 100))^2))
}

colnames(heat_expan) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_expan) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_expan

heat_gauss = matrix(NA, ncol = 11, nrow = 11)

for (i in 1:11) {
  heat_gauss[i,1] = sqrt(sum((p0_adaptive[,i] - rep(0, 100))^2))
  heat_gauss[i,2] = sqrt(sum((p01_adaptive[,i] - gauss(x0[2], 100, t[i]))^2))
  heat_gauss[i,3] = sqrt(sum((p02_adaptive[,i] - gauss(x0[3], 100, t[i]))^2))
  heat_gauss[i,4] = sqrt(sum((p03_adaptive[,i] - gauss(x0[4], 100, t[i]))^2))
  heat_gauss[i,5] = sqrt(sum((p04_adaptive[,i] - gauss(x0[5], 100, t[i]))^2))
  heat_gauss[i,6] = sqrt(sum((p05_adaptive[,i] - gauss(x0[6], 100, t[i]))^2))
  heat_gauss[i,7] = sqrt(sum((p06_adaptive[,i] - gauss(x0[7], 100, t[i]))^2))
  heat_gauss[i,8] = sqrt(sum((p07_adaptive[,i] - gauss(x0[8], 100, t[i]))^2))
  heat_gauss[i,9] = sqrt(sum((p08_adaptive[,i] - gauss(x0[9], 100, t[i]))^2))
  heat_gauss[i,10] = sqrt(sum((p09_adaptive[,i] - gauss(x0[10], 100, t[i]))^2))
  heat_gauss[i,11] = sqrt(sum((p1_adaptive[,i] - rep(0, 100))^2))
}

colnames(heat_gauss) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_gauss) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_gauss

library(RColorBrewer)

heatmap(heat_expan, Colv = NA, Rowv = NA, scale = 'column', xlab = expression(p[0]),
        ylab = "t", main = 'Asymptotical expansion', 
        col= colorRampPalette(brewer.pal(8, "Oranges"))(25))


##beta parametrizada como media y varianza
media = matrix(NA, ncol= 11, nrow = 11)
for (i in 1:11) {
  media[i,1] = mean(p0_inicial0[,i])
  media[i,2] = mean(p0_inicial01[,i])
  media[i,3] = mean(p0_inicial02[,i])
  media[i,4] = mean(p0_inicial03[,i])
  media[i,5] = mean(p0_inicial04[,i])
  media[i,6] = mean(p0_inicial05[,i])
  media[i,7] = mean(p0_inicial06[,i])
  media[i,8] = mean(p0_inicial07[,i])
  media[i,9] = mean(p0_inicial08[,i])
  media[i,10] = mean(p0_inicial09[,i])
  media[i,11] = mean(p0_inicial1[,i])
}

colnames(media) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(media) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
media


variance = matrix(NA, ncol= 11, nrow = 11)
for (i in 1:11) {
  variance[i,1] = var(p0_inicial0[,i])
  variance[i,2] = var(p0_inicial01[,i])
  variance[i,3] = var(p0_inicial02[,i])
  variance[i,4] = var(p0_inicial03[,i])
  variance[i,5] = var(p0_inicial04[,i])
  variance[i,6] = var(p0_inicial05[,i])
  variance[i,7] = var(p0_inicial06[,i])
  variance[i,8] = var(p0_inicial07[,i])
  variance[i,9] = var(p0_inicial08[,i])
  variance[i,10] = var(p0_inicial09[,i])
  variance[i,11] = var(p0_inicial1[,i])
}

for (i in 1:11) {
  for (j in 1:11){
    if (variance[i,j] == 0){variance[i,j] = min( variance[variance!=min(variance)] )} else
    {variance[i,j] = variance[i,j]}
  }
}

#cambia los 0 por el segundo más chico
colnames(variance) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(variance) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
variance

##parametros de forma en función a la media y la varianza
sh1 = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  for (j in 1:11){
    sh1[i,j] = abs((((media[i,j] * (1 - media[i,j])) / variance[i,j]) - 1) * media[i,j])
  }
}

colnames(sh1) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(sh1) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
sh1

sh2 = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  for (j in 1:11){
    sh2[i,j] = abs(((media[i,j] * (1 - media[i,j])) / variance[i,j]) - 1) * (1 - media[i,j])
  }
}

colnames(sh2) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(sh2) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
sh2


##probabilidades p0

prob_p0 = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  prob_p0[i,1] = sum((p0_inicial0[,i] == 0)*1) / length(p0_inicial0[,i])
  prob_p0[i,2] = sum((p0_inicial01[,i] == 0)*1) / length(p0_inicial01[,i])
  prob_p0[i,3] = sum((p0_inicial02[,i] == 0)*1) / length(p0_inicial02[,i])
  prob_p0[i,4] = sum((p0_inicial03[,i] == 0)*1) / length(p0_inicial03[,i])
  prob_p0[i,5] = sum((p0_inicial04[,i] == 0)*1) / length(p0_inicial04[,i])
  prob_p0[i,6] = sum((p0_inicial05[,i] == 0)*1) / length(p0_inicial05[,i])
  prob_p0[i,7] = sum((p0_inicial06[,i] == 0)*1) / length(p0_inicial06[,i])
  prob_p0[i,8] = sum((p0_inicial07[,i] == 0)*1) / length(p0_inicial07[,i])
  prob_p0[i,9] = sum((p0_inicial08[,i] == 0)*1) / length(p0_inicial08[,i])
  prob_p0[i,10] = sum((p0_inicial09[,i] == 0)*1) / length(p0_inicial09[,i])
  prob_p0[i,11] = sum((p0_inicial1[,i] == 0)*1) / length(p0_inicial1[,i])
}

colnames(prob_p0) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(prob_p0) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
prob_p0

prob_p1 = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  prob_p1[i,1] = sum((p0_inicial0[,i] == 1)*1) / length(p0_inicial0[,i])
  prob_p1[i,2] = sum((p0_inicial01[,i] == 1)*1) / length(p0_inicial01[,i])
  prob_p1[i,3] = sum((p0_inicial02[,i] == 1)*1) / length(p0_inicial02[,i])
  prob_p1[i,4] = sum((p0_inicial03[,i] == 1)*1) / length(p0_inicial03[,i])
  prob_p1[i,5] = sum((p0_inicial04[,i] == 1)*1) / length(p0_inicial04[,i])
  prob_p1[i,6] = sum((p0_inicial05[,i] == 1)*1) / length(p0_inicial05[,i])
  prob_p1[i,7] = sum((p0_inicial06[,i] == 1)*1) / length(p0_inicial06[,i])
  prob_p1[i,8] = sum((p0_inicial07[,i] == 1)*1) / length(p0_inicial07[,i])
  prob_p1[i,9] = sum((p0_inicial08[,i] == 1)*1) / length(p0_inicial08[,i])
  prob_p1[i,10] = sum((p0_inicial09[,i] == 1)*1) / length(p0_inicial09[,i])
  prob_p1[i,11] = sum((p0_inicial1[,i] == 1)*1) / length(p0_inicial1[,i])
}

colnames(prob_p1) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(prob_p1) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
prob_p1


par(mfrow = c(1,1))
x = seq(0,1, length.out = 100)
x[1] = x[100] = 0
#x = x[2:101]


heat_beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_beta[i,1] = sqrt(sum((p0_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,1], shape2 = sh2[i,1]), 0))^2))
  heat_beta[i,2] = sqrt(sum((p01_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,2], shape2 = sh2[i,2]), 0))^2))
  heat_beta[i,3] = sqrt(sum((p02_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,3], shape2 = sh2[i,3]), 0))^2))
  heat_beta[i,4] = sqrt(sum((p03_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,4], shape2 = sh2[i,4]), 0))^2))
  heat_beta[i,5] = sqrt(sum((p04_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,5], shape2 = sh2[i,5]), 0))^2))
  heat_beta[i,6] = sqrt(sum((p05_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,6], shape2 = sh2[i,6]), 0))^2))
  heat_beta[i,7] = sqrt(sum((p06_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,7], shape2 = sh2[i,7]), 0))^2))
  heat_beta[i,8] = sqrt(sum((p07_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,8], shape2 = sh2[i,8]), 0))^2))
  heat_beta[i,9] = sqrt(sum((p08_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,9], shape2 = sh2[i,9]), 0))^2))
  heat_beta[i,10] = sqrt(sum((p09_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,10], shape2 = sh2[i,10]), 0))^2))
  heat_beta[i,11] = sqrt(sum((p1_adaptive[,i] - c(0, dbeta(x[2:99], shape1 = sh1[i,11], shape2 = sh2[i,11]), 0))^2))
}

colnames(heat_beta) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_beta) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_beta



##beta con spikes

## step funciton

step_fun = function(t, t0){
  u = rep(NA, length(t))
  for (i in 1: length(t)) {
    if (t[i] < t0){
      u[i] = 0
    } else {u[i] = 1}
  }
  return(u)
}


t = seq(0,1, 1/100)
plot(step_fun(t, 0.5))
plot(step_fun(t, 0.5), type = 'l')

step_fun(t, 0.3) - step_fun(t, 0.32)
plot(t, step_fun(t, 0.3) - step_fun(t, 0.32))
plot(t, step_fun(t, 0.3) - step_fun(t, 0.32), type = 'l')


plot(ecdf(step_fun(t, 0.3) - step_fun(t, 0.32)))

##delta dirac usando la función step
#valor de epsilon 
eps = 0.001

step_fun(t, (0.3-(eps/2)))
step_fun(t, (0.3+(eps/2)))

#cómo se ve la función delta dirac
plot(t, (1/eps) * (step_fun(t, (0.01-(eps/2))) - step_fun(t, (0.01+(eps/2)))), type = 'l')

##delta dirac como función pulso rectangular
delta_dir = function(t, t0, eps){
  del_dir = (1/eps) * (step_fun(t, (t0-(eps/2))) - step_fun(t, (t0+(eps/2))))

  return(del_dir)
}

plot(delta_dir(t, 0.99, 0.001), type = 'l')
plot(density(delta_dir(t, 0.99, 0.0001)))

#viendo la función de distribución asociada
plot(ecdf(delta_dir(t, 0.99, 0.001)))

#parámetros de la beta con spikes m* y v*
m_star = (media - prob_p1) / (1 - prob_p0 - prob_p1)
m_star[,1] = 0
m_star[,11] = 1
v_star = abs((variance - media^2 - prob_p1) / (1 - prob_p0 - prob_p1) * m_star^2)
v_star[,c(1,11)] = min(v_star[,2:10])

#quito la primera y última columna que tienen NaN
a_star = (m_star*(1-m_star)/v_star)*m_star
a_star = ifelse(a_star == 0, 0.07011922, a_star)
b_star = (m_star*(1-m_star)/v_star)*(1-m_star)
b_star = ifelse(b_star == 0, 1.376186e-02, b_star)

#esto es una densidad que va cambiando para cada valores de (t,p0) (11x11)
#11x11x100 11 (valores iniciales de p0) matrices de dimension 100x11 (obs x gen)

##p0
x = seq(0, 1, length.out=102)
beta_spikes_p0 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p0[,j] = prob_p0[j,1]*delta_dir(x, 0.01, 0.001) + prob_p1[j,1]*delta_dir(x, 0.99, 0.001) +
     dbeta(x,shape1 = a_star[j,1], shape2 = b_star[j,1])
}

beta_spikes_p0 = beta_spikes_p0[2:101,]

##p=0.1
beta_spikes_p01 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p01[,j] = prob_p0[j,2]*delta_dir(x, 0.01, 0.001) + prob_p1[j,2]*delta_dir(x, 0.99, 0.001) +
     dbeta(x,shape1 = a_star[j,2], shape2 = b_star[j,2])
}

beta_spikes_p01 = beta_spikes_p01[2:101,]

##p=0.2
beta_spikes_p02 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p02[,j] = prob_p0[j,3]*delta_dir(x, 0.01, 0.001) + prob_p1[j,3]*delta_dir(x, 0.99, 0.001) +
    dbeta(x,shape1 = a_star[j,3], shape2 = b_star[j,3])
}

beta_spikes_p02 = beta_spikes_p02[2:101,]


##p=0.3
beta_spikes_p03 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p03[,j] = prob_p0[j,4]*delta_dir(x, 0.01, 0.001) + prob_p1[j,4]*delta_dir(x, 0.99, 0.001) +
    dbeta(x,shape1 = a_star[j,4], shape2 = b_star[j,4])
}

beta_spikes_p03 = beta_spikes_p03[2:101,]


##p=0.4
beta_spikes_p04 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p04[,j] = prob_p0[j,5]*delta_dir(x, 0.01, 0.001) + prob_p1[j,5]*delta_dir(x, 0.99, 0.001) +
   dbeta(x,shape1 = a_star[j,5], shape2 = b_star[j,5])
}

beta_spikes_p04 = beta_spikes_p04[2:101,]

##p=0.5
beta_spikes_p05 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p05[,j] = prob_p0[j,6]*delta_dir(x, 0.01, 0.001) + prob_p1[j,6]*delta_dir(x, 0.99, 0.001) +
    dbeta(x,shape1 = a_star[j,6], shape2 = b_star[j,6])
}

beta_spikes_p05 = beta_spikes_p05[2:101,]

##p=0.6
beta_spikes_p06 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p06[,j] = prob_p0[j,7]*delta_dir(x, 0.01, 0.001) + prob_p1[1,7]*delta_dir(x, 0.99, 0.001) +
   dbeta(x,shape1 = a_star[j,7], shape2 = b_star[j,7])
}

beta_spikes_p06 = beta_spikes_p06[2:101,]

##p=0.7
beta_spikes_p07 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p07[,j] = prob_p0[j,8]*delta_dir(x, 0.01, 0.001) + prob_p1[j,8]*delta_dir(x, 0.99, 0.001) +
    dbeta(x,shape1 = a_star[j,8], shape2 = b_star[j,8])
}

beta_spikes_p07 = beta_spikes_p07[2:101,]


##p=0.8
beta_spikes_p08 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p08[,j] = prob_p0[j,9]*delta_dir(x, 0.01, 0.001) + prob_p1[j,9]*delta_dir(x, 0.99, 0.001) +
   dbeta(x,shape1 = a_star[j,9], shape2 = b_star[j,9])
}

beta_spikes_p08 = beta_spikes_p08[2:101,]

##p=0.9
beta_spikes_p09 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p09[,j] = prob_p0[j,10]*delta_dir(x, 0.01, 0.001) + prob_p1[j,10]*delta_dir(x, 0.99, 0.001) +
  dbeta(x,shape1 = a_star[j,10], shape2 = b_star[j,10])
}

beta_spikes_p09 = beta_spikes_p09[2:101,]

##p=1
beta_spikes_p1 = matrix(0, ncol = 11, nrow = 102)

for (j in 1:11){
  beta_spikes_p1[,j] = prob_p0[j,11]*delta_dir(x, 0.01, 0.001) + prob_p1[j,11]*delta_dir(x, 0.99, 0.001) +
   dbeta(x,shape1 = a_star[j,11], shape2 = b_star[j,11])
}

beta_spikes_p1 = beta_spikes_p1[2:101,]


### heatmpap beta with spikes
heat_beta_sp = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_beta_sp[i,1] = sqrt(sum((p0_beta[,i] - beta_spikes_p0[,i])^2))
  heat_beta_sp[i,2] = sqrt(sum((p01_beta[,i] - beta_spikes_p01[,i])^2))
  heat_beta_sp[i,3] = sqrt(sum((p02_beta[,i] - beta_spikes_p02[,i])^2))
  heat_beta_sp[i,4] = sqrt(sum((p03_beta[,i] - beta_spikes_p03[,i])^2))
  heat_beta_sp[i,5] = sqrt(sum((p04_beta[,i] - beta_spikes_p04[,i])^2))
  heat_beta_sp[i,6] = sqrt(sum((p05_beta[,i] - beta_spikes_p05[,i])^2))
  heat_beta_sp[i,7] = sqrt(sum((p06_beta[,i] - beta_spikes_p06[,i])^2))
  heat_beta_sp[i,8] = sqrt(sum((p07_beta[,i] - beta_spikes_p07[,i])^2))
  heat_beta_sp[i,9] = sqrt(sum((p08_beta[,i] - beta_spikes_p08[,i])^2))
  heat_beta_sp[i,10] = sqrt(sum((p09_beta[,i] - beta_spikes_p09[,i])^2))
  heat_beta_sp[i,11] = sqrt(sum((p1_beta[,i] - beta_spikes_p1[,i])^2))
}

colnames(heat_beta_sp) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_beta_sp) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_beta_sp


###heatmap Adaptive Density Estimator
p0_adaptive = cbind(p0_bt1_adaptive[,1], p0_bt50_adaptive[,1], p0_bt100_adaptive[,1],
                     p0_bt150_adaptive[,1], p0_bt200_adaptive[,1], p0_bt250_adaptive[,1],
                     p0_bt300_adaptive[,1], p0_bt350_adaptive[,1], p0_bt400_adaptive[,1],
                     p0_bt450_adaptive[,1], p0_bt500_adaptive[,1])
p0_adaptive[1,] = 0

p01_adaptive = cbind(p01_bt1_adaptive[,1], p01_bt50_adaptive[,1], p01_bt100_adaptive[,1],
                     p01_bt150_adaptive[,1], p01_bt200_adaptive[,1], p01_bt250_adaptive[,1],
                     p01_bt300_adaptive[,1], p01_bt350_adaptive[,1], p01_bt400_adaptive[,1],
                     p01_bt450_adaptive[,1], p01_bt500_adaptive[,1])

p02_adaptive = cbind(p02_bt1_adaptive[,1], p02_bt50_adaptive[,1], p02_bt100_adaptive[,1],
                     p02_bt150_adaptive[,1], p02_bt200_adaptive[,1], p02_bt250_adaptive[,1],
                     p02_bt300_adaptive[,1], p02_bt350_adaptive[,1], p02_bt400_adaptive[,1],
                     p02_bt450_adaptive[,1], p02_bt500_adaptive[,1])

p03_adaptive = cbind(p03_bt1_adaptive[,1], p03_bt50_adaptive[,1], p03_bt100_adaptive[,1],
                     p03_bt150_adaptive[,1], p03_bt200_adaptive[,1], p03_bt250_adaptive[,1],
                     p03_bt300_adaptive[,1], p03_bt350_adaptive[,1], p03_bt400_adaptive[,1],
                     p03_bt450_adaptive[,1], p03_bt500_adaptive[,1])

p04_adaptive = cbind(p04_bt1_adaptive[,1], p04_bt50_adaptive[,1], p04_bt100_adaptive[,1],
                     p04_bt150_adaptive[,1], p04_bt200_adaptive[,1], p04_bt250_adaptive[,1],
                     p04_bt300_adaptive[,1], p04_bt350_adaptive[,1], p04_bt400_adaptive[,1],
                     p04_bt450_adaptive[,1], p04_bt500_adaptive[,1])

p05_adaptive = cbind(p05_bt1_adaptive[,1], p05_bt50_adaptive[,1], p05_bt100_adaptive[,1],
                     p05_bt150_adaptive[,1], p05_bt200_adaptive[,1], p05_bt250_adaptive[,1],
                     p05_bt300_adaptive[,1], p05_bt350_adaptive[,1], p05_bt400_adaptive[,1],
                     p05_bt450_adaptive[,1], p05_bt500_adaptive[,1])

p06_adaptive = cbind(p06_bt1_adaptive[,1], p06_bt50_adaptive[,1], p06_bt100_adaptive[,1],
                     p06_bt150_adaptive[,1], p06_bt200_adaptive[,1], p06_bt250_adaptive[,1],
                     p06_bt300_adaptive[,1], p06_bt350_adaptive[,1], p06_bt400_adaptive[,1],
                     p06_bt450_adaptive[,1], p06_bt500_adaptive[,1])

p07_adaptive = cbind(p07_bt1_adaptive[,1], p07_bt50_adaptive[,1], p07_bt100_adaptive[,1],
                     p07_bt150_adaptive[,1], p07_bt200_adaptive[,1], p07_bt250_adaptive[,1],
                     p07_bt300_adaptive[,1], p07_bt350_adaptive[,1], p07_bt400_adaptive[,1],
                     p07_bt450_adaptive[,1], p07_bt500_adaptive[,1])

p08_adaptive = cbind(p08_bt1_adaptive[,1], p08_bt50_adaptive[,1], p08_bt100_adaptive[,1],
                     p08_bt150_adaptive[,1], p08_bt200_adaptive[,1], p08_bt250_adaptive[,1],
                     p08_bt300_adaptive[,1], p08_bt350_adaptive[,1], p08_bt400_adaptive[,1],
                     p08_bt450_adaptive[,1], p08_bt500_adaptive[,1])

p09_adaptive = cbind(p09_bt1_adaptive[,1], p09_bt50_adaptive[,1], p09_bt100_adaptive[,1],
                     p09_bt150_adaptive[,1], p09_bt200_adaptive[,1], p09_bt250_adaptive[,1],
                     p09_bt300_adaptive[,1], p09_bt350_adaptive[,1], p09_bt400_adaptive[,1],
                     p09_bt450_adaptive[,1], p09_bt500_adaptive[,1])

p1_adaptive = cbind(p1_bt1_adaptive[,1], p1_bt50_adaptive[,1], p1_bt100_adaptive[,1],
                     p1_bt150_adaptive[,1], p1_bt200_adaptive[,1], p1_bt250_adaptive[,1],
                     p1_bt300_adaptive[,1], p1_bt350_adaptive[,1], p1_bt400_adaptive[,1],
                     p1_bt450_adaptive[,1], p1_bt500_adaptive[,1])

p1_adaptive[100,] = 0

heat_adap = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_adap[i,1] = sqrt(sum((p0_beta[,i] - p0_adaptive[,i])^2))
  heat_adap[i,2] = sqrt(sum((p01_beta[,i] - p01_adaptive[,i])^2))
  heat_adap[i,3] = sqrt(sum((p02_beta[,i] - p02_adaptive[,i])^2))
  heat_adap[i,4] = sqrt(sum((p03_beta[,i] - p03_adaptive[,i])^2))
  heat_adap[i,5] = sqrt(sum((p04_beta[,i] - p04_adaptive[,i])^2))
  heat_adap[i,6] = sqrt(sum((p05_beta[,i] - p05_adaptive[,i])^2))
  heat_adap[i,7] = sqrt(sum((p06_beta[,i] - p06_adaptive[,i])^2))
  heat_adap[i,8] = sqrt(sum((p07_beta[,i] - p07_adaptive[,i])^2))
  heat_adap[i,9] = sqrt(sum((p08_beta[,i] - p08_adaptive[,i])^2))
  heat_adap[i,10] = sqrt(sum((p09_beta[,i] - p09_adaptive[,i])^2))
  heat_adap[i,11] = sqrt(sum((p1_beta[,i] - p1_adaptive[,i])^2))
}

colnames(heat_adap) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_adap) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_adap

###heatmap normal bruta
x= seq(0,1, length.out=100)
p0_norm = cbind(dnorm(x, mean = 0, sd = sqrt((1 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((50 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((100 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((150 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((200 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((250 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((300 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((350 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((400 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((450 / (200)) * (0*1) )),
                dnorm(x, mean = 0, sd = sqrt((500 / (200)) * (0*1) )))
p0_norm[1,] = 0

p01_norm = cbind(dnorm(x, mean = 0.1, sd = sqrt((1 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((50 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((100 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((150 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((200 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((250 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((300 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((350 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((400 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((450 / (200)) * (0.1*0.9) )),
                dnorm(x, mean = 0.1, sd = sqrt((500 / (200)) * (0.1*0.9) )))

p02_norm = cbind(dnorm(x, mean = 0.2, sd = sqrt((1 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((50 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((100 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((150 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((200 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((250 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((300 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((350 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((400 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((450 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.2, sd = sqrt((500 / (200)) * (0.2*0.8) )))

p03_norm = cbind(dnorm(x, mean = 0.3, sd = sqrt((1 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((50 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((100 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((150 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((200 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((250 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((300 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((350 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((400 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((450 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.3, sd = sqrt((500 / (200)) * (0.3*0.7) )))

p04_norm = cbind(dnorm(x, mean = 0.4, sd = sqrt((1 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((50 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((100 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((150 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((200 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((250 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((300 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((350 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((400 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((450 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.4, sd = sqrt((500 / (200)) * (0.4*0.6) )))

p05_norm = cbind(dnorm(x, mean = 0.5, sd = sqrt((1 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((50 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((100 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((150 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((200 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((250 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((300 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((350 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((400 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((450 / (200)) * (0.5*0.5) )),
                 dnorm(x, mean = 0.5, sd = sqrt((500 / (200)) * (0.5*0.5) )))

p06_norm = cbind(dnorm(x, mean = 0.6, sd = sqrt((1 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((50 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((100 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((150 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((200 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((250 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((300 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((350 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((400 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((450 / (200)) * (0.4*0.6) )),
                 dnorm(x, mean = 0.6, sd = sqrt((500 / (200)) * (0.4*0.6) )))

p07_norm = cbind(dnorm(x, mean = 0.7, sd = sqrt((1 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((50 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((100 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((150 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((200 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((250 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((300 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((350 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((400 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((450 / (200)) * (0.3*0.7) )),
                 dnorm(x, mean = 0.7, sd = sqrt((500 / (200)) * (0.3*0.7) )))

p08_norm = cbind(dnorm(x, mean = 0.8, sd = sqrt((1 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((50 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((100 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((150 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((200 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((250 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((300 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((350 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((400 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((450 / (200)) * (0.2*0.8) )),
                 dnorm(x, mean = 0.8, sd = sqrt((500 / (200)) * (0.2*0.8) )))

p09_norm = cbind(dnorm(x, mean = 0.9, sd = sqrt((1 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((50 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((100 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((150 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((200 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((250 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((300 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((350 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((400 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((450 / (200)) * (0.1*0.9) )),
                 dnorm(x, mean = 0.9, sd = sqrt((500 / (200)) * (0.1*0.9) )))

p1_norm = cbind(rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100),
                rep(0, 100))

heat_norm = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_norm[i,1] = sqrt(sum((p0_adaptive[,i] - p0_norm[,i])^2))
  heat_norm[i,2] = sqrt(sum((p01_adaptive[,i] - p01_norm[,i])^2))
  heat_norm[i,3] = sqrt(sum((p02_adaptive[,i] - p02_norm[,i])^2))
  heat_norm[i,4] = sqrt(sum((p03_adaptive[,i] - p03_norm[,i])^2))
  heat_norm[i,5] = sqrt(sum((p04_adaptive[,i] - p04_norm[,i])^2))
  heat_norm[i,6] = sqrt(sum((p05_adaptive[,i] - p05_norm[,i])^2))
  heat_norm[i,7] = sqrt(sum((p06_adaptive[,i] - p06_norm[,i])^2))
  heat_norm[i,8] = sqrt(sum((p07_adaptive[,i] - p07_norm[,i])^2))
  heat_norm[i,9] = sqrt(sum((p08_adaptive[,i] - p08_norm[,i])^2))
  heat_norm[i,10] = sqrt(sum((p09_adaptive[,i] - p09_norm[,i])^2))
  heat_norm[i,11] = sqrt(sum((p1_adaptive[,i] - p1_norm[,i])^2))
}

colnames(heat_norm) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_norm) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_norm

######## Distancia de Hellinger

##1 asymtotical expansion
heat_hell_expan = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_hell_expan[i,1] = (1/2)*sqrt(sum((sqrt(p0_adaptive[,i]) - sqrt(rep(0, 100)))^2))
  heat_hell_expan[i,2] = (1/2)*sqrt(sum((sqrt(p01_adaptive[,i]) - sqrt(expan(x0[2], 100, t[i])))^2))
  heat_hell_expan[i,3] = (1/2)*sqrt(sum((sqrt(p02_adaptive[,i]) - sqrt(expan(x0[3], 100, t[i])))^2))
  heat_hell_expan[i,4] = (1/2)*sqrt(sum((sqrt(p03_adaptive[,i]) - sqrt(expan(x0[4], 100, t[i])))^2))
  heat_hell_expan[i,5] = (1/2)*sqrt(sum((sqrt(p04_adaptive[,i]) - sqrt(expan(x0[5], 100, t[i])))^2))
  heat_hell_expan[i,6] = (1/2)*sqrt(sum((sqrt(p05_adaptive[,i]) - sqrt(expan(x0[6], 100, t[i])))^2))
  heat_hell_expan[i,7] = (1/2)*sqrt(sum((sqrt(p06_adaptive[,i]) - sqrt(expan(x0[7], 100, t[i])))^2))
  heat_hell_expan[i,8] = (1/2)*sqrt(sum((sqrt(p07_adaptive[,i]) - sqrt(expan(x0[8], 100, t[i])))^2))
  heat_hell_expan[i,9] = (1/2)*sqrt(sum((sqrt(p08_adaptive[,i]) - sqrt(expan(x0[9], 100, t[i])))^2))
  heat_hell_expan[i,10] = (1/2)*sqrt(sum((sqrt(p09_adaptive[,i]) - sqrt(expan(x0[10], 100, t[i])))^2))
  heat_hell_expan[i,11] = (1/2)*sqrt(sum((sqrt(p1_adaptive[,i]) - sqrt(rep(0, 100)))^2))
}


colnames(heat_hell_expan) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_hell_expan) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_hell_expan

##2 beta con parametros adecuados
heat_hell_beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_hell_beta[i,1] = (1/2)*sqrt(sum((sqrt(p0_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,1], shape2 = sh2[i,1]), 0)))^2))
  heat_hell_beta[i,2] = (1/2)*sqrt(sum((sqrt(p01_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,2], shape2 = sh2[i,2]), 0)))^2))
  heat_hell_beta[i,3] = (1/2)*sqrt(sum((sqrt(p02_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,3], shape2 = sh2[i,3]), 0)))^2))
  heat_hell_beta[i,4] = (1/2)*sqrt(sum((sqrt(p03_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,4], shape2 = sh2[i,4]), 0)))^2))
  heat_hell_beta[i,5] = (1/2)*sqrt(sum((sqrt(p04_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,5], shape2 = sh2[i,5]), 0)))^2))
  heat_hell_beta[i,6] = (1/2)*sqrt(sum((sqrt(p05_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,6], shape2 = sh2[i,6]), 0)))^2))
  heat_hell_beta[i,7] = (1/2)*sqrt(sum((sqrt(p06_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,7], shape2 = sh2[i,7]), 0)))^2))
  heat_hell_beta[i,8] = (1/2)*sqrt(sum((sqrt(p07_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,8], shape2 = sh2[i,8]), 0)))^2))
  heat_hell_beta[i,9] = (1/2)*sqrt(sum((sqrt(p08_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,9], shape2 = sh2[i,9]), 0)))^2))
  heat_hell_beta[i,10] = (1/2)*sqrt(sum((sqrt(p09_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,10], shape2 = sh2[i,10]), 0)))^2))
  heat_hell_beta[i,11] = (1/2)*sqrt(sum((sqrt(p1_adaptive[,i]) - sqrt(c(0, dbeta(x[2:99], shape1 = sh1[i,11], shape2 = sh2[i,11]), 0)))^2))
}

colnames(heat_hell_beta) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_hell_beta) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_hell_beta


##3 normal con parametros adecuados
heat_hell_nor = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  heat_hell_nor[i,1] = (1/2)*sqrt(sum((sqrt(p0_adaptive[,i]) - sqrt(p0_norm[,i]))^2))
  heat_hell_nor[i,2] = (1/2)*sqrt(sum((sqrt(p01_adaptive[,i]) - sqrt(p01_norm[,i]))^2))
  heat_hell_nor[i,3] = (1/2)*sqrt(sum((sqrt(p02_adaptive[,i]) - sqrt(p02_norm[,i]))^2))
  heat_hell_nor[i,4] = (1/2)*sqrt(sum((sqrt(p03_adaptive[,i]) - sqrt(p03_norm[,i]))^2))
  heat_hell_nor[i,5] = (1/2)*sqrt(sum((sqrt(p04_adaptive[,i]) - sqrt(p04_norm[,i]))^2))
  heat_hell_nor[i,6] = (1/2)*sqrt(sum((sqrt(p05_adaptive[,i]) - sqrt(p05_norm[,i]))^2))
  heat_hell_nor[i,7] = (1/2)*sqrt(sum((sqrt(p06_adaptive[,i]) - sqrt(p06_norm[,i]))^2))
  heat_hell_nor[i,8] = (1/2)*sqrt(sum((sqrt(p07_adaptive[,i]) - sqrt(p07_norm[,i]))^2))
  heat_hell_nor[i,9] = (1/2)*sqrt(sum((sqrt(p08_adaptive[,i]) - sqrt(p08_norm[,i]))^2))
  heat_hell_nor[i,10] = (1/2)*sqrt(sum((sqrt(p09_adaptive[,i]) - sqrt(p09_norm[,i]))^2))
  heat_hell_nor[i,11] = (1/2)*sqrt(sum((sqrt(p1_adaptive[,i]) - sqrt(p1_norm[,i]))^2))
}

colnames(heat_hell_nor) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_hell_nor) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
heat_hell_nor

##4 Gaussian approximation
heat_hell_gauss = matrix(NA, ncol = 11, nrow = 11)

for (i in 1:11) {
  heat_hell_gauss[i,1] = (1/2)*sqrt(sum((sqrt(p0_adaptive[,i]) - sqrt(rep(0, 100)))^2))
  heat_hell_gauss[i,2] = (1/2)*sqrt(sum((sqrt(p01_adaptive[,i]) - sqrt(gauss(x0[2], 100, t[i])))^2))
  heat_hell_gauss[i,3] = (1/2)*sqrt(sum((sqrt(p02_adaptive[,i]) - sqrt(gauss(x0[3], 100, t[i])))^2))
  heat_hell_gauss[i,4] = (1/2)*sqrt(sum((sqrt(p03_adaptive[,i]) - sqrt(gauss(x0[4], 100, t[i])))^2))
  heat_hell_gauss[i,5] = (1/2)*sqrt(sum((sqrt(p04_adaptive[,i]) - sqrt(gauss(x0[5], 100, t[i])))^2))
  heat_hell_gauss[i,6] = (1/2)*sqrt(sum((sqrt(p05_adaptive[,i]) - sqrt(gauss(x0[6], 100, t[i])))^2))
  heat_hell_gauss[i,7] = (1/2)*sqrt(sum((sqrt(p06_adaptive[,i]) - sqrt(gauss(x0[7], 100, t[i])))^2))
  heat_hell_gauss[i,8] = (1/2)*sqrt(sum((sqrt(p07_adaptive[,i]) - sqrt(gauss(x0[8], 100, t[i])))^2))
  heat_hell_gauss[i,9] = (1/2)*sqrt(sum((sqrt(p08_adaptive[,i]) - sqrt(gauss(x0[9], 100, t[i])))^2))
  heat_hell_gauss[i,10] = (1/2)*sqrt(sum((sqrt(p09_adaptive[,i]) - sqrt(gauss(x0[10], 100, t[i])))^2))
  heat_hell_gauss[i,11] = (1/2)*sqrt(sum((sqrt(p1_adaptive[,i]) - sqrt(rep(0, 100)))^2))
}

heat_hell_gauss
colnames(heat_hell_gauss) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(heat_hell_gauss) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")


###heatmaps
library("lattice")
levelplot(t(scale(heat_expan)), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Asymptotical expansion')
levelplot(t(heat_expan), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Asymptotical expansion')

levelplot(t(scale(heat_gauss)), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Gaussian approximation')
levelplot(t(heat_gauss), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Gaussian approximation')

levelplot(t(scale(heat_beta)), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Beta Distribution')
levelplot(t(heat_beta), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Beta Distribution')

levelplot(t(scale(heat_beta_sp)), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Beta with \n spikes distribution')

levelplot(t(heat_beta_sp), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Beta with \n spikes distribution')

library(gridExtra)
g1 = levelplot(t(scale(heat_gauss)), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Gaussian \n approximation', cex.lab=0.75)

g2 = levelplot(t(scale(heat_expan)), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Asymptotical \n expansion', cex.lab=0.75)

g3 = levelplot(t(scale(heat_beta)), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Beta \n Distribution', cex.lab=0.75)

g4 = levelplot(t(scale(heat_norm)), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Beta with \n spikes distribution', cex.lab=0.75)

g5 = levelplot(t(heat_gauss), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Gaussian \n approximation', cex.lab=0.75)

g6 = levelplot(t(heat_expan), col.regions=heat.colors(200), xlab = expression(p[0]),
              ylab = "t", main = 'Asymptotical \n expansion', cex.lab=0.75)

g7 = levelplot(t(heat_beta), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Beta \n Distribution', cex.lab=0.75)

g8 = levelplot(t(heat_norm), col.regions=heat.colors(200), xlab = expression(p[0]),
               ylab = "t", main = 'Gaussian \n distribution', cex.lab=0.75)

levelplot(t(heat_adap), col.regions=heat.colors(200), xlab = expression(p[0]),
          ylab = "t", main = 'Beta with \n spikes distribution', cex.lab=0.75)

plot(g5)
plot(g6)
plot(g7)
plot(g8)


grid.arrange(g6, g5, ncol=2)
grid.arrange(g7, g8, ncol=2)

##lines plot para tres valores de t fijos y valor inicial p0=0.5
#t=100

m = matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3, byrow = T)
layout(mat = m, heights = c(0.8,0.2))

x = seq(from=0.001, to=1, by=0.001)
par(mar=c(2,2,1,1))
hist(p0_inicial01[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.1, 1000, 100), col = "green", lwd = 2)
lines(x, gauss(0.1, 1000, 100), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.1, sd = sqrt((100 / (1000)) * (0.1*0.9) )) , col = "red", lwd = 2)
lines(p01_bt100_adaptive[,2], p01_bt100_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial05[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F)
lines(x, expan(0.5, 1000, 100), col = "green", lwd = 2)
lines(x, gauss(0.5, 1000, 100), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.5, sd = sqrt((100 / (1000)) * (0.5*0.5) )) , col = "red", lwd = 2)
lines(p05_bt100_adaptive[,2], p05_bt100_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial09[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.9, 1000, 100), col = "green", lwd = 2)
lines(x, gauss(0.9, 1000, 100), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.9, sd = sqrt((100 / (1000)) * (0.9*0.1) )) , col = "red", lwd = 2)
lines(p09_bt100_adaptive[,2], p09_bt100_adaptive[,1], lwd = 2, col = 'darkblue')

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors = c('darkblue', 'red', 'green', 'purple')
legend(x = "top",inset = 0,
       legend = c("ADE", "Gauss", "AE", "GaussA"), 
       col=plot_colors, lwd=5, cex=1.2, horiz = TRUE)


#t=200
hist(p0_inicial01[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.1, 1000, 200), col = "green", lwd = 2)
lines(x, gauss(0.1, 1000, 200), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.1, sd = sqrt((200 / (1000)) * (0.1*0.9) )) , col = "red", lwd = 2)
lines(p01_bt200_adaptive[,2], p01_bt200_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial05[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F)
lines(x, expan(0.5, 1000, 200), col = "green", lwd = 2)
lines(x, gauss(0.5, 1000, 200), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.5, sd = sqrt((200 / (1000)) * (0.5*0.5) )) , col = "red", lwd = 2)
lines(p05_bt200_adaptive[,2], p05_bt200_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial09[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.9, 1000, 200), col = "green", lwd = 2)
lines(x, gauss(0.9, 1000, 200), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.9, sd = sqrt((200 / (1000)) * (0.9*0.1) )) , col = "red", lwd = 2)
lines(p09_bt200_adaptive[,2], p09_bt200_adaptive[,1], lwd = 2, col = 'darkblue')

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors = c('darkblue', 'red', 'green', 'purple')
legend(x = "top",inset = 0,
       legend = c("ADE", "Gauss", "AE", "GaussA"), 
       col=plot_colors, lwd=5, cex=1.2, horiz = TRUE)

#t=300
hist(p0_inicial01[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.1, 1000, 300), col = "green", lwd = 2)
lines(x, gauss(0.1, 1000, 300), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.1, sd = sqrt((300 / (1000)) * (0.1*0.9) )) , col = "red", lwd = 2)
lines(p01_bt300_adaptive[,2], p01_bt300_adaptive[,1], lwd = 2, col = 'darkblue')

hist(p0_inicial05[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F)
lines(x, expan(0.5, 1000, 300), col = "green", lwd = 2)
lines(x, gauss(0.5, 1000, 300), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.5, sd = sqrt((300 / (1000)) * (0.5*0.5) )) , col = "red", lwd = 2)
lines(p05_bt300_adaptive[,2], p05_bt300_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial09[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.9, 1000, 300), col = "green", lwd = 2)
lines(x, gauss(0.9, 1000, 300), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.9, sd = sqrt((300 / (1000)) * (0.9*0.1) )) , col = "red", lwd = 2)
lines(p09_bt300_adaptive[,2], p09_bt300_adaptive[,1], lwd = 2, col = 'darkblue')

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors = c('darkblue', 'red', 'green', 'purple')
legend(x = "top",inset = 0,
       legend = c("ADE", "Gauss", "AE", "GaussA"), 
       col=plot_colors, lwd=5, cex=1.2, horiz = TRUE)

#t=400
hist(p0_inicial01[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.1, 1000, 400), col = "green", lwd = 2)
lines(x, gauss(0.1, 1000, 400), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.1, sd = sqrt((400 / (1000)) * (0.1*0.9) )) , col = "red", lwd = 2)
lines(p01_bt400_adaptive[,2], p01_bt400_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial05[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F)
lines(x, expan(0.5, 1000, 400), col = "green", lwd = 2)
lines(x, gauss(0.5, 1000, 400), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.5, sd = sqrt((400 / (1000)) * (0.5*0.5) )) , col = "red", lwd = 2)
lines(p05_bt400_adaptive[,2], p05_bt400_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial09[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.9, 1000, 400), col = "green", lwd = 2)
lines(x, gauss(0.9, 1000, 400), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.9, sd = sqrt((400 / (1000)) * (0.9*0.1) )) , col = "red", lwd = 2)
lines(p09_bt400_adaptive[,2], p09_bt400_adaptive[,1], lwd = 2, col = 'darkblue')

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors = c('darkblue', 'red', 'green', 'purple')
legend(x = "top",inset = 0,
       legend = c("ADE", "Gauss", "AE", "GaussA"), 
       col=plot_colors, lwd=5, cex=1.2, horiz = TRUE)

#t=500
hist(p0_inicial01[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.1, 1000, 500), col = "green", lwd = 2)
lines(x, gauss(0.1, 1000, 500), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.1, sd = sqrt((500 / (1000)) * (0.1*0.9) )) , col = "red", lwd = 2)
lines(p01_bt500_adaptive[,2], p01_bt500_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial05[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F)
lines(x, expan(0.5, 1000, 500), col = "green", lwd = 2)
lines(x, gauss(0.5, 1000, 500), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.5, sd = sqrt((500 / (1000)) * (0.5*0.5) )) , col = "red", lwd = 2)
lines(p05_bt500_adaptive[,2], p05_bt500_adaptive[,1], lwd = 2, col = 'darkblue')


hist(p0_inicial09[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,6),freq = F)
lines(x, expan(0.9, 1000, 500), col = "green", lwd = 2)
lines(x, gauss(0.9, 1000, 500), col = "purple", lwd = 2)
lines(x, dnorm(x, mean = 0.9, sd = sqrt((500 / (1000)) * (0.9*0.1) )) , col = "red", lwd = 2)
lines(p09_bt500_adaptive[,2], p09_bt500_adaptive[,1], lwd = 2, col = 'darkblue')

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors = c('darkblue', 'red', 'green', 'purple')
legend(x = "top",inset = 0,
       legend = c("ADE", "Gauss", "AE", "GaussA"), 
       col=plot_colors, lwd=5, cex=1.2, horiz = TRUE)

#### visualizacion de la mejor del asymptotical expansion
round(heat_beta, 2)
round(heat_norm, 2)
round(heat_gauss, 2)
round(heat_expan, 2)




plot(t, heat_norm[,6], type = 'o', lwd=2, ylim = c(0,60) ,ylab = expression(L^2-norm), main = bquote(paste(p[0], "=0.5")), col = 'red1')
lines(t, heat_beta[,6], lwd=2, type = 'o', col = 'yellow2')
lines(t, heat_expan[,6], lwd=2, type = 'o', col = 'red4')
lines(t, heat_gauss[,6], lwd=2, type = 'o', col = 'darkorange2')
lines(t, heat_norm[,6], type = 'o', lwd=2, col = 'red1')
legend("topleft",
       legend =  c("Beta", "Gauss", "GaussA", "AE"),
       pch = c(1, 1, 1, 1),
       col = c('yellow2', 'red1',  'darkorange2','red4'),
       lwd = 2,
       horiz = T,
       cex = 0.8)


plot(t, heat_norm[,10], type = 'o', lwd=2, ylim = c(0,60) ,ylab = expression(L^2-norm), main = bquote(paste(p[0], "=0.9")), col = 'red1')
lines(t, heat_beta[,10], lwd=2, type = 'o', col = 'yellow2')
lines(t, heat_expan[,10], lwd=2, type = 'o', col = 'red4')
lines(t, heat_gauss[,10], lwd=2, type = 'o', col = 'darkorange2')
legend("topleft",
       legend =  c("Beta", "Gauss", "GaussA", "AE"),
       pch = c(1, 1, 1, 1),
       col = c('yellow2', 'red1',  'darkorange2','red4'),
       lwd = 2,
       horiz = T,
       cex = 0.8)


plot(t, heat_norm[,2], type = 'o', lwd=2, ylim = c(0,60) ,ylab = expression(L^2-norm), main = bquote(paste(p[0], "=0.1")), col = 'red1')
lines(t, heat_beta[,2], lwd=2, type = 'o', col = 'yellow2')
lines(t, heat_expan[,2], lwd=2, type = 'o', col = 'red4')
lines(t, heat_gauss[,2], lwd=2, type = 'o', col = 'darkorange2')
legend("topleft",
       legend =  c("Beta", "Gauss", "GaussA", "AE"),
       pch = c(1, 1, 1, 1),
       col = c('yellow2', 'red1',  'darkorange2','red4'),
       lwd = 2,
       horiz = T,
       cex = 0.8)

##otros gráficos
library('plot.matrix')
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(t(heat_gauss), border=NA, xlab ='t', ylab = expression(p[0]), main = '')


### heatmaps con ggplot 

##hacer un data frame con 484 filas y 4 columnas

t_2N = rep(t/1000, 11)
p0_df = c(rep(0,11), rep(0.1, 11), rep(0.2, 11), rep(0.3, 11), rep(0.4, 11), rep(0.5, 11), rep(0.6, 11), rep(0.7, 11), rep(0.8, 11), rep(0.9, 11), rep(1.0, 11))

expan_df = cbind.data.frame(t_2N, p0_df, rep("AE", 121), as.vector(heat_expan))
colnames(expan_df) = c("1","2","3","4")

gauss_df = cbind.data.frame(t_2N, p0_df, rep("GaussA", 121), as.vector(heat_gauss))
colnames(gauss_df) = c("1","2","3","4")

beta_df = cbind.data.frame(t_2N, p0_df, rep("Beta", 121), as.vector(heat_beta))
colnames(beta_df) = c("1","2","3","4")

norm_df = cbind.data.frame(t_2N, p0_df, rep("Gauss", 121), as.vector(heat_norm))
colnames(norm_df) = c("1","2","3","4")

heat_df = rbind(expan_df, gauss_df, norm_df,  beta_df)
colnames(heat_df) = c("t/2N", "p0", "dist", "L2")

dim(heat_df)


##librerías
library(ggplot2)       #
library(scales)        # working with ggplot2 for label formatting
library(gridExtra)     # working with ggplots for arranging plots
library(ggthemes)      # clean theme for ggplot2
library(viridis)       # color palette
library(DT)

gg = ggplot(heat_df[1:121,], aes(x=heat_df[1:121,]$p0, y=heat_df[1:121,]$`t/2N`, fill=heat_df[1:121,]$L2))
gg = gg + geom_tile(color="white", size=0.1)
gg = gg + scale_fill_viridis(name=expression(L^2), label=comma)
gg = gg + coord_equal()
gg = gg + labs(x=NULL, y=NULL, title="Titulo")
gg = gg + theme_tufte(base_family="Helvetica")
gg = gg + theme(plot.title=element_text(hjust=0))
#gg = gg + theme(axis.ticks = element_blank())
#gg = gg + theme(axis.text = element_text(size=7))

gg = gg + theme(legend.title = element_text(size=8))
gg = gg + theme(legend.text = element_text(size=6))

gg


###juntando 
gg <- ggplot(heat_df, aes(x=p0, y=`t/2N`, fill=L2))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y=expression(t/2*N), x = expression(p[0]))
gg <- gg + scale_fill_viridis(name=expression(L^2))
gg <- gg + coord_equal()
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



######## Distancia de Hellinger
## distancia de Hellinger con paquete
#install.packages('statip ')
library(statip)
library(plot.matrix)

##1 asymtotical expansion
x0
t = t/200; t

hh_expan = matrix(ncol = 11, nrow = 11)
for (i in 1:11) {
  hh_expan[i,1] = hellinger(p0_adaptive[,i], expan(x0[1], 100, t[i]), 0, 1) 
  hh_expan[i,2] = hellinger(p01_adaptive[,i], expan(x0[2], 100, t[i]), 0, 1)
  hh_expan[i,3] = hellinger(p02_adaptive[,i], expan(x0[3], 100, t[i]), 0, 1)
  hh_expan[i,4] = hellinger(p03_adaptive[,i], expan(x0[4], 100, t[i]), 0, 1)
  hh_expan[i,5] = hellinger(p04_adaptive[,i], expan(x0[5], 100, t[i]), 0, 1)
  hh_expan[i,6] = hellinger(p05_adaptive[,i], expan(x0[6], 100, t[i]), 0, 1)
  hh_expan[i,7] = hellinger(p06_adaptive[,i], expan(x0[7], 100, t[i]), 0, 1)
  hh_expan[i,8] = hellinger(p07_adaptive[,i], expan(x0[8], 100, t[i]), 0, 1)
  hh_expan[i,9] = hellinger(p08_adaptive[,i], expan(x0[9], 100, t[i]), 0, 1)
  hh_expan[i,10] = hellinger(p09_adaptive[,i], expan(x0[10], 100, t[i]), 0, 1)
  hh_expan[i,11] = hellinger(p1_adaptive[,i], expan(x0[11], 100, t[i]), 0, 1)
}

colnames(hh_expan) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_expan) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_expan


##2 beta con parametros adecuados
hh_beta = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh_beta[i,1] = hellinger(p0_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,1], shape2 = sh2[i,1]), 0), 0, 1)
  hh_beta[i,2] = hellinger(p01_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,2], shape2 = sh2[i,2]), 0), 0, 1)
  hh_beta[i,3] = hellinger(p02_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,3], shape2 = sh2[i,3]), 0), 0, 1)
  hh_beta[i,4] = hellinger(p03_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,4], shape2 = sh2[i,4]), 0), 0, 1)
  hh_beta[i,5] = hellinger(p04_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,5], shape2 = sh2[i,5]), 0), 0, 1)
  hh_beta[i,6] = hellinger(p05_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,6], shape2 = sh2[i,6]), 0), 0, 1)
  hh_beta[i,7] = hellinger(p06_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,7], shape2 = sh2[i,7]), 0), 0, 1)
  hh_beta[i,8] = hellinger(p07_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,8], shape2 = sh2[i,8]), 0), 0, 1)
  hh_beta[i,9] = hellinger(p08_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,9], shape2 = sh2[i,9]), 0), 0, 1)
  hh_beta[i,10] = hellinger(p09_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,10], shape2 = sh2[i,10]), 0), 0, 1)
  hh_beta[i,11] = hellinger(p1_adaptive[,i], c(0, dbeta(x[2:99], shape1 = sh1[i,11], shape2 = sh2[i,11]), 0), 0, 1)
}

colnames(hh_beta) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_beta) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_beta


##3 normal con parametros adecuados
hh_nor = matrix(NA, ncol = 11, nrow = 11)
for (i in 1:11) {
  hh_nor[i,1] = hellinger(p0_adaptive[,i],  p0_norm[,i], 0, 1)
  hh_nor[i,2] = hellinger(p01_adaptive[,i], p01_norm[,i], 0, 1)
  hh_nor[i,3] = hellinger(p02_adaptive[,i], p02_norm[,i], 0, 1)
  hh_nor[i,4] = hellinger(p03_adaptive[,i], p03_norm[,i], 0, 1)
  hh_nor[i,5] = hellinger(p04_adaptive[,i], p04_norm[,i], 0, 1)
  hh_nor[i,6] = hellinger(p05_adaptive[,i], p05_norm[,i], 0, 1)
  hh_nor[i,7] = hellinger(p06_adaptive[,i], p06_norm[,i], 0, 1)
  hh_nor[i,8] = hellinger(p07_adaptive[,i], p07_norm[,i], 0, 1)
  hh_nor[i,9] = hellinger(p08_adaptive[,i], p08_norm[,i], 0, 1)
  hh_nor[i,10] = hellinger(p09_adaptive[,i],p09_norm[,i], 0, 1)
  hh_nor[i,11] = hellinger(p1_adaptive[,i], p1_norm[,i], 0, 1)
}

colnames(hh_nor) = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
rownames(hh_nor) = c("1", "50", "100", "150", "200", "250","300", "350", "400", "450", "500")
hh_nor

##4 Gaussian approximation
hh_gauss = matrix(NA, ncol = 11, nrow = 11)

for (i in 1:11) {
  hh_gauss[i,1] = hellinger(p0_adaptive[,i], rep(0, 100), 0, 1)
  hh_gauss[i,2] = hellinger(p01_adaptive[,i], gauss(x0[2], 100, t[i]), 0, 1)
  hh_gauss[i,3] = hellinger(p02_adaptive[,i], gauss(x0[3], 100, t[i]), 0, 1)
  hh_gauss[i,4] = hellinger(p03_adaptive[,i], gauss(x0[4], 100, t[i]), 0, 1)
  hh_gauss[i,5] = hellinger(p04_adaptive[,i], gauss(x0[5], 100, t[i]), 0, 1)
  hh_gauss[i,6] = hellinger(p05_adaptive[,i], gauss(x0[6], 100, t[i]), 0, 1)
  hh_gauss[i,7] = hellinger(p06_adaptive[,i], gauss(x0[7], 100, t[i]), 0, 1)
  hh_gauss[i,8] = hellinger(p07_adaptive[,i], gauss(x0[8], 100, t[i]), 0, 1)
  hh_gauss[i,9] = hellinger(p08_adaptive[,i], gauss(x0[9], 100, t[i]), 0, 1)
  hh_gauss[i,10] = hellinger(p09_adaptive[,i], gauss(x0[10], 100, t[i]), 0, 1)
  hh_gauss[i,11] = hellinger(p1_adaptive[,i], rep(0, 100), 0, 1)
}

hh_gauss

##### con la distancia de Hellinger
t_2N = rep(t/1000, 11)
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
hell_heat_df = rbind(hell_expan_df, hell_gauss_df, hell_norm_df, hell_beta_df)
colnames(hell_heat_df) = c("t/2N", "p0", "dist", "L2")

dim(hell_heat_df)

gg = ggplot(hell_heat_df[1:121,], aes(x=hell_heat_df[1:121,]$p0, y=hell_heat_df[1:121,]$`t/2N`, fill=hell_heat_df[1:121,]$L2))
gg = gg + geom_tile(color="white", size=0.1)
gg = gg + scale_fill_viridis(name=expression(L^2), label=comma)
gg = gg + coord_equal()
gg = gg + labs(x=NULL, y=NULL, title="Titulo")
gg = gg + theme_tufte(base_family="Helvetica")
gg = gg + theme(plot.title=element_text(hjust=0))
#gg = gg + theme(axis.ticks = element_blank())
#gg = gg + theme(axis.text = element_text(size=7))

gg = gg + theme(legend.title = element_text(size=8))
gg = gg + theme(legend.text = element_text(size=6))

gg


###juntando 
gg <- ggplot(hell_heat_df, aes(x=p0, y=`t/2N`, fill=L2))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + labs(y=expression(t/2*N), x = expression(p[0]))
gg <- gg + scale_fill_viridis(name="Hellinger")
gg <- gg + coord_equal()
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

