###### W-F Difussion 

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

###para graficar las trayectorias
rep = 1000
repeticiones05 = replicate(rep, trayectoria(1000, 0.5, 500))
repeticiones01 = replicate(rep, trayectoria(1000, 0.1, 500))
repeticiones09 = replicate(rep, trayectoria(1000, 0.9, 500))

### agregar media y desviacion estandar para cada t
#p0=0.5
plot(repeticiones05[,1], type = "l", ylim = c(0,1), xlab = "t", ylab = "p(t)", col=rgb(0,0,0,0.1))
for (i in 2:rep) {
  lines(repeticiones05[,i], col=rgb(0,0,0,0.3))
}


##using ggplot
library(ggplot2)
t=rep(1:500, 1000)
dataf_trayec = data.frame(as.vector(repeticiones05), t)

ggplot(dataf_trayec, aes(x = t, y = repeticiones05)) +
  labs(x = "t", y = expression(x[0])) + 
  #scale_color_manual(values=rep(rgb(0,0,0,0.1), 1000)) + 
  theme(legend.position='none') +
  theme_bw() +
  geom_line(col=rgb(0,0,0,0.15))
#######


t_med = apply(repeticiones05, 1, mean)
t_sd = apply(repeticiones05, 1, sd)

lines(t_med, col = "green", lwd = 2)
lines(t_med + t_sd, col = "red", lwd = 2)
lines(t_med - t_sd, col = "red", lwd = 2)
lines(t_med + 2*t_sd, col = "orange", lwd = 2)
lines(t_med - 2*t_sd, col = "orange", lwd = 2)
lines(t_med + 3*t_sd, col = "blue", lwd = 2)
lines(t_med - 3*t_sd, col = "blue", lwd = 2)
legend("bottomleft", legend=c("Mean", "Standard deviation","2*Standard deviation","3*Standard deviation"),
       col=c("green", "red", "orange","blue"), lty=1, cex=0.5)

#p0=0.1
plot(repeticiones01[,1], type = "l", ylim = c(0,1), xlab = "t", ylab = "p(t)", col=rgb(0,0,0,0.1))
for (i in 2:rep) {
  lines(repeticiones01[,i], col=rgb(0,0,0,0.3))
}

t_med01 = apply(repeticiones01, 1, mean)
t_sd01 = apply(repeticiones01, 1, sd)

lines(t_med01, col = "green", lwd = 2)
lines(t_med01 + t_sd01, col = "red", lwd = 2)
lines(t_med01 + 2*t_sd01, col = "orange", lwd = 2)
lines(t_med01 + 3*t_sd01, col = "blue", lwd = 2)
legend("topleft", legend=c("Mean", "Standard deviation","2*Standard deviation","3*Standard deviation"),
       col=c("green", "red", "orange","blue"), lty=1, cex=0.5)


#p0=0.9
plot(repeticiones09[,1], type = "l", ylim = c(0,1), xlab = "t", ylab = "p(t)", col=rgb(0,0,0,0.1))
for (i in 2:rep) {
  lines(repeticiones09[,i], col=rgb(0,0,0,0.3))
}

t_med09 = apply(repeticiones09, 1, mean)
t_sd09 = apply(repeticiones09, 1, sd)

lines(t_med09, col = "green", lwd = 2)
lines(t_med09 - t_sd09, col = "red", lwd = 2)
lines(t_med09 - 2*t_sd09, col = "orange", lwd = 2)
lines(t_med09 - 3*t_sd09, col = "blue", lwd = 2)
legend("bottomleft", legend=c("Mean", "Standard deviation","2*Standard deviation","3*Standard deviation"),
       col=c("green", "red", "orange","blue"), lty=1, cex=0.5)


## extraer una generacion t en particular
#t = c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)

### para mejorar la simulacion se tomaron 1e6 repeticiones
### pero solo para 
### t = c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)

repe = 1e6 
## p0 = 0.5
start.time = proc.time()
rep_trayec05 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec05[i,] = trayectoria(1000, 0.5, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec05, file = "rep_05.txt", sep = " ") #para guardar los resultados
proc.time() - start.time 


## p0 = 0.1
start.time = proc.time()
p0_inicial01 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec01[i,] = trayectoria(1000, 0.1, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec01, file = "rep_01.txt", sep = " ") #para guardar los resultados
proc.time() - start.time 

## p0 = 0.9
start.time = proc.time()
rep_trayec09 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec09[i,] = trayectoria(1000, 0.9, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec09, file = "rep_09.txt", sep = " ") #para guardar los resultados
proc.time() - start.time 

###m??s valores iniciales para hacer el heatmap
repe = 1e6 
## p0 = 0
start.time = proc.time()
rep_trayec0 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec0[i,] = trayectoria(1000, 0, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec0, file = "rep_0.txt", sep = " ") #para guardar los resultados
proc.time() - start.time 


## p0 = 0.2
start.time = proc.time()
rep_trayec02 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec02[i,] = trayectoria(1000, 0.2, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec02, file = "rep_02.txt", sep = " ") #para guardar los resultados
proc.time() - start.time

## p0 = 0.3
start.time = proc.time()
rep_trayec03 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec03[i,] = trayectoria(1000, 0.3, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec03, file = "rep_03.txt", sep = " ") #para guardar los resultados
proc.time() - start.time


## p0 = 0.4
start.time = proc.time()
rep_trayec04 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec04[i,] = trayectoria(1000, 0.4, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec04, file = "rep_04.txt", sep = " ") #para guardar los resultados
proc.time() - start.time

## p0 = 0.6
start.time = proc.time()
rep_trayec06 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec06[i,] = trayectoria(1000, 0.6, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec06, file = "rep_06.txt", sep = " ") #para guardar los resultados
proc.time() - start.time

## p0 = 0.7
start.time = proc.time()
rep_trayec07 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec07[i,] = trayectoria(1000, 0.7, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec07, file = "rep_07.txt", sep = " ") #para guardar los resultados
proc.time() - start.time


## p0 = 0.8
start.time = proc.time()
rep_trayec08 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec08[i,] = trayectoria(1000, 0.8, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec08, file = "rep_08.txt", sep = " ") #para guardar los resultados
proc.time() - start.time

## p0 = 1
start.time = proc.time()
rep_trayec1 = matrix(0, ncol = 11, nrow = repe)
for (i in 1:repe) {
  rep_trayec1[i,] = trayectoria(1000, 1, 500)[c(1,50,100,150,200,250,300,350,400,450,500)]
}
write.table(rep_trayec1, file = "rep_1.txt", sep = " ") #para guardar los resultados
proc.time() - start.time

################# Esto es considerando  t'/N
#rep_trayec01 = read.table(file.choose(), header = T) #rep_01.txt
#rep_trayec05 = read.table(file.choose(), header = T) #rep_05.txt
#rep_trayec09 = read.table(file.choose(), header = T) #rep_09.txt


x01= 0.1
x0 = 0.5
x09= 0.9
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


## LA formula (asymptotical expansion), p0=0.1
expan_01_1 = (1/sqrt(2*pi*(1/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*1)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_50 = (1/sqrt(2*pi*(50/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*50)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_100 = (1/sqrt(2*pi*(100/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*100)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_150 = (1/sqrt(2*pi*(150/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*150)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_200 = (1/sqrt(2*pi*(200/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*200)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_250 = (1/sqrt(2*pi*(250/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*250)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_300 = (1/sqrt(2*pi*(300/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*300)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_350 = (1/sqrt(2*pi*(350/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*350)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_400 = (1/sqrt(2*pi*(400/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*400)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_450 = (1/sqrt(2*pi*(450/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*450)) * (asin(1-2*x01) - asin(1-2*x))^(2))
expan_01_500 = (1/sqrt(2*pi*(500/1000))) * (x01 * (1-x01))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*500)) * (asin(1-2*x01) - asin(1-2*x))^(2))


## LA formula (asymptotical expansion), p0=0.9
expan_09_1 = (1/sqrt(2*pi*(1/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*1)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_50 = (1/sqrt(2*pi*(50/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*50)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_100 = (1/sqrt(2*pi*(100/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*100)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_150 = (1/sqrt(2*pi*(150/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*150)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_200 = (1/sqrt(2*pi*(200/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*200)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_250 = (1/sqrt(2*pi*(250/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*250)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_300 = (1/sqrt(2*pi*(300/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*300)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_350 = (1/sqrt(2*pi*(350/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*350)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_400 = (1/sqrt(2*pi*(400/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*400)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_450 = (1/sqrt(2*pi*(450/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*450)) * (asin(1-2*x09) - asin(1-2*x))^(2))
expan_09_500 = (1/sqrt(2*pi*(500/1000))) * (x09 * (1-x09))^(1/4) / (x*(1-x))^(3/4) * exp(-(1000/(2*500)) * (asin(1-2*x09) - asin(1-2*x))^(2))



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

###LA formula (gaussian approximation), p0=0.1
gauss_app_01_1 = (1/sqrt(2*pi*(1/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*1)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_50 = (1/sqrt(2*pi*(50/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*50)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_100 = (1/sqrt(2*pi*(100/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*100)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_150 = (1/sqrt(2*pi*(150/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*150)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_200 = (1/sqrt(2*pi*(200/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*200)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_250 = (1/sqrt(2*pi*(250/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*250)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_300 = (1/sqrt(2*pi*(300/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*300)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_350 = (1/sqrt(2*pi*(350/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*350)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_400 = (1/sqrt(2*pi*(400/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*400)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_450 = (1/sqrt(2*pi*(450/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*450)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )
gauss_app_01_500 = (1/sqrt(2*pi*(500/1000))) * (1/ sqrt((x01*(1-x01)))) * exp(-(1000/(2*500)) * ((x-x01)^2) / (sqrt(x01*(1-x01))) )

###LA formula (gaussian approximation), p0=0.9
gauss_app_09_1 = (1/sqrt(2*pi*(1/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*1)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_50 = (1/sqrt(2*pi*(50/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*50)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_100 = (1/sqrt(2*pi*(100/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*100)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_150 = (1/sqrt(2*pi*(150/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*150)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_200 = (1/sqrt(2*pi*(200/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*200)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_250 = (1/sqrt(2*pi*(250/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*250)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_300 = (1/sqrt(2*pi*(300/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*300)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_350 = (1/sqrt(2*pi*(350/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*350)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_400 = (1/sqrt(2*pi*(400/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*400)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_450 = (1/sqrt(2*pi*(450/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*450)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )
gauss_app_09_500 = (1/sqrt(2*pi*(500/1000))) * (1/ sqrt((x09*(1-x09)))) * exp(-(1000/(2*500)) * ((x-x09)^2) / (sqrt(x09*(1-x09))) )


######t=1 
#p0=0.1
par(mfrow=c(1,3))
hist(p0_inicial01[,1], main = "", xlab = "p(t)",breaks = 15, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_01_1, col = "green", lwd = 2)
lines(x, gauss_app_1, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((1 / (N)) * (x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,1], main = "", xlab = "p(t)",breaks = 15, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_1, col = "green", lwd = 2)
lines(x, gauss_app_1, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((1 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,1], main = "", xlab = "p(t)",breaks = 15, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_09_1, col = "green", lwd = 2)
lines(x, gauss_app_09_1, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((1 / (N)) * (x09*(1-x09)) )) , col = "red", lwd = 2)



#######t=50
#p0=0.1
hist(p0_inicial01[,2], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_50, col = "green", lwd = 2)
lines(x, gauss_app_01_50, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((50 / (N)) * (x01*(1-x01)) )) , col = "red", lwd = 2)


#p0=0.5
hist(p0_inicial05[,2], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,4),freq = F, col = "navy")
lines(x, expan_50, col = "green", lwd = 2)
lines(x, gauss_app_50, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((50 / (N)) * (x0*x0) )) , col = "red", lwd = 2)


#p0=0.9
hist(p0_inicial09[,2], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_50, col = "green", lwd = 2)
lines(x, gauss_app_09_50, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((50 / (N)) * (x09*(1-x09)) )) , col = "red", lwd = 2)


#####t=100
#p0=0.1
hist(p0_inicial01[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_100, col = "green", lwd = 2)
lines(x, gauss_app_01_100, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((100 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_100, col = "green", lwd = 2)
lines(x, gauss_app_100, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((100 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,3], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_100, col = "green", lwd = 2)
lines(x, gauss_app_09_100, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((100 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)

######t=150
#p0=0.1
hist(p0_inicial01[,4], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_150, col = "green", lwd = 2)
lines(x, gauss_app_01_150, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((150 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,4], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_150, col = "green", lwd = 2)
lines(x, gauss_app_150, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((100 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,4], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_150, col = "green", lwd = 2)
lines(x, gauss_app_09_150, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((150 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=200
#p0=0.1
hist(p0_inicial01[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_200, col = "green", lwd = 2)
lines(x, gauss_app_01_200, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((200 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_200, col = "green", lwd = 2)
lines(x, gauss_app_200, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((200 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,5], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_200, col = "green", lwd = 2)
lines(x, gauss_app_09_200, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((200 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=250
#p0=0.1
hist(p0_inicial01[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_250, col = "green", lwd = 2)
lines(x, gauss_app_01_250, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((250 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_250, col = "green", lwd = 2)
lines(x, gauss_app_250, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((250 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,6], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_250, col = "green", lwd = 2)
lines(x, gauss_app_09_250, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((250 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=300
#p0=0.1
hist(p0_inicial01[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_300, col = "green", lwd = 2)
lines(x, gauss_app_01_300, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((300 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_300, col = "green", lwd = 2)
lines(x, gauss_app_300, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((300 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,7], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_300, col = "green", lwd = 2)
lines(x, gauss_app_09_300, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((300 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=350
#p0=0.1
hist(p0_inicial01[,8], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_350, col = "green", lwd = 2)
lines(x, gauss_app_01_350, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((350 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,8], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_350, col = "green", lwd = 2)
lines(x, gauss_app_350, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((350 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,8], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_350, col = "green", lwd = 2)
lines(x, gauss_app_09_350, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((350 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=400
#p0=0.1
hist(p0_inicial01[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_400, col = "green", lwd = 2)
lines(x, gauss_app_01_400, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((400 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_400, col = "green", lwd = 2)
lines(x, gauss_app_400, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((400 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,9], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_400, col = "green", lwd = 2)
lines(x, gauss_app_09_400, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((400 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=450
#p0=0.1
hist(p0_inicial01[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_450, col = "green", lwd = 2)
lines(x, gauss_app_01_450, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((450 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_450, col = "green", lwd = 2)
lines(x, gauss_app_450, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((450 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,10], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_450, col = "green", lwd = 2)
lines(x, gauss_app_09_450, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((450 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)


######t=500
#p0=0.1
hist(p0_inicial01[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_01_500, col = "green", lwd = 2)
lines(x, gauss_app_01_500, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x01, sd = sqrt((500 / (N)) *(x01*(1-x01)) )) , col = "red", lwd = 2)

#p0=0.5
hist(p0_inicial05[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,3),freq = F, col = "navy")
lines(x, expan_500, col = "green", lwd = 2)
lines(x, gauss_app_500, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x0, sd = sqrt((500 / (N)) * (x0*x0) )) , col = "red", lwd = 2)

#p0=0.9
hist(p0_inicial09[,11], main = "", xlab = "p(t)",breaks = 45, xlim = c(0,1), ylim = c(0,7),freq = F, col = "navy")
lines(x, expan_09_500, col = "green", lwd = 2)
lines(x, gauss_app_09_500, col = "purple", lwd = 2)
lines(x, dnorm(x, mean = x09, sd = sqrt((500 / (N)) *(x09*(1-x09)) )) , col = "red", lwd = 2)

######## Lectura de los datos simulados
p0_inicial01 = read.table(file.choose(), header = T) #rep01.txt
p0_inicial05 = read.table(file.choose(), header = T) #rep05.txt
p0_inicial09 = read.table(file.choose(), header = T) #rep09.txt

#hist(p0_inicial05[,8])
#head(p0_inicial05)

#### construir los kernels (beta y boundary kernel)
##0. definir una malla en el intervalo [0,1]
n = 50
y = seq(0, 1, 1/n)

##1. beta kernel 
# se calcula una matriz para un t fijo (columna de p0_inicial01)
# se calcula el beta kernel (Bertin) para cada elemento de la trayectoria
# y una malla del intervalo [0,1] y un valor de b en especial 
# finalmente se calcula el beta kernel estimator
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

###hay que ordenarlos y asi queda bien el grafico
bk_t50_b09 = beta_kernel(p0_inicial01[,5], y, 0.9); bk_t50_b09
dim(bk_t50_b09)
colMeans(bk_t50_b09)


plot(y,colMeans(bk_t50_b09), type = "l", xlab = "grid")
##2. estimador beta kernel
beta_est = function(x,b,n){
  y = seq(0,1,length.out = n)
  beta_e = rep(NA, length(y))
  
  for (i in 1:length(y)) {
    beta_e[i] = mean(beta_kernel(x, y[i], b))
  }
  return(cbind(beta_e, y))
}

#generacion t=100, b=0.1
beta_e_b01 = beta_est(p0_inicial05[,3], 0.1, 100)

#generacion t=100, b=0.3
beta_e_b03 = beta_est(p0_inicial05[,3], 0.3, 100)

#generacion t=100, b=0.05
beta_e_b005 = beta_est(p0_inicial05[,3], 0.05, 100)

#generacion t=100, b=0.01
beta_e_b001 = beta_est(p0_inicial05[,3], 0.01, 100)

hist(p0_inicial05[,3], breaks = 35, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(beta_e_b01[,2], beta_e_b01[,1], col = "red", lwd=2)
lines(beta_e_b03[,2], beta_e_b03[,1], col = "green", lwd=2)
lines(beta_e_b005[,2], beta_e_b005[,1], col = "blue", lwd=2)
lines(beta_e_b001[,2], beta_e_b001[,1], col = "orange", lwd=2)
legend("topleft", legend=c("b=0.1", "b=0.3","b=0.05","b=0.01"),
       col=c("red", "green","blue","orange"), lty=1, cex=0.8)



######### ADAPTIVE ESTIMATION (para el valor ??ptimo de b) ######
#definir valores iniciales
n = 1e6
varep_1 = 10^(-4) 
gamma_01 = varep_1
gamma_02 = (1/2) * ((1 + varep_1)/(1 - varep_1))

##log en base 10 o log natural? (voy a usar log natural)
Kn = as.integer(log10(n))^2 #36
Kn = as.integer(log(n))^2 #169

k_values = seq(0, Kn, 1)

##gamma_k1 (p=1)  
gamma_k1 = c()
for (i in k_values) {
  gamma_k1[i] = gamma_01 + (k_values[i] / Kn) * (2 - gamma_01)
}
## los valores de b tienen que ser menores a 1
gamma_k1 = gamma_k1[gamma_k1 <= 1]
gamma_k1 = gamma_k1[2:85]

#gamma_k1
range(gamma_k1)


## estimadores con gamma_k1
est_gamma_k1_t1 = matrix(NA, ncol = 84, nrow = 70)
ptm <- proc.time()
for (j in 1:84) {
  est_gamma_k1_t1[,j] = beta_est(p0_inicial05[,1], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t50 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t50[,j] = beta_est(p0_inicial05[,2], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t100 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t100[,j] = beta_est(p0_inicial05[,3], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t150 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t150[,j] = beta_est(p0_inicial05[,4], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t200 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t200[,j] = beta_est(p0_inicial05[,5], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t250 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t250[,j] = beta_est(p0_inicial05[,6], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t300 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t300[,j] = beta_est(p0_inicial05[,7], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t350 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t350[,j] = beta_est(p0_inicial05[,8], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t400 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t400[,j] = beta_est(p0_inicial05[,9], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t450 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t450[,j] = beta_est(p0_inicial05[,10], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

ptm <- proc.time()
est_gamma_k1_t500 = matrix(NA, ncol = 84, nrow = 70)
for (j in 1:84) {
  est_gamma_k1_t500[,j] = beta_est(p0_inicial05[,11], gamma_k1[j], 70)[,1]
}
proc.time() - ptm

write.table(est_gamma_k1_t1, file = "est_gamma_k1_t1.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t50, file = "est_gamma_k1_t50.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t100, file = "est_gamma_k1_t100.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t150, file = "est_gamma_k1_t150.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t200, file = "est_gamma_k1_t200.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t250, file = "est_gamma_k1_t250.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t300, file = "est_gamma_k1_t300.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t350, file = "est_gamma_k1_t350.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t400, file = "est_gamma_k1_t400.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t450, file = "est_gamma_k1_t450.txt", sep = " ") #para guardar los resultados
write.table(est_gamma_k1_t500, file = "est_gamma_k1_t500.txt", sep = " ") #para guardar los resultados


##gamma_k2 (p=2)  
gamma_k2 = c()
for (i in k_values) {
  gamma_k2[i] = gamma_02 + (k_values[i] / Kn) * (2 - gamma_02)
}

#gamma_k2
range(gamma_k2)

##usando la funci??n beta con los valores de gamma_k1 y gamma_k2 (oracle)
x =seq(0,1, length.out =  1000000)
y_t100 = sort(p0_inicial05[,3]) #t=100 (generacion)
plot(x, y_t100, type = "l")

#argmin_{k in Kn} || f_{gamma_k} - f ||_2^2
#f es mis valores de p0??--> son mis densidades aproximadas anteriores

#valores iniciales
x01= 0.1
x0 = 0.5
x09= 0.9
N = 1000
t = 500
x = seq(0,1, length.out = 1000) #para superponer sobre el histograma






### no funciona para la primera componente de gamma_k1[1]

beta_e_gammak1_2 =  beta_est(p0_inicial05[,3], gamma_k1[2], 100)

hist(p0_inicial05[,3], breaks = 35, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
lines(beta_e_gammak1_2[,2], beta_e_gammak1_2[,1], col = "red", lwd=2)

hist(p0_inicial05[,3], breaks = 35, freq = F, ylim = c(0,3.5), xlab = expression(p[0]), main = "")
p05_est_gamma_k1_t100

