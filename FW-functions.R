### Created functions

#1. trayectorias
trayectoria = function(N,p,t){
  frec = rep(0, t)
  frec[1] = p
  
  for (i in 2:t) {
    p = sum( rbinom(n=N, size = 1, prob = p) ) /N
    frec[i] = p
  }
  return(frec)
}


#2. Beta Kernel
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

#3. Beta Kernel Estimator
beta_est = function(x,b,n){
  y = seq(0,1,length.out = n)
  beta_e = rep(NA, length(y))
  
  for (i in 1:length(y)) {
    beta_e[i] = mean(beta_kernel(x, y[i], b))
  }
  return(cbind(beta_e, y))
}


#4. Asymptotical expansion
expan = function(x0, N, t){
  x = seq(0,1, length.out = N) #para superponer sobre el histograma
  a_expan = (1/sqrt(2*pi*(t/(2*N)))) * (x0 * (1-x0))^(1/4) / (x*(1-x))^(3/4) * exp(-((2*N)/(2*t)) * (asin(1-2*x0) - asin(1-2*x))^(2))
  a_expan[1] = a_expan[N] = 0
  return(a_expan)
}


#5. Gaussian approximation
gauss = function(x0, N, t){
  x = seq(0,1, length.out = N) #para superponer sobre el histograma
  gauss_app = (1/sqrt(2*pi*(t/(2*N)))) * (1/ sqrt((x0*(1-x0)))) * exp(-((2*N)/(2*t)) * ((x-x0)^2) / (sqrt(x0*(1-x0))) )
  return(gauss_app)
}


## No he creado una función para el estimador beta adaptive
