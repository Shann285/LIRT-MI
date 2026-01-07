
rm(list=ls()) #clear screen
setwd("C:/Users/dell/Desktop/Data")

###############################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##read data
mydata <- read.csv(file="fy11.csv", header = TRUE)
#head(mydata)
bt<-proc.time()

set.seed(123)
#data 
N=dim(mydata)[1]       #sample size
q=2         #latent intercept/slope
T=3         #time points
I=7        #items

K<-matrix(c(
1,0,
1,1,
1,2),nrow=T,ncol=q,byr=TRUE)


# compile the stan model
model_code1 = "
data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=2> T;
  int<lower=1> q;
  matrix[T,q] K;
  int<lower=0,upper=1> Y[N,T,I];
}

parameters {
  real mua0;
  vector<lower=0>[q] taua;
  cholesky_factor_corr[q] Cova;
  vector[T] psx[N];
  vector[I] A0[T];
  vector<lower=0>[I] L0[T];
  vector<lower=0>[I] Ls;
  vector<lower=0>[I] Lamdai;
  vector<lower=0>[I] Lamdas;
  vector[q+I] fac_dist_helper[N];
}

transformed parameters {
  real<lower=0> tau;
  vector[q] mua;
  vector[T] athre[N];
  vector[q] ath[N];
  vector[I] A00[T-1];
  vector[I] L00[T-1];
  matrix[q,q] Siga;
  mua[1] = 0;
  mua[2] = mua0;
  tau = 1-taua[1]*taua[1];
  Siga = diag_pre_multiply(taua,Cova) * diag_pre_multiply(taua,Cova)';
  for(n in 1:N){
    ath[n] = mua + diag_matrix(taua) * Cova * fac_dist_helper[n,1:q];
    athre[n] = K*ath[n] + psx[n];
  }
  for(t in 2:T){
    A00[t-1] = A0[t]-A0[1];
    L00[t-1] = L0[t]-L0[1];
  }
}

model {
  for(n in 1:N){
    psx[n] ~ normal(0, sqrt(tau));
    for(t in 1:T){
      Y[n,t] ~ bernoulli_logit(A0[t]+athre[n,t]*L0[t]+Ls.*fac_dist_helper[n,(q+1):(q+I)]);   
    }
    fac_dist_helper[n] ~ normal(0,1);
  }
  A0[1] ~ normal(0,2);
  A0[2] ~ double_exponential(A0[1], 1 ./ sqrt(Lamdai) ); 
  A0[3] ~ double_exponential(A0[1], 1 ./ sqrt(Lamdai) ); 
  L0[1] ~ normal(0,2);
  L0[2] ~ double_exponential(L0[1], 1 ./ sqrt(Lamdas) ); 
  L0[3] ~ double_exponential(L0[1], 1 ./ sqrt(Lamdas) ); 
  Lamdai ~ gamma(9, 3);
  Lamdas ~ gamma(9, 3);
  Ls ~ normal(0,2);
  mua0 ~ normal(0,2);
  taua ~ cauchy(0,2.5);
  Cova ~ lkj_corr_cholesky(2);
}
"

##data convert
Y<-array(0, dim=c(N, T, I))

Y[,,1] <- as.matrix(mydata[,c(1,(I+1),(2*I+1))])
Y[,,2] <- as.matrix(mydata[,c(2,(I+2),(2*I+2))])
Y[,,3] <- as.matrix(mydata[,c(3,(I+3),(2*I+3))])
Y[,,4] <- as.matrix(mydata[,c(4,(I+4),(2*I+4))])
Y[,,5] <- as.matrix(mydata[,c(5,(I+5),(2*I+5))])
Y[,,6] <- as.matrix(mydata[,c(6,(I+6),(2*I+6))])
Y[,,7] <- as.matrix(mydata[,c(7,(I+7),(2*I+7))])
  

data0<-list(N=N, I=I, T=T, q=q, K=K, Y=Y)

inits = function() {
  init.values<-list(
    mua0 = runif((q-1),0,.1),
    taua = runif(q,0,0.1),
    Cova = matrix(c(1,0,0,0.9),nrow=q),
    L0 = matrix(runif(I*T,1,1.1),nrow=T),             
    Ls = runif(I,0,.1),
    A0 = matrix(runif(I*T,0,.1),nrow=T)
    )
  return(init.values);
}

irt_1pl <- stan(model_code = model_code1, pars = c("mua0","Siga","A00","A0","L00","L0","Ls"), data=data0, iter=10000, init=inits, chains=3, cores=3)


get_num_divergent(irt_1pl)
summary(irt_1pl)


