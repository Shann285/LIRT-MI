
rm(list=ls()) #clear screen
setwd("C:/Users/dell/Desktop/NU500L1H")


###############################
library(MASS)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(123)
#data 
N=500       #sample size
q=2         #latent intercept/slope
T=4         #time points
I=8        #items
CNUM=50
YFILE="Y.txt"

###true parameter values
##correlation=0.201 
siga<-matrix(c(
0.5,0.045,
0.045,0.1),nrow=q,ncol=q)

K<-matrix(c(
1,0,
1,1,
1,2,
1,3),nrow=T,ncol=q,byr=T)

L<-matrix(c(
1,1.25,0.95,0.85,1.0,1.15,1.05,0.75,
1,1.25,0.95,0.85,1.0,1.15,1.05,0.75,
1,1.25,0.95-0.4,0.85,1.0,1.15,1.05,0.75,
1,1.25,0.95-0.4,0.85,1.0,1.15,1.05,0.75
),nrow=T,ncol=I,byr=T)

A<-matrix(c(
0,0.25,-0.15,-0.1,0.1,0.15,-0.15,-0.1,
0,0.25,-0.15,-0.1,0.1,0.15,-0.15,-0.1,
0,0.25,-0.15+0.6,-0.1,0.1,0.15,-0.15,-0.1,
0,0.25,-0.15+0.6,-0.1,0.1,0.15,-0.15,-0.1
),nrow=T,ncol=I,byr=T)

psx<-0.5

## simulated data
for(CIR in 1:CNUM){
  aa0<-mvrnorm(N, c(0,0.2), Sigma=siga )
  the0<-matrix(0,nrow=N,ncol=T)
  sf0<-mvrnorm(N, rep(0,I), Sigma=diag(1,I))
  Y<-array(0, dim=c(N,T,I))
  for(n in 1:N){
    the0[n,1]= t(K[1,])%*%aa0[n,] + rnorm(1, 0, sqrt(psx))
    the0[n,2]= t(K[2,])%*%aa0[n,] + rnorm(1, 0, sqrt(psx))
    the0[n,3]= t(K[3,])%*%aa0[n,] + rnorm(1, 0, sqrt(psx))
    the0[n,4]= t(K[4,])%*%aa0[n,] + rnorm(1, 0, sqrt(psx))
    for(t in 1:T){
      for(i in 1:I){
        temp=exp(A[t,i]+the0[n,t]*L[t,i]+1.0*sf0[n,i])/(1+exp(A[t,i]+the0[n,t]*L[t,i]+1.0*sf0[n,i]))
        Y[n,t,i]=rbinom(1,1,temp)
      }
    }
  }
write(Y, file="Y.txt", ncol=dim(Y)[1], append=T)
print(CIR)
}

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
  A0[4] ~ double_exponential(A0[1], 1 ./ sqrt(Lamdai) ); 
  L0[1] ~ normal(0,2);
  L0[2] ~ double_exponential(L0[1], 1 ./ sqrt(Lamdas) ); 
  L0[3] ~ double_exponential(L0[1], 1 ./ sqrt(Lamdas) ); 
  L0[4] ~ double_exponential(L0[1], 1 ./ sqrt(Lamdas) ); 
  Lamdai ~ gamma(9, 3);
  Lamdas ~ gamma(9, 3);
  Ls ~ normal(0,2);
  mua0 ~ normal(0,2);
  taua ~ cauchy(0,2.5);
  Cova ~ lkj_corr_cholesky(2);
}
"

## 50 replications
for(CIR in 1:CNUM){
  bt<-proc.time()
  Y<-array(0, dim=c(N, T, I))
  Y<-array(scan(YFILE, skip=(CIR-1)*T*I, nlines=T*I), dim=c(N, T, I))

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

irt_1pl<-stan(model_code = model_code1, pars = c("mua0","Siga","A00","A0","L00","L0","Ls"), data=data0, iter=8000, init=inits, chains=3, cores=3)

aa<-summary(irt_1pl, probs = c(0.025, 0.975), pars = c("mua0","Siga","A00","A0","L00","L0","Ls"))

write(get_num_divergent(irt_1pl), file="diver.txt", ncol=1, append=TRUE, sep="\t")

write(aa$summary[,c(7)], file="Rhh.txt", ncol=length(aa$summary[,c(7)]), append=TRUE, sep="\t")
write(aa$summary[,c(1)], file="amean.txt", ncol=length(aa$summary[,c(1)]), append=TRUE, sep="\t")
write(aa$summary[,c(3)], file="asd.txt", ncol=length(aa$summary[,c(3)]), append=TRUE, sep="\t")
write(aa$summary[,c(4)], file="ainl.txt", ncol=length(aa$summary[,c(4)]), append=TRUE, sep="\t")
write(aa$summary[,c(5)], file="ainr.txt", ncol=length(aa$summary[,c(5)]), append=TRUE, sep="\t")

print(CIR)
et<-proc.time()
print((et-bt)[3])
write((et-bt)[3], file="time.txt", ncol=1, append=TRUE, sep="\t")
}

date()
save.image(paste("Rthres",".RData",sep=""))


aR<-matrix(scan("Rhh.txt",nlines=CNUM),nrow=CNUM, byr=T)

ad<-matrix(scan("diver.txt",nlines=CNUM),nrow=CNUM, byr=T)

aain<-which( (rowSums(aR[,c(1:3,5, 30:61, 86:117,118:125)]>1.05)==0) & (ad==0) )

length(aain)/CNUM  ##convergence rate



amean<-matrix(scan("amean.txt",nlines=CNUM),nrow=CNUM, byr=T)

truepa<-c(0.2, 0.5,0.045,0.1, as.vector(t(A)), as.vector(t(L)), rep(1.0,I))


bias1 <- apply(amean[aain, c(1:3,5, 30:61, 86:117,118:125)], MARGIN=2, FUN=mean)-truepa


Rmse1 <- colMeans((amean[aain, c(1:3,5, 30:61, 86:117,118:125)] - matrix(rep(truepa, each=length(aain)), nrow=length(aain)))^2)


colMeans(abs(cbind(bias1, sqrt(Rmse1))[(4*I+5):(5*I+4),])) ##a1 

colMeans(abs(cbind(bias1, sqrt(Rmse1))[c((5*I+5):(6*I+4), (6*I+5),(6*I+6),(6*I+8):(7*I+4), (7*I+5),(7*I+6),(7*I+8):(8*I+4) ),])) ##a2-a4 non-DIF increments

colMeans(abs(cbind(bias1, sqrt(Rmse1))[c((6*I+7),(7*I+7)),])) ##a2-a4 DIF increments


colMeans(abs(cbind(bias1, sqrt(Rmse1))[5:(I+4),])) ## d1

colMeans(abs(cbind(bias1, sqrt(Rmse1))[c((I+5):(2*I+4), (2*I+5),(2*I+6),(2*I+8):(3*I+4), (3*I+5),(3*I+6),(3*I+8):(4*I+4) ),])) ##d2-d4 non-DIF increments

colMeans(abs(cbind(bias1, sqrt(Rmse1))[c((2*I+7),(3*I+7)),])) ##d2-d4 DIF increments


colMeans(abs(cbind(bias1, sqrt(Rmse1))[(8*I+5):(9*I+4),])) ##specifc slopes


cbind(bias1, sqrt(Rmse1))[1:4,] ##growth parameters
