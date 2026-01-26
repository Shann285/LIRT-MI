
rm(list=ls()) #clear screen
setwd("C:/Users/dell/Desktop/NU200S3M")

####
library(lavaan)
library(sirt)


set.seed(123)
#data 
N=200       #sample size
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
1,1.25-0.2,0.95-0.2,0.85-0.2,1.0,1.15,1.05,0.75,
1,1.25-0.2,0.95-0.2,0.85-0.2,1.0,1.15,1.05,0.75
),nrow=T,ncol=I,byr=T)

A<-matrix(c(
0,0.25,-0.15,-0.1,0.1,0.15,-0.15,-0.1,
0,0.25,-0.15,-0.1,0.1,0.15,-0.15,-0.1,
0,0.25+0.3,-0.15+0.3,-0.1+0.3,0.1,0.15,-0.15,-0.1,
0,0.25+0.3,-0.15+0.3,-0.1+0.3,0.1,0.15,-0.15,-0.1
),nrow=T,ncol=I,byr=T)

psx<-0.5

### configural function
fit_with_cog <- function(n_items, n_waves){  
  yindices <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  ynames <- t(yindices)
  ynames[] <- paste0("y", ynames)
  dp <- t(yindices)
  dp[] <- paste0("d", dp)
  lp <- t(yindices)
  lp[] <- paste0("l", lp)

  load_lines <- apply(`[<-`(ynames, paste0(lp, " * ", ynames)), 1, paste0, collapse = " + ")
  conf_model <- paste0(
    c(
     paste0("theta", seq_len(n_waves), " =~ NA * ", ynames[ , 1], " + ",
             load_lines,
             collapse = "\n  "),
    paste0("theta", seq_len(n_waves), " ~~ 1*", "theta", seq_len(n_waves), collapse = "\n  "),
    paste0("theta", seq_len(n_waves), " ~ 0*", "1", collapse = "\n  "),
    paste0(ynames, " | ", "0*t1", collapse = "\n  "),
    paste0(ynames, " ~ ", dp, "*1", collapse = "\n  ")
   ),
   collapse = "\n  "
  ) 
  conf_model
}

## growth function
fit_with_aliglgm <- function(n_items, n_waves, dp, lp){
  yindices <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  ynames <- t(yindices)
  ynames[] <- paste0("y", ynames)
  
  load_lines <- apply(`[<-`(ynames, paste0(lp, " * ", ynames)), 1, paste0, collapse = " + ")
  conf_model <- paste0(
    c(
     paste0("theta", seq_len(n_waves), " =~ ",
             load_lines,
             collapse = "\n  "),
    paste0("int =~ ", paste0("1 * theta", seq_len(n_waves),
                             collapse = " + ")),
    paste0("slo =~ ", paste0(seq_len(n_waves) - 1, " * theta",
                             seq_len(n_waves), collapse = " + ")),
    paste0("int ", "~~ a*", "int"),
    paste0("slo ", "~~ ", "slo"),
    paste0("int ", "~~ ", "slo"),
    paste0("int ", "~ ", "1"),
    paste0("slo ", "~ ", "1"),
    paste0("theta", seq_len(n_waves), " ~~ b*", "theta", seq_len(n_waves),collapse = "\n  "),
    paste0(ynames, " | ", "0*t1", collapse = "\n  "),
    paste0(ynames, " ~ ", dp, "*1", collapse = "\n  ")
   ),
   collapse = "\n  "
  ) 
  conf_model
}


## 50 replications
for(CIR in 1:CNUM){
  bt<-proc.time()
  Y<-array(0, dim=c(N, T, I))
  Y<-array(scan(YFILE, skip=(CIR-1)*T*I, nlines=T*I), dim=c(N, T, I))

y11<-Y[,1,1]
y21<-Y[,1,2]
y31<-Y[,1,3]
y41<-Y[,1,4]
y51<-Y[,1,5]
y61<-Y[,1,6]
y71<-Y[,1,7]
y81<-Y[,1,8]

y12<-Y[,2,1]
y22<-Y[,2,2]
y32<-Y[,2,3]
y42<-Y[,2,4]
y52<-Y[,2,5]
y62<-Y[,2,6]
y72<-Y[,2,7]
y82<-Y[,2,8]

y13<-Y[,3,1]
y23<-Y[,3,2]
y33<-Y[,3,3]
y43<-Y[,3,4]
y53<-Y[,3,5]
y63<-Y[,3,6]
y73<-Y[,3,7]
y83<-Y[,3,8]

y14<-Y[,4,1]
y24<-Y[,4,2]
y34<-Y[,4,3]
y44<-Y[,4,4]
y54<-Y[,4,5]
y64<-Y[,4,6]
y74<-Y[,4,7]
y84<-Y[,4,8]


data1<-cbind(y11,y21,y31,y41,y51,y61,y71,y81,y12,y22,y32,y42,y52,y62,y72,y82,y13,y23,y33,y43,y53,y63,y73,y83,y14,y24,y34,y44,y54,y64,y74,y84)


fit2 <- sem(fit_with_cog(I,T), data = data1, ordered = c("y11","y21","y31","y41","y51","y61","y71","y81","y12","y22","y32","y42","y52","y62","y72","y82",
          "y13","y23","y33","y43","y53","y63","y73","y83","y14","y24","y34","y44","y54","y64","y74","y84"), parameterization = "theta")

Converg <- inspect(fit2, what = "converged")

fit2ind <- fitMeasures(fit2, c("chisq", "df", "pvalue", "cfi", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper"))
write(fit2ind, file="fit2indM2.txt", ncol=length(fit2ind), append=TRUE, sep="\t")

# Extract loadings and intercepts for alignment
lam_mat <- (lavInspect(fit2, what = "est")$lambda)[,1:T]
nu_vec <- lavInspect(fit2, what = "est")$nu
# Put them into T x p matrices

lam_config <- crossprod(lam_mat, rep(1, T) %x% diag(I))
nu_config <- matrix(nu_vec, nrow = T, ncol = I, byrow = TRUE)
# Add indicator names
colnames(lam_config) <- colnames(nu_config) <- c("Y1", "Y2" , "Y3", "Y4", "Y5", "Y6", "Y7", "Y8")

# Alignment optimization
aligned_pars <- sirt::invariance.alignment(
lambda = lam_config,
nu = nu_config,
fixed = TRUE
)


ind <- 1   ###which(colMeans(aligned_pars$lambda.aligned) == max(colMeans(aligned_pars$lambda.aligned))) 

  yindices <- outer(seq_len(I), seq_len(T), paste0)
  dp <- t(yindices)
  dp[] <- paste0("d", dp)
  lp <- t(yindices)
  lp[] <- paste0("l", lp)

lp[,ind] <- aligned_pars$lambda.aligned[,ind]

dp[,ind] <- aligned_pars$nu.aligned[,ind]


fit3 <- sem(fit_with_aliglgm(I,T,dp,lp), data = data1, ordered = c("y11","y21","y31","y41","y51","y61","y71","y81","y12","y22","y32","y42","y52","y62","y72","y82",
          "y13","y23","y33","y43","y53","y63","y73","y83","y14","y24","y34","y44","y54","y64","y74","y84"), parameterization = "theta")


fit3ind <- fitMeasures(fit3, c("chisq", "df", "pvalue", "cfi", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper"))
write(fit3ind, file="fit3indM2.txt", ncol=length(fit3ind), append=TRUE, sep="\t")

Converg <- c(Converg, inspect(fit3, what = "converged"))

estr <- parameterEstimates(fit3)$est[c(44:45, 41:43, 82:113, 1:32, 46 )]
write(estr, file="estrM2.txt", ncol=length(estr), append=TRUE, sep="\t")


ser <- parameterEstimates(fit3)$se[c(44:45, 41:43, 82:113, 1:32, 46 )]
write(ser, file="serM2.txt", ncol=length(ser), append=TRUE, sep="\t")


pvaluer <- parameterEstimates(fit3)$pvalue[c(44:45, 41:43, 82:113, 1:32, 46 )]
write(pvaluer, file="pvaluerM2.txt", ncol=length(pvaluer), append=TRUE, sep="\t")


cil <- parameterEstimates(fit3)$ci.lower[c(44:45, 41:43, 82:113, 1:32, 46 )]
write(cil, file="cilM2.txt", ncol=length(cil), append=TRUE, sep="\t")

ciu <- parameterEstimates(fit3)$ci.upper[c(44:45, 41:43, 82:113, 1:32, 46 )]
write(ciu, file="ciuM2.txt", ncol=length(ciu), append=TRUE, sep="\t")


pva <- NULL

for(item in 2:I){
  for(t in 2:T){
    dptemp <- dp
    dptemp[t,item] <- dptemp[1,item]
    fit3t <- sem(fit_with_aliglgm(I,T,dptemp,lp), data = data1, ordered = c("y11","y21","y31","y41","y51","y61","y71","y81","y12","y22","y32","y42","y52","y62","y72","y82",
          "y13","y23","y33","y43","y53","y63","y73","y83","y14","y24","y34","y44","y54","y64","y74","y84"), parameterization = "theta")
   
    Converg <- c(Converg, inspect(fit3t, what = "converged"))
    resl <- try(lavTestLRT(fit3, fit3t), silent=TRUE)
    if(length(class(resl)) == 1){
      pva <- c(pva, 99)
    }else{
      pva <- c(pva, resl[2,7])
    }
  }
}

write(pva, file="pva3M2.txt", ncol=length(pva), append=TRUE, sep="\t")


pval <- NULL

for(item in 2:I){
  for(t in 2:T){
    lptemp <- lp
    lptemp[t,item] <- lptemp[1,item]
    fit3t <- sem(fit_with_aliglgm(I,T,dp,lptemp), data = data1, ordered = c("y11","y21","y31","y41","y51","y61","y71","y81","y12","y22","y32","y42","y52","y62","y72","y82",
          "y13","y23","y33","y43","y53","y63","y73","y83","y14","y24","y34","y44","y54","y64","y74","y84"), parameterization = "theta")
   
    Converg <- c(Converg, inspect(fit3t, what = "converged"))
    resl <- try(lavTestLRT(fit3, fit3t), silent=TRUE)
    if(length(class(resl)) == 1){
      pval <- c(pval, 99)
    }else{
      pval <- c(pval, resl[2,7])
    }
  }
}


write(pval, file="pva31M2.txt", ncol=length(pval), append=TRUE, sep="\t")


write(Converg, file="ConvergM2.txt", ncol=length(Converg), append=TRUE, sep="\t")

print(CIR)
et<-proc.time()
print((et-bt)[3])
write((et-bt)[3], file="timelavM2.txt", ncol=1, append=TRUE, sep="\t")

}



date()
save.image(paste("RAOlaM2",".RData",sep=""))


#warnings()

Convergr<- rowSums( read.table("ConvergM2.txt") )   

indcat <- Convergr==(23+21)


estres<-matrix(scan("estrM2.txt",nlines=CNUM),nrow=CNUM, byr=T)


truepa<-c(0, 0.2, 0.5, 0.1, 0.045, as.vector(A)/1.7, as.vector(L)/1.7, psx)


bias1 <- apply(estres[indcat,], MARGIN=2, FUN=mean)-truepa


Rmse1 <- colMeans((estres[indcat,] - matrix(rep(truepa, each=length(indcat)), nrow=length(indcat)))^2)


#colMeans(cbind(abs(bias1[38:69]), sqrt(Rmse1[38:69])))
colMeans(cbind(abs(bias1[seq(38,69,4)]), sqrt(Rmse1[seq(38,69,4)]))) ##a1 

colMeans(cbind(abs(bias1[c(39:41,43:45,47,51:53,55:57,59:61,63:65,67:69)]), sqrt(Rmse1[c(39:41,43:45,47,51:53,55:57,59:61,63:65,67:69)]))) ##a2-a4 non-DIF

colMeans(cbind(abs(bias1[c(48,49)]), sqrt(Rmse1[c(48,49)]))) ##a2-a4 DIF


#colMeans(cbind(abs(bias1[6:37]), sqrt(Rmse1[6:37])))

colMeans(cbind(abs(bias1[seq(6,37,4)]), sqrt(Rmse1[seq(6,37,4)])))  ## d1

colMeans(cbind(abs(bias1[c(7:9,11:13,15,19:21,23:25,27:29,31:33,35:37)]), sqrt(Rmse1[c(7:9,11:13,15,19:21,23:25,27:29,31:33,35:37)]))) ##d2-d4 non-DIF

colMeans(cbind(abs(bias1[c(16,17)]), sqrt(Rmse1[c(16,17)]))) ##d2-d4 DIF


cbind(bias1, sqrt(Rmse1))[c(1:5,70),] ##growth parameters



