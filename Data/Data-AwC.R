
rm(list=ls()) #clear screen
setwd("C:/Users/dell/Desktop/Data")

##
library(lavaan)
library(sirt)

##read data
mydata <- read.csv(file="fy11.csv", header = TRUE)
#dim(mydata)
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
1,2),nrow=T,ncol=q,byr=T)

##data convert
Y<-array(0, dim=c(N, T, I))

Y[,,1] <- as.matrix(mydata[,c(1,(I+1),(2*I+1))])
Y[,,2] <- as.matrix(mydata[,c(2,(I+2),(2*I+2))])
Y[,,3] <- as.matrix(mydata[,c(3,(I+3),(2*I+3))])
Y[,,4] <- as.matrix(mydata[,c(4,(I+4),(2*I+4))])
Y[,,5] <- as.matrix(mydata[,c(5,(I+5),(2*I+5))])
Y[,,6] <- as.matrix(mydata[,c(6,(I+6),(2*I+6))])
Y[,,7] <- as.matrix(mydata[,c(7,(I+7),(2*I+7))])

y11<-Y[,1,1]
y21<-Y[,1,2]
y31<-Y[,1,3]
y41<-Y[,1,4]
y51<-Y[,1,5]
y61<-Y[,1,6]
y71<-Y[,1,7]


y12<-Y[,2,1]
y22<-Y[,2,2]
y32<-Y[,2,3]
y42<-Y[,2,4]
y52<-Y[,2,5]
y62<-Y[,2,6]
y72<-Y[,2,7]


y13<-Y[,3,1]
y23<-Y[,3,2]
y33<-Y[,3,3]
y43<-Y[,3,4]
y53<-Y[,3,5]
y63<-Y[,3,6]
y73<-Y[,3,7]

 
##configural model 
fit_with_cog <- function(n_items, n_waves){  
  yindices <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  ynames <- t(yindices)
  ynames[] <- paste0("y", ynames)
  ynames0 <- yindices
  ynames0[] <- paste0("y", ynames0)
  dp <- t(yindices)
  dp[] <- paste0("d", dp)
  lp <- t(yindices)
  lp[] <- paste0("l", lp)

  load_lines <- apply(`[<-`(ynames, paste0(lp, " * ", ynames)), 1, paste0, collapse = " + ")
  lambda_prefix0 <- matrix(rep(paste0("s",seq_len(n_items)), n_waves),nrow=n_items)
  load_lines0 <- apply(`[<-`(ynames0, paste0(lambda_prefix0, " * ", ynames0)), 1, paste0, collapse = " + ")
  conf_model <- paste0(
    c(
     paste0("theta", seq_len(n_waves), " =~ NA * ", ynames[ , 1], " + ",
             load_lines,
             collapse = "\n  "),
     paste0("c", seq_len(n_items), " =~ NA * ", ynames0[ , 1], " + ",
             load_lines0,
             collapse = "\n  "),
    paste0("theta", seq_len(n_waves), " ~~ 1*", "theta", seq_len(n_waves), collapse = "\n  "),
    paste0("c", seq_len(n_items), " ~~ 1*", "c", seq_len(n_items), collapse = "\n  "), 
    paste0(combn(paste0("c", seq_len(n_items)),2,FUN = paste0, collapse = " ~~ 0*"), collapse = "\n  "), 
    paste0(outer(paste0("c", seq_len(n_items), " ~~ 0*"), paste0("theta", seq_len(n_waves)),  FUN= paste0, sep = ""), collapse = "\n  "),
    paste0("c", seq_len(n_items), " ~ 0*", "1", collapse = "\n  "),
    paste0("theta", seq_len(n_waves), " ~ 0*", "1", collapse = "\n  "),
    paste0(ynames, " | ", "0*t1", collapse = "\n  "),
    paste0(ynames, " ~ ", dp, "*1", collapse = "\n  ")
   ),
   collapse = "\n  "
  ) 
  conf_model
}


##growth model
fit_with_aliglgm <- function(n_items, n_waves, dp, lp){
  yindices <- outer(seq_len(n_items), seq_len(n_waves), paste0)
  ynames <- t(yindices)
  ynames[] <- paste0("y", ynames)
  ynames0 <- yindices
  ynames0[] <- paste0("y", ynames0)
  
  load_lines <- apply(`[<-`(ynames, paste0(lp, " * ", ynames)), 1, paste0, collapse = " + ")
  lambda_prefix0 <- matrix(rep(paste0("s",seq_len(n_items)), n_waves),nrow=n_items)
  load_lines0 <- apply(`[<-`(ynames0, paste0(lambda_prefix0, " * ", ynames0)), 1, paste0, collapse = " + ")
  conf_model <- paste0(
    c(
     paste0("theta", seq_len(n_waves), " =~ ",
             load_lines,
             collapse = "\n  "),
     paste0("c", seq_len(n_items), " =~ NA * ", ynames0[ , 1], " + ",
             load_lines0,
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
    paste0("int ", "~~ 0*", paste0("c",seq_len(n_items)), collapse = "\n  "),
    paste0("slo ", "~~ 0*", paste0("c",seq_len(n_items)), collapse = "\n  "),
    paste0("c", seq_len(n_items), " ~~ 1*", "c", seq_len(n_items), collapse = "\n  "), 
    paste0(combn(paste0("c", seq_len(n_items)),2,FUN = paste0, collapse = " ~~ 0*"), collapse = "\n  "), 
    paste0(outer(paste0("c", seq_len(n_items), " ~~ 0*"), paste0("theta", seq_len(n_waves)),  FUN= paste0, sep = ""), collapse = "\n  "),
    paste0("c", seq_len(n_items), " ~ 0*", "1", collapse = "\n  "),
    paste0(ynames, " | ", "0*t1", collapse = "\n  "),
    paste0(ynames, " ~ ", dp, "*1", collapse = "\n  ")
   ),
   collapse = "\n  "
  ) 
  conf_model
}


##lavaan fit
data1<-cbind(y11,y21,y31,y41,y51,y61,y71,y12,y22,y32,y42,y52,y62,y72,y13,y23,y33,y43,y53,y63,y73)

fit2 <- sem(fit_with_cog(I,T), data = data1, ordered = c("y11","y21","y31","y41","y51","y61","y71","y12","y22","y32","y42","y52","y62","y72",
          "y13","y23","y33","y43","y53","y63","y73"), parameterization = "theta")

inspect(fit2, what = "converged")

fitMeasures(fit2, c("chisq", "df", "pvalue", "cfi", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper"))


# Extract loadings and intercepts for alignment
lam_mat <- (lavInspect(fit2, what = "est")$lambda)[,1:T]
nu_vec <- lavInspect(fit2, what = "est")$nu
# Put them into T x p matrices

lam_config <- crossprod(lam_mat, rep(1, T) %x% diag(I))
nu_config <- matrix(nu_vec, nrow = T, ncol = I, byrow = TRUE)
# Add indicator names
colnames(lam_config) <- colnames(nu_config) <- c("Y1", "Y2" , "Y3", "Y4", "Y5", "Y6", "Y7")

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


fit3 <- sem(fit_with_aliglgm(I,T,dp,lp), data = data1, ordered = c("y11","y21","y31","y41","y51","y61","y71","y12","y22","y32","y42","y52","y62","y72",
          "y13","y23","y33","y43","y53","y63","y73"), parameterization = "theta")


inspect(fit3, what = "converged")

fitMeasures(fit3, c("chisq", "df", "pvalue", "cfi", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper"))

parameterEstimates(fit3)


