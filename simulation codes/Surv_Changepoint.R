library(survival)
library(mgcv.taps)
library(mgcv)
N=c(1:10)*500
A=c(0:4)*0.1
B=B0=array(0,c(1000,5,10))
filename="~/mgcv_maps_simulation/simulation/Surv_Changepoint.rds"

for(i in 1:5){
for(j in 1:10){
cat(paste0("n = ",N[j]," a = ",A[i],"\n"))
pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for(iter in 1:1000){
indicator <- FALSE
tryCatch({
setTxtProgressBar(pb, iter)
n=N[j]
X=MASS::mvrnorm(n,rep(0,4),matrix(0.25,4,4)+0.75*diag(4))
x1=qbeta(pnorm(X[,1]),1.5,1.5)
x2=qbeta(pnorm(X[,2]),1.5,1.5)
x3=qbeta(pnorm(X[,3]),1.5,1.5)
x4=qbeta(pnorm(X[,4]),1.5,1.5)
a=A[i]
if(a!=0){
X1_design=smoothed_piecewise_linearity(x1,list(0.5,a))[,-1]
}else{
X1_design=piecewise_linearity(x1,0.5)[,-1]
}
f1=c(X1_design%*%c(6,-18))
t2=2*pi*x2
f2=0.4*sin(t2)+0.8*cos(t2)+1.2*sin(t2)^2+1.6*cos(t2)^3+2*sin(t2)^3
t3=2*(x3-0.5)
f3=3*sin(3*t3)+6*exp(-36*t3^2)
f4=0*x4
eta=f1+f2+f3
lambda <- 0.1
shape <- 2
scale <- 1 / (lambda*exp(eta-2+rnorm(n,0,sqrt(2))))
time <- rweibull(n, shape = shape, scale = scale)

censor_time <- rexp(n, rate = 0.05)

observed_time <- pmin(time, censor_time)
status <- as.numeric(time <= censor_time)

data <- data.frame(
time = observed_time,
status = status,
x1 = x1,
x2 = x2,
x3 = x3,
x4 = x4
)

fit <- gam(Surv(time, status) ~ s(x1,bs="AMatern",k=3+i,m=100,xt=list(getA=piecewise_linearity,para=0.5))
           +s(x2, bs="cr", k=10)+s(x3, bs="cr", k=15)+s(x4, bs="cr", k=10),
           data = data,
           method="REML",
           family = cox.ph())

test1=taps_wald_test(fit,test.component=1)
B[iter,i,j]=ifelse(test$smooth.pvalue<=0.05,1,0)

}, error = function(e) {
# Error handling block
cat("Error occurred: ", e$message, "\n")
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
iter <<- iter - 1  # Decrement the iteration counter to retry
})
if (indicator) {
next  # Retry the current iteration
}
}
close(pb)
print(apply(B,c(2,3),mean))
saveRDS(B,filename)
}
}
B1=apply(B,c(2,3),mean)
