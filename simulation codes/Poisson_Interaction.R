library(mgcv.taps)
library(mgcv)
N=c(1:10)*500
A=c(0:9)
B1=B2=array(0,c(1000,5,10))
Pvalue1=Pvalue2=matrix(0,1000,10)
filename="~/mgcv_maps_simulation/simulation/Poisson_Interaction.rds"

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
f12=x1+x2+x1*x2+a*exp(-abs(x1-x2))
t3=2*pi*x3
f3=0.4*sin(t3)+0.8*cos(t3)+1.2*sin(t3)^2+1.6*cos(t3)^3+2*sin(t3)^3
t4=2*(x4-0.5)
f4=4*sin(3*t4)+6*exp(-36*t4^2)
eta=f12+f3+f4
eta=eta/sd(eta)
y=rpois(n,exp(eta))
fit=gam(y~s(x1,x2,bs="A2Matern",k=10,m=100)+s(x3,bs="cr",k=10)+s(x4,bs="cr",k=15),method="REML",family=poisson())
test1=taps_wald_test(fit,test.component=1)
test2=taps_score_test(fit,test.component=1)
B1[iter,i,j]=ifelse(test1$smooth.pvalue<=0.05,1,0)
B2[iter,i,j]=ifelse(test2$smooth.pvalue<=0.05,1,0)
if(i==1){
  Pvalue1[iter,j]=test1$smooth.pvalue
  Pvalue2[iter,j]=test2$smooth.pvalue
}
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
print(apply(B1,c(2,3),mean))
print(apply(B2,c(2,3),mean))
saveRDS(list(B1=B1,B2=B2,P1=Pvalue1,P2=Pvalue2),filename)
}
}

