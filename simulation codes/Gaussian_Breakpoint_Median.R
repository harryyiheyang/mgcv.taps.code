library(mgcv.taps)
library(mgcv)
library(qgam)
N=c(1:10)*500
A=c(0:9)*0.02
B1=B2=array(0,c(1000,5,10))
filename="~/mgcv_maps_simulation/simulation/Gaussian_Breakpoint_Median.rds"

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
X1_design=smoothed_linearity_discontinuity(x1,list(0.5,a))[,-1]
f1=c(X1_design%*%c(4,4,-0))
}else{
X1_design=linearity_discontinuity(x1,0.5)[,-1]
f1=c(X1_design%*%c(4,4,-0))
}
t2=2*pi*x2
f2=0.4*sin(t2)+0.8*cos(t2)+1.2*sin(t2)^2+1.6*cos(t2)^3+2*sin(t2)^3
t3=2*(x3-0.5)
f3=3*sin(3*t3)+6*exp(-36*t3^2)
f4=0*x4
eta=f1+f2+f3
outlier=eta*0
outlier[sample(n,round(0.01*n))]=sign(rnorm(round(0.01*n),0,5))*10
y=1+eta+rnorm(n,0,1)*sd(eta/2)+outlier
fit1=gam(y~s(x1,bs="AMatern",k=10,m=100,xt=list(getA=linearity_discontinuity,para=0.5))+s(x2,bs="cr",k=10)+s(x3,bs="cr",k=15)+s(x4,bs="cr",k=10),method="REML")
test1=taps_wald_test(fit1,test.component=1)
dat=data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4)
fit2=qgam(y~s(x1,bs="AMatern",k=10,m=100,xt=list(getA=linearity_discontinuity,para=0.5))+s(x2,bs="cr",k=10)+s(x3,bs="cr",k=15)+s(x4,bs="cr",k=10),qu=0.5,data=dat)
test2=taps_wald_test(fit2,test.component=1)
B1[iter,i,j]=ifelse(test1$smooth.pvalue<=0.05,1,0)
B2[iter,i,j]=ifelse(test2$smooth.pvalue<=0.05,1,0)
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
saveRDS(list(B1=B1,B2=B2),filename)
}
}

