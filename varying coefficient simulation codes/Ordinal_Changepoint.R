library(mgcv.taps)
library(mgcv)
N=c(1:10)*500
A=c(0:4)*0.15
B=array(0,c(1000,5,10))
filename="~/mgcv_taps_varying/simulation/Ordinal_Changepoint.rds"

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
Z=MASS::mvrnorm(n,rep(0,4),matrix(0.25,4,4)+0.75*diag(4))
x1=qbeta(pnorm(X[,1]),1.5,1.5)
x2=qbeta(pnorm(X[,2]),1.5,1.5)
x3=qbeta(pnorm(X[,3]),1.5,1.5)
x4=qbeta(pnorm(X[,4]),1.5,1.5)
z1=qbeta(pnorm(Z[,1]),1.5,1.5)
z2=qbeta(pnorm(Z[,2]),1.5,1.5)
z3=qbeta(pnorm(Z[,3]),1.5,1.5)
z4=qbeta(pnorm(Z[,4]),1.5,1.5)
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
eta = (f1*z1 + f2*z2 + f3*z3)/sqrt(3.2)
u = eta -mean(eta) + rnorm(n, 0, sqrt(2))
alpha = c(-Inf, -3, 3, Inf)
R = length(alpha) - 1  # R=3
y = numeric(n)
for (ii in 1:R) {
y[u > alpha[ii] & u <= alpha[ii+1]] <- ii
}

fit = gam(y ~ s(x1,by=z1, bs="AMatern", k=10, m=100,
              xt=list(getA=piecewise_linearity, para=0.5)) +
            s(x2,by=z2, bs="cr", k=10) +
            s(x3,by=z3, bs="cr", k=15) +
            s(x4,by=z4, bs="cr", k=10),
          method="REML", family=ocat(R=R))

test=taps_wald_test(fit)
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

