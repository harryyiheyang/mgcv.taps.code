library(mgcv.taps)
library(mgcv)
N=c(1:10)*500
A=c(0:9)
B=array(0,c(1000,5,10))
filename="~/mgcv_maps_simulation/simulation/Ordinal_Interaction.rds"

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
eta = eta/sqrt(3.2)
u = eta -mean(eta) + rnorm(n, 0, sqrt(2))
alpha = c(-Inf, -3, 3, Inf)
R = length(alpha) - 1  # R=3
y = numeric(n)
for (ii in 1:R) {
  y[u > alpha[ii] & u <= alpha[ii+1]] <- ii
}
fit=gam(y~s(x1,x2,bs="A2Matern",k=10,m=100)+s(x3,bs="cr",k=10)+s(x4,bs="cr",k=15),method="REML",family=ocat(R=R))
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
