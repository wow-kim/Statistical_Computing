## E|X_1 - X_2|의 추정, X_1, X_2 ~ N(0,1)
m <- 1000
g <- numeric(m)
for (i in 1:m) { 
	x <- rnorm(2)
	g[i] <- abs(x[1]-x[2])
}
hist(g, prob = TRUE)
est <- mean(g)
est


## 절사평균의 MSE 
n <- 20
m <- 1000
tmean <- numeric(m)
for (i in 1:m) {
	x <- sort(rnorm(n))
      tmean[i] <- sum(x[2:(n-1)])/(n-2)
}
g <- tmean^2
mse <- mean(g)
se <- sqrt(sum((g-mean(g))^2))/m # se of mse
mse
se

 
# 중앙값의 MSE
n <- 20
m <- 1000
tmean <- numeric(m)
for (i in 1:m) {
	x <- sort(rnorm(n))
	tmean[i] <- median(x)
}
g <- tmean^2
mse <- mean(g)
se <- sqrt(sum((g-mean(g))^2))/m  # se of mse
mse
se

# 오염된 정규분포에서 k차 절사평균의 MSE
n <- 20
K <- n/2-1
m <- 1000
mse <- matrix(0,n/2,6)
trimmed.mse <- function(n, m, k, p) {
	tmean <- numeric(m)
	for (i in 1:m) {
		sigma <- sample(c(1,10), size = n,
				replace = TRUE, prob = c(p,1-p))
		x <- sort(rnorm(n, 0, sigma))
		tmean[i] <- sum(x[(k+1):(n-k)])/(n-2*k)
	}
	g <- tmean^2
	mse.est <- mean(g)
	se.mse <- sqrt(mean((g-mean(g))^2))/sqrt(m)
	return(c(mse.est, se.mse))
}

for (k in 0:K) {
	mse[k+1, 1:2] <- trimmed.mse(n=n, m=m, k=k, p=1.0)
	mse[k+1, 3:4] <- trimmed.mse(n=n, m=m, k=k, p=.95)
	mse[k+1, 5:6] <- trimmed.mse(n=n, m=m, k=k, p=.9)
}


## N(mu, 1)에서 평균의 신뢰구간과 신뢰수준
n <- 20
alpha <- .05
x <- rnorm(n, mean=3, sd=1)
UCL <- mean(x)-qt(1-alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)
LCL <- mean(x)+qt(1-alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)
c(LCL, UCL)
 
n <- 20
alpha <- .05
CL <- replicate(1000, expr = {
		x <- rnorm(n, mean=3, sd=1)
		LCL <- mean(x)-qt(1-alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)
		UCL <- mean(x)+qt(1-alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)
		c(LCL,UCL)
	} )
sum(CL[1,] < 3 & CL[2,] > 3)
mean(CL[1,] < 3 & CL[2,] > 3)

 
## 경험적 1종의 오류율 계산
# 정규분포의 랜덤표본에 대하여 alpha=0.05일 때 H_0: mu=500 vs H_1: mu>500에 대한 t-검정
n <- 20
alpha <- .05
mu0 <- 500
sigma <- 100
m <- 10000          #number of replicates
p <- numeric(m)     #storage for p-values
for (j in 1:m)
{
	x <- rnorm(n,mu0,sigma)
	ttest <- t.test(x, alternative="greater", mu=mu0)
	p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat,se.hat))


## 경험적 검정력 추정
n <- 20
m <- 1000
mu0 <- 500
sigma <- 100
mu <- c(seq(450, 650, 10))    # alternatives
M <- length(mu)
power <- numeric(M)
for (i in 1:M) {
	mu1 <- mu[i]
	pvalues <- replicate(m, expr={
	# simulate under alternative mu1
	x <- rnorm(n, mean=mu1, sd=sigma)
	ttest <- t.test(x,alternative="greater", mu=mu0)
      ttest$p.value} )
   	power[i] <- mean(pvalues <= .05)
}

library(Hmisc)   # for errbar
plot(mu, power)
abline(v=mu0, lty=1)
abline(h=.05, lty=1)
# add standard errors
se <- sqrt(power*(1-power)/m)
errbar(mu, power, yplus=power+se, yminus=power-se, xlab=bquote(theta))
lines(mu, power, lty=3)
detach(package:Hmisc)