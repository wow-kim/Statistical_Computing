## r 함수

# Bernoulli r.v.
set.seed(23207)
guesses = runif(20)
correct.answers = (guesses < 0.2)
correct.answers
table(correct.answers)

# Binomial r.v.
dbinom(x=4, size=6, prob=0.5)		# Pr(X=4)
pbinom(4, 6, 0.5)				# Pr(X <= 4)
qbinom(0.89, 6, 0.5)			# 89 percentile

defectives = rbinom(24, 15, 0.1)
defectives
any(defectives > 5)

# Poisson r.v.

dpois(x=3, lambda=0.5)
rpois(10, 3.7)

# exponential r.v.
pexp(1, rate = 3)

# normal r.v.
qnorm(0.95, mean=2.7, sd=3.3)
rnorm(10,-3,0.5)
# x ~ N(0,1) conditional on 0<x<3
x = rnorm(10000)
x = x[(0<x) & (x<3)]
hist(x, probability=T)

# multivariate normal r.v.
library(MASS)
mu = c(0, 1)
sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2, byrow = T)
mvrnorm(10, mu, sigma)


## 일반적인 난수 발생법
# 역변환법
rexponential = function(n, lambda)
{
	if (lambda <= 0) stop("lambda should be positive")
	return(-log(runif(n))/lambda)
}

rBernoulli = function(n, p)
{
	if (p<0 || p>1) stop("p should be in [0,1]")
	q = 1 - p
	return(ceiling(runif(n)-q))
}

# 기각법
rnormal = function(n, mu=0, std=1)
{
	if (std<=0) stop("std should be positive")
	r = rep(0, n)
	for (i in 1:n) 
	{
		repeat
		{
			u = runif(2)
			u = - log(u)
			if (u[2] > (u[1]-1)^2)
			{
				r[i] = u[1]
				if (runif(1) < 0.5) r[i] = -r[i]
				break
			}
		}
	}
	r = std*r + mu
	return(r)
}


rbinomial.rejection = function(n, size, prob)
{
	if (prob<0 || prob>1) stop("prob must be in [0,1]")
	p = ifelse(prob<=0.5, prob, 1-prob)
	a = sqrt(2*size*p*(1-p))
	x0 = size*p
	rbinom = rep(0, n)
	for (i in 1:n)
	{
		repeat
		{
			u12 = runif(2)
			yx0 = a * tan(pi*u12[1])
			y = yx0 + x0
			gx = 1.2*a*(1+yx0^2/(a^2))*dbinom(as.integer(y), size, p)
			if (u12[2] <= gx) break
		}
		rbinom[i] = as.integer(y)
	}
	rbinom = ifelse(prob<=0.5, rbinom, size - rbinom)
	
	return(rbinom)
}

## 이산 확률난수의 발생
rdiscunif = function(n, a, b)
{
	if (as.integer(a) == n && as.integer(b) == b)
		return(as.integer(a+(b-a+1)*runif(n)))
	else stop("a and b should be integers.")
}

rbinomial = function(n, size, prob)
{
	if (prob>1 || prob<0) stop("p should be in [0,1]")
	p = ifelse(prob<=0.5, prob, 1-prob)
	f0 = (1-p)^size
	p1p = p/(1-p)
	ranbin = runif(n)
	for (i in 1:n)
	{
		x = 0
		fx = f0
		repeat
		{
			if (ranbin[i] <= fx)
			{
				ranbin[i] = x
				break
			} else
			{
				ranbin[i] = ranbin[i] - fx
				fx = (size-x)/(x+1)*p1p*fx
				x = x + 1
			}
		}
	}
	ranbin = ifelse(prob<=0.5, ranbin, size-ranbin)
	return(ranbin)
}

rpoisson = function(n, lambda)
{
	ep = exp(-lambda)
	rpoiss = rep(0, n)
	for (i in 1:n)
	{
		tr = 1
		repeat
		{
			tr = tr*runif(1)
			if (tr <= ep) break
			else rpoiss[i] = rpoiss[i] + 1
		}
	}
	return(rpoiss)
}

# 다변량 정규난수: 분광분해법
mu <- c(0,0)
Sigma <- matrix(c(1,.9,.9,1), nrow=2, ncol=2)

rmvn.eigen <- function(n,mu,Sigma) {
	d <- length(mu)
	ev <- eigen(Sigma, symmetric=TRUE)
    	lambda <- ev$values
    	V <- ev$vectors
    	R <- V%*%diag(sqrt(lambda))%*%t(V)
    	Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    	X <- Z%*%R+matrix(mu, n, d, byrow=TRUE)
    	X
}

X <- rmvn.eigen(1000,mu,Sigma)
plot(X, xlab="x", ylab="y", pch=20)
print(colMeans(X))
print(cor(X))

# 다변량 정규난수: 촐레스키 분해법
rmvn.Choleski <- function(n,mu,Sigma) {
   	d <- length(mu)
 	Q <- chol(Sigma) # Choleski factorization of Sigma
     	Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
     	X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
     	X
}

y <- subset(x=iris, Species=="virginica")[, 1:4]
mu <- colMeans(y)
Sigma <- cov(y)
mu
Sigma

X <- rmvn.Choleski(200, mu, Sigma)
pairs(X)
 
# 다변량 정규난수 생성기 비교
library(MASS)
library(mvtnorm)
n <- 100            # sample size
d <- 30             # dimension
N <- 2000           # iterations
mu <- numeric(d)
set.seed(100)
system.time(for (i in 1:N)
			rmvn.eigen(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
			rmvn.Choleski(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
			mvrnorm(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
			rmvnorm(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
			cov(matrix(rnorm(n*d), n, d)))
 