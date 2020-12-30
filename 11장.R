## Rayliegh(sigma) 분포에 대한 M-H 표집기
# 제안분포는 chi^2(X_t)

f <- function(x, sigma) {
	if (any(x < 0)) return (0)
	stopifnot(sigma > 0)
	return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

m <- 10000
sigma <- 4
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)

for (i in 2:m) {
	xt <- x[i-1]
	y <- rchisq(1, df = xt)
      num <- f(y, sigma) * dchisq(xt, df = y)
      den <- f(xt, sigma) * dchisq(y, df = xt)
      if (u[i] <= num/den) x[i] <- y else {
      	x[i] <- xt
            k <- k+1     #y is rejected
      }
}
print(k)

index <- 5000:5500
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")


# qq plot과 히스토그램

b <- 2001      #discard the burnin sample
y <- x[b:m]
a <- ppoints(100)
QR <- sigma * sqrt(-2 * log(1 - a))  #quantiles of Rayleigh
Q <- quantile(y, a)
    
qqplot(QR, Q, main="", cex=.5,
	xlab="Rayleigh Quantiles", ylab="Sample Quantiles")
abline(0, 1)
hist(y, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR, f(QR, 4))
    

## 이변량 정규분포에 대한 깁스 표집기
N <- 5000               #length of chain
burn <- 1000            #burn-in length
X <- matrix(0, N, 2)    #the chain, a bivariate sample

rho <- -.75             #correlation
mu1 <- 0
mu2 <- 2
sigma1 <- 1
sigma2 <- .5
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

X[1, ] <- c(mu1, mu2)            #initialize
for (i in 2:N) {
	x2 <- X[i-1, 2]
      m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
      X[i, 1] <- rnorm(1, m1, s1)
      x1 <- X[i, 1]
      m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
      X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ]

# compare sample statistics to parameters
colMeans(x)
cov(x)
cor(x)

plot(x, main="", cex=.5, xlab=bquote(X[1]),
	ylab=bquote(X[2]), ylim=range(x[,2]))
         


## 수렴진단

# 데이터 생성
if (u[i] <= (dt(y,n)/dt(x[i-1],n))  {x[i] <- y} else
	{x[i] <- x[i-1]}

rw.Metropolis <- function(n,sigma,x0,N) {
	x <- numeric(N)
      x[1] <- x0
      u <- runif(N)
      k <- 0
      for (i in 2:N) {
		y <- rnorm(1, x[i-1], sigma)
            if (u[i] <= (dt(y,n)/dt(x[i-1],n))) x[i] <- y  else
            {
			x[i] <- x[i-1]
                  k <- k+1
            }
	}
	return(list(x=x, k=k))
}

n <- 4    # degrees of freedom for target Student t dist.
N <- 2000
sigma <- c(.05,.5,2,16)

x0 <- 25
rw1 <- rw.Metropolis(n, sigma[1], x0, N)
rw2 <- rw.Metropolis(n, sigma[2], x0, N)
rw3 <- rw.Metropolis(n, sigma[3], x0, N)
rw4 <- rw.Metropolis(n, sigma[4], x0, N)

# 몬테칼로 오차
MCE <- function(batch,x){
	v <- N/batch  # N : sample size after burn-in period
	batch.m <- numeric(batch)
      for(b in 1:batch){
      	batch.m[b] <- mean(x[((b-1)*v+1):(b*v)])
 	}
	sample.m <- mean(x)
      return(MCE <- sqrt(sum((batch.m-sample.m)^2)/(batch*(batch-1))))
}

N <- 1500  # burn-in : 500
b.rw1 <- rw1$x[501:2000];b.rw2<-rw2$x[501:2000]
b.rw3 <- rw3$x[501:2000];b.rw4<-rw4$x[501:2000]
cbind(MCE(30,b.rw1),MCE(30,b.rw2),MCE(30,b.rw3),MCE(30,b.rw4))

# 그래프를 이용하는 방법
acf.plot <- function(x){
	index <- c(1:N)
	acf(x, type="correlation", plot=TRUE, main="Autocorrelation")
}

par(mfrow=c(2,2))
acf.plot(b.rw1); acf.plot(b.rw2); acf.plot(b.rw3); acf.plot(b.rw4)

N <- 15500
rw3 <- rw.Metropolis(n,sigma[3],x0,N)
index <- seq(501,15500,by=10)  # sampling lag: 10
b.rw3 <- rw3$x[index]   # sample size: 1,500
acf.plot(b.rw3)