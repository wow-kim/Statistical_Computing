## 적중법
# 예 7-1
myPi.1 = function(n)
{
	x = matrix(runif(n*2), ncol=2)
	r = mean(apply(x^2, 1, sum) < 1)
	return(4*r)
}

myPi.1(4000000)


## 표본평균 몬테칼로
# pi
Sample.Mean = function(n) 
{
	x = runif(n)
	return(4*mean(sqrt(1-x^2)))
}

Sample.Mean(4000000)


## 분산감소: 대조변수법
# Phi(x) = E[xe^{-(xU)^2/2}], U ~ U[0,1]
MC.Phi <- function(x, R=10000, antithetic=TRUE) {
	u <- runif(R/2)
 	if (!antithetic) v <- runif(R/2) else
      	v <- 1-u
     	u <- c(u, v)
	cdf <- numeric(length(x))
	for (i in 1:length(x)) {
      	g <- x[i]*exp(-(u*x[i])^2/2)            
            cdf[i] <- 0.5+mean(g)/sqrt(2*pi)        
    	}
     	cdf
}

x <- seq(.1, 2.5, length=5)
Phi <- pnorm(x)

set.seed(123)
MC1 <- MC.Phi(x, anti = FALSE)

set.seed(123)
MC2 <- MC.Phi(x)
# 추정값 비교
print(round(rbind(x, MC1, MC2, Phi), 5))

# 분산 비교
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.95
for (i in 1:m) {
	MC1[i] <- MC.Phi(x, R=1000, anti=FALSE)
	MC2[i] <- MC.Phi(x, R=1000)
}
print(sd(MC1))
print(sd(MC2))
print(100*(var(MC1)-var(MC2))/var(MC1))      


## 분산감소 기법: 주표본 기법
# int_0^1 exp(-x)/(1+x^2) dx

m = 10000
theta.hat = se = numeric(3)
g = function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

# candidate 1: I(0<x<1)
x = runif(m)     
fg = g(x)
theta.hat[1] = mean(fg)
se[1] = sd(fg)

# candiate 2: exp(-x)
x = rexp(m, 1) 
fg = g(x) / exp(-x)
theta.hat[2] = mean(fg)
se[2] = sd(fg)

# candiate 3: 4*(pi *(1+x^2))^{-1} I(0<x<1) (Cauchy dist on (0,1))
u = runif(m)    # inverse transform method
x = tan(pi * u / 4)
fg = g(x) / (4 / ((1 + x^2) * pi))
theta.hat[3] = mean(fg)
se[3] = sd(fg)

rbind(theta.hat, se)

## 분산감소 기법: 제어변량
# E[e^U]
m <- 10000
a <- -12+6*(exp(1)-1)    # a^*
U <- runif(m)

T1 <- exp(U)                        # simple MC
T2 <- exp(U)+a*(U-1/2)              # controlled

mean(T1)    
mean(T2)    

100*(var(T1)-var(T2))/var(T1)

# a^* 추정하기
# int_0^1 e^{-x}/(1+x^2) dx
# g(x) = e^{-x}/(1+x^2)과 가까운 f(x)=e^{-.5}/(1+x^2)

f <- function(u)  exp(-.5)/(1+u^2)
g <- function(u)  exp(-u)/(1+u^2)
x <- seq(0,1,by=0.01)
plot(x, f(x), type="l", ylim=c(0,1), ylab="")
lines(x, g(x), ylim=c(0,1), lty=2)
legend("topright", 1, c("g(x)", "f(x)"), lty = c(1:2), inset = .02)


# 파일럿 스터디에 의한 a^* 추정
set.seed(510)    # needed later
u <- runif(10000)

B <- f(u)
A <- g(u)

cor(A,B)
a <- -cov(A,B)/var(B)      
a

# 본 모의실험
m <- 100000
u <- runif(m)

T1 <- g(u)
T2 <- T1+a*(f(u)-exp(-.5)*pi/4)
c(mean(T1), mean(T2))
c(var(T1), var(T2))
100*(var(T1)-var(T2))/var(T1)
