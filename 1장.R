# f(x) = (1-cos x)/x^2: 0근방에서 부동소수점연산에 의해 오차가 생김

evaluatefunction1 = function(xmin,xmax,n){
	x = c(0)
	f = c(0)
	for (i in (0:n)){
		x[i+1] = xmin + i*(xmax-xmin)/n
		f[i+1] = (1-cos(x[i+1]))/(x[i+1])^2
	}
	plot(x,f,type="l",col="blue",xlab="x",ylab="function")
}

evaluatefunction1(-1,1,100)

evaluatefunction1(-10^-7, 10^-7, 200)


# 동치인 함수 f(x) = (sin^2(x/2))/(x^2/2)를 고려하면 0근처에서 
# 오차의 영향을 줄일 수 있지만 0에서 불연속
evaluatefunction2 = function(xmin,xmax,n){
	x = c(0)
	f = c(0)
	for (i in (0:n)){
		x[i+1] = xmin + i*(xmax-xmin)/n
		f[i+1] = (sin(x[i+1]/2))^2/((x[i+1])^2/2)
	}
	plot(x,f,type="l",col="blue",xlab="x",ylab="function")
}

evaluatefunction2(-10^-7, 10^-7, 200)


# f(x) = (sin^2(x/2))/(x^2/2)을 0에서 연속이 되도록 정의
evaluatefunction2withcheck = function(xmin,xmax,n,epsilon){
	x = c(0)
	f = c(0)
	for (i in (0:n)){
		x[i+1] = xmin + i*(xmax-xmin)/n
		if (abs(x[i+1]) > epsilon){
			f[i+1] = (sin(x[i+1]/2))^2/((x[i+1])^2/2)
		}
		else{
			f[i+1] = 1/2
		}
	}
	plot(x,f,type="l",col="blue",xlab="x",ylab="function")
}

evaluatefunction2withcheck(-10^-7, 10^-7, 200, 10^-10)


## n차 다항식의 계산
# 직접 계산
polyeval = function(a,x0){
	n = length(a) - 1
	polyvalue = a[1]
	for (j in (1:n)){
		polyvalue = polyvalue + a[j+1]*(x0)^j
	}
	return(polyvalue)
}

# 멱급수 계산을 하지 않도록 수정
polyevalimproved = function(a,x0){
	n = length(a) - 1
	polyvalue = a[1]
	powersofx0 = 1
	for (j in (1:n)){
		powersofx0 = x0*powersofx0
		polyvalue = polyvalue + a[j+1]*powersofx0
	}
	return(polyvalue)
}

# Horner의 알고리즘
polyevalhorner = function(a,x0){
	n = length(a) - 1
	polyeval = a[n+1]
	for (j in (1:n)){
		polyeval = a[n+1-j] + x0*polyeval
	}
	return(polyeval)
}

# p(x) = (x-2)^{10}를 $x=2.01$에서 계산하는 경우 계산결과가 
# rounding error에 의해 달라짐
a = c(1024, -5120, 11520, -15360, 13440, -8064, 3360, -960, 180, -20, 1)
polyeval(a, 2.01)
polyevalimproved(a, 2.01)
polyevalhorner(a, 2.01)

# p(x) = (x-1)^3 = x^3 - 3x^2 + 3x - 1를 x=1에서 계산
polyevalhornermultiple = function(a,xmin,xmax,n){
	x = c(0)
	y = c(0)
	for (i in (0:n)){
		x[i+1] = xmin + i*(xmax-xmin)/n
		y[i+1] = polyevalhorner(a,x[i+1])
	}
	plot(x,y,type="l",col="blue",xlab="x",ylab="function",cex.axis=1.5,cex.lab=1.5)
}

polyevalhornermultiple(-c(1,-3,3,-1), 0.99999, 1.00001, 200)
curve((x-1)^3, from=0.99999, to=1.00001)

## 행렬계산
# epsilon이 작은 경우
epsilon = 10^-10
A = matrix(c(1, epsilon, 0, 1, 0, epsilon), nrow=3)
A
B = t(A) %*% A
B

epsilon = 10^-10
v = c(1, epsilon, -1)
w = c(1, epsilon, 1)
t(v) %*% w


## 실수의 비교시 주의

isTRUE(all.equal(.2, .3 - .1))
all.equal(.2, .3)         
isTRUE(all.equal(.2, .3)) 
