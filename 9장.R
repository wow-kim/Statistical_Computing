## boot 함수
library(boot)
library(bootstrap)

## 표준오차에 대한 붓스트랩 추정
data(law)
data(law82)
cor(law$LSAT, law$GPA)
cor(law82$LSAT, law82$GPA)

B = 200
n = nrow(law)
R = numeric(B)

# 직접 붓스트랩 구현
for (b in 1:B) {
	i = sample(1:n, size=n, replace=T)
	LSAT = law$LSAT[i]
	GPA = law$GPA[i]
	R[b] = cor(LSAT, GPA)
}
print(se.R <- sd(R))
hist(R, prob=TRUE)

# boot 함수 이용
r = function(x, i) {
	cor(x[i,1], x[i,2])
}

obj = boot(data=law, statistic=r, R=2000) # R=B붓스트랩 반복수
obj

## 편의에 대한 붓스트랩 추정
theta.hat = cor(law$LSAT, law$GPA)
B = 2000
n = nrow(law)
theta.b = numeric(B)

for (b in 1:B) {
	i = sample(1:n, size=n, replace=T)
	LSAT = law$LSAT[i]
	GPA = law$GPA[i]
	theta.b[b] = cor(LSAT, GPA)
}
bias = mean(theta.b) - theta.hat
bias


## interval estimation for sd
dice = c(1,2,3,2,6,6,5,1,1,1,4,2,4,1,4,5,6,6,3,2,
	 5,6,4,1,2,3,2,2,5,3,5,6,1,4,4,4,3,5,5,1,
	 6,1,3,3,2,5,2,2,1,4)

sd(dice)
boot.sample = matrix(0, 200, 50)
std = rep(0, 200)
set.seed(12345)

for(n in 1:200) {
	boot.sample[n,] = sample(dice, replace=T, 50)
	std[n] = sd(boot.sample[n,])
}

summary(std)
hist(std, xlim=c(1.2,2.2), ylim=c(0,50), xlab="bootst.sd")
quantile(std, probs=c(0.025, 0.975))



## jackknife for patch data
data(patch, package="bootstrap")
n = nrow(patch)
y = patch$y
z = patch$z
theta.hat = mean(y)/mean(z)
theta.hat

# 편의
theta.jack = numeric(n)
for (i in 1:n)
	theta.jack[i] = mean(y[-i])/mean(z[-i])
bias = (n-1) * (mean(theta.jack) - theta.hat)
bias

# 표준오차
se = sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
print(se)

