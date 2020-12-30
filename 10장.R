### EM algorithm

## Example 1: multinomial 
eps = 1e-5
theta = diff = 0.5
k = 0
result = c(k, theta, diff)
while(diff > eps) {
	y = 125 * theta / (theta+2)
	th_hat = (34+y) / (38+34+y)
	diff = abs(theta - th_hat)
	theta = th_hat
	k = k + 1
	result = rbind(result, c(k, theta, diff))
}
round(result, 8)

## Example 2: censored data(심장병환자들의 생존시간)
y1 = c(39, 2, 101, 20, 31, 15) # 절단 x
y0 = c(118, 91, 427)           # 절단 o
n = length(c(y1,y0))
r = length(y1)

eps = 1e-5
maxiter = 1000
mu.old = mean(y1) # 초기값

for (i in 1:maxiter) {
	mu.new = (sum(y1) + sum(y0) + (n-r)*mu.old) / n
      diff = abs(mu.old - mu.new) / abs(mu.old)
	if (diff < eps) break
	mu.old = mu.new
}
cat("iter = ", i, "estimate= ", mu.new,"\n")




## Example 3: normal mixture

Log.lik = function(x, R=2, mu, sigma, prior)
{   
   lik = 0
   for (r in 1:R)
      lik = lik + prior[r] * dnorm(x, mean = mu[r], sd = sigma[r]) 
   return(sum(log(lik)))
}

Normal.Mixture = function(X, R=2, maxiter=100, eps=1e-5)
{
   X = as.vector(X)
   N = length(X)
   mu = sigma = prior = rep(0, R)
   gama = matrix(0, R, N)
   # find initial centroids using K-means clustering
   prior = rep(1/R, R)
   kmfit = kmeans(X, R)
   mu = kmfit$centers
   sigma = sqrt(kmfit$withinss /(kmfit$size - 1))     
   old.lik = Log.lik(X, R, mu, sigma, prior) 
   track.lik = as.vector(NULL)
   track.lik = c(old.lik)
   for (i in 1:maxiter)
   {
      for (r in 1:R)
             gama[r, ] = prior[r] * dnorm(X, mean = mu[r], sd = sigma[r])
      denom = apply(gama, 2, sum)
      for (r in 1:R)
      {
         gama[r, ] = gama[r, ] / denom
         mu[r] = t(gama[r, ]) %*% X / sum(gama[r, ])
         sigma[r] = sqrt(t(gama[r, ]) %*% (X - mu[r])^2 / sum(gama[r, ]))
      }
      prior = apply(gama, 1, sum) / N
      new.lik = Log.lik(X, R, mu, sigma, prior)
      if (abs(old.lik - new.lik) < eps * abs(old.lik))  break
      old.lik = new.lik
      track.lik = c(track.lik, old.lik)
   }
   return(list(mu = mu, sigma = sigma, prior = prior, 
                  track = track.lik, resp = gama[r, ]))
}

Mixture.prob = function(x, mu, sigma, prior)
{
   R = length(mu)
   prob = 0
   for (r in 1:R)
      prob = prob + prior[r] * dnorm(x, mu[r], sigma[r]) 
   return(prob)
}

X = c(-0.39, 0.12, 0.94, 1.67, 1.76, 2.44, 3.72, 4.28, 4.92, 5.53,
         0.06, 0.48, 1.01, 1.68, 1.80, 3.25, 4.12, 4.60, 5.28, 6.22)
fit = Normal.Mixture(X)     

hist(X, nclass = 20, xlab = "X", freq = FALSE, col = 2, ylim = c(0, 1))

X.grid = seq(min(X), max(X), length=100)
plot(X.grid, Mixture.prob(X.grid, fit$mu, fit$sigma, fit$prior), 
   type = "l", ylab = "density", xlab = "X",
   ylim = c(0, 1), col = 2)
lines(sort(X), 1-fit$resp[sort(X, index.return = T, decreasing=T)$ix], 
   type = "b", col = 3)

plot(1:length(fit$track), fit$track, 
   xlab = "iteration", ylab = "Obs. log-likelihood", type = "b", 
   col = 3)