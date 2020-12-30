# pseudo random number
runif(5)
runif(10, min=-3, max=-1)
set.seed(32789)
runif(5)
set.seed(32789)
runif(5)

# sample function
x = sample(1:3, size=100, replace=T, prob=c(.2, .3, .5))
table(x)

# Buffon's needle
Buffon = function(n, lofneedle, distance)
{
	lofneedle = lofneedle / 2
	distance = distance / 2
	r1 = runif(n)
	r2 = runif(n)
	prob = mean(r1*distance < lofneedle*sin(r2*pi))
	return(prob)
}

Buffon(5000,15,20)

