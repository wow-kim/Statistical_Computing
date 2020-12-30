# constructing matrix
H3 = 1/cbind(1:3, 2:4, 3:5)
matrix(seq(1,12), nrow = 3)
x = seq(1,3)
x2 = x^2
X = cbind(x, x2)
X
X = matrix(c(1,2,3,1,4,9), ncol = 2)
X

# accessing matrix elements
X[3,2]
X[3,]
X[,2]
colnames(X)
rownames(X)
rownames(X) = c("obs1","obs2","obs3")
X

# matrix properties
mat = rep(1:4, rep(3,4))
dim(mat) = c(3,4)
mat
dim(mat)

# diagonal matrix
D = diag(c(1,2,3))
D
diag(D)
I = diag(3)
I

# triangluar matrix
lower.tri(H3)
Hnew = H3
Hnew[upper.tri(H3)] = 0
Hnew


# transpose
A = matrix(c(rep(1,5), -1, 1, -1), byrow = T, ncol = 4)
A
A1 = t(A)
A1


# matrix arithmetic
Y = 2 * X
Y
Y + X
t(Y) + X
X * Y

# matrix multiplication
t(Y) %*% X
Y %*% X
crossprod(Y,X)


# determinant and inverse
A = matrix(10:13, 2)
B = matrix(c(3,5,8,9), 2)
det(A)
det(t(A))
det(A %*% B) 
det(A) * det(B)


H3inv = solve(H3)
H3inv
H3inv %*% H3 
zapsmall(H3inv %*% H3)

# rank
library(fBasics)
set.seed(45)
A = sample(11:26)
A = matrix(A, 4)
rk(A)
B = A
B[1,] = rep(0, ncol(B))
rk(B)

# generalized inverse
library(limSolve)
A = matrix(c(1:8, 6, 8, 10, 12), nrow=4, ncol=3)
B = 0:3
X = Solve(A, B)
A %*% X - B
(gA = Solve(A))



# solving linear equations
A = matrix(c(4,2,1,-2,-1,3,3,1,-1), ncol = 3)
A
b = c(9,3,4)
solve(A,b)



# Singular value decomposition
H3.svd = svd(H3)
H3.svd
H3.svd$u %*% diag(H3.svd$d) %*% t(H3.svd$v)
H3.svd$v %*% diag(1/H3.svd$d) %*% t(H3.svd$u) # inverse

# Choleski decomposition
H3.chol = chol(H3)
H3.chol
crossprod(H3.chol, H3.chol)
chol2inv(H3.chol)	# inverse
# H_3 x = b, b = (1,2,3)'
b = 1:3
y = forwardsolve(t(H3.chol), b) 	# forward elimination (lower)
backsolve(H3.chol,y)			# back substitution (upper)

# QR decomposition
H3.qr = qr(H3)
H3.qr
Q = qr.Q(H3.qr)
Q
R = qr.R(H3.qr)
R
Q %*% R
qr.solve(R) %*% t(Q)	# inverse

# LU decomposition
set.seed(1)
library(Matrix)
mm = Matrix(round(rnorm(9),2), nrow=3)
mm

lum = lu(mm)
str(lum)

elu = expand(lum)
elu

with(elu, P %*% L %*% U)

# spectral decomposition
A = matrix(c(5,25,35,25,155,175,35,175,325), ncol = 3)
A
EA = eigen(A, symmetric = T)
EA
det(A)
prod(EA$values)


# rank: SV and eigenvalue
rk(A)
A.svd = svd(A)
sum(A.svd$d != 0)
sum(EA$values != 0)


# inner and outer product
x1 = 1:5
outer(x1, x1, "/")
outer(x1, x1, "-")
y = 5:10
outer(x1, y, "+")

# real symmetric and positive definite matrix
A1 = matrix(c(4, 345, 193, 297), 2)
eigen(A1)$value
library(matrixcalc)
is.positive.definite(A1)
A2 = (A1 + t(A1)) / 2
eigen(A2)$value
is.positive.definite(A2)


## 해법
# 전방 대입법 
forwardsub = function(L,b){
	x = c(0)
	n = nrow(L)
	for (i in (1:n)){
		x[i] = b[i]
		if (i > 1){
			for (j in (1:(i-1))){
				x[i] = x[i] - L[i,j]*x[j]
			}
		}
		x[i] = x[i]/L[i,i]
	}
	return(cbind(x))
}

# 후방 대입법
backwardsub = function(U,b){
	x = c(0)
	n = nrow(U)
	for (i in (n:1)){
		x[i] = b[i]
		if (i < n){
			for (j in ((i+1):n)){
				x[i] = x[i] - U[i,j]*x[j]
			}
		}
		x[i] = x[i]/U[i,i]
	}
	return(cbind(x))
}

# 가우스 소거법

gaussianelimination = function(Ab){
	n = nrow(Ab)
	for (k in (1:(n-1))){
		for (i in ((k+1):n)){
			mik = Ab[i,k]/Ab[k,k]
			Ab[i,k]=0
			for (j in ((k+1):(n+1))){
				Ab[i,j] = Ab[i,j] - mik*Ab[k,j]
			}
		}
	}
	return(Ab)
}

# 가우스 소거법(부분 피버팅)
gaussianeliminationpartial = function(Ab){
	n = nrow(Ab)
	for (k in (1:(n-1))){
		pivotindex = k
		for (i in ((k+1):n)){
			if (abs(Ab[i,k]) > abs(Ab[pivotindex,k])){
				pivotindex = i
			}
		}
		if (pivotindex != k){
			for (j in (k:(n+1))){
				buffer = Ab[k,j]
				Ab[k,j] = Ab[pivotindex,j]
				Ab[pivotindex,j] = buffer
			}
		}
		for (i in ((k+1):n)){
			mik = Ab[i,k]/Ab[k,k]
			Ab[i,k] = 0
			for (j in ((k+1):(n+1))){
				Ab[i,j] = Ab[i,j] - mik*Ab[k,j]
			}
		}
	}
	return(Ab)
}


# LU 분해 
lufactorization = function(A){
	n = nrow(A)
	L = matrix(0,nrow=n,ncol=n)
	for (k in (1:(n-1))){
		for (i in ((k+1):n)){
			L[i,k] = A[i,k]/A[k,k]
			A[i,k] = 0
			for (j in ((k+1):n)){
				A[i,j] = A[i,j] - L[i,k]*A[k,j]
			}
		}
	}
	for (k in (1:n)) L[k,k] = 1
	return(cbind(L,A))
}

# 촐레스키 분해
choleskyfactorization = function(A){
	n = nrow(A)
	L = matrix(0,nrow=n,ncol=n)
	for (i in (1:n)){
		L[i,i] = A[i,i]
		if (i > 1){
			for (k in (1:(i-1))){
				L[i,i] = L[i,i] - L[i,k]*L[i,k]
			}
		}
		L[i,i] = (L[i,i])^(1/2)
		if (i < n){
			for (j in ((i+1):n)){
				L[j,i] = A[j,i]
				if (i > 1){
					for (k in (1:(i-1))){
						L[j,i] = L[j,i] - L[j,k]*L[i,k]
					}
				}
				L[j,i] = L[j,i]/L[i,i]
			}
		}
	}
	return(L)
}

# 가우스 조던 소거법
gaussjordanelimination = function(Ab){
	n = nrow(Ab)
	for (k in (1:n)){
		if (k > 1){
			for (i in (1:(k-1))){
				mik = Ab[i,k]/Ab[k,k]
				Ab[i,k] = 0
				for (j in ((k+1):(n+1))){
					Ab[i,j] = Ab[i,j] - mik*Ab[k,j]
				}
			}
		}
		if (k < n){
			for (i in ((k+1):n)){
				mik = Ab[i,k]/Ab[k,k]
				Ab[i,k] = 0
				for (j in ((k+1):(n+1))){
					Ab[i,j] = Ab[i,j] - mik*Ab[k,j]
				}
			}			
		}
	}
	return(Ab)
}