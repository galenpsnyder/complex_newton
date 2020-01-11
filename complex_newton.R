# test functions
rosenbrock <- function(x){
  y <- x[1]
  x <- x[2]
  (1 - x)^2 + 100 * (y - x^2)^2
}

f <- function(x) (x[1]^2 * x[2] * x[3] * x[4]^2) + (x[2]^2 * x[3]^3 * x[4])


# this function estimates the first derivative of a function using the complex-step method
fpc <- function(x, h, f, ...){
  k <- length(x)                      # number of parameters specified
  jac <- numeric(length = k)          # initializing jacobian vetor
  pmat <- diag(1, k)                  # identity matrix used to perturb single parameters to estimate partial derivatives
  for(j in seq_len(k)){               # loop over all parameters for partial derivatives
    im <- h*pmat[, j]                 # perturbing single parameter -- becomes imaginary part of complex number
    comp.x <- complex(length.out = k, real = x, imaginary = im) # creating complex number vector
    jac[j] <- Im(f(comp.x, ...)) / h  # placing first derivative estimates in jacobian
  }
  jac                                 # return vector of first derivatives
}

# testing function
fpc(x = c(1, -1.2), h = 1e-8, f = rosenbrock)
fpc(x = c(5, 3, 6, 4), h = 1e-8, f = f)


# this function estimates the second derivative using an extension of the complex-step method
fpc_sec <- function(x, h, f, ...){
  k <- length(x)                # number of parameters specified
  hes <- matrix(0, k, k)        # initializing empty hessian matrix
  pmat2 <- pmat1 <- diag(1, k)  # creating perturbation matrices
  for(m in 1:k){                # loop over first derivatives
    for(n in 1:k){              # loop over parameters
      comp.x <- complex(length.out = k, real = x, imaginary = h*pmat1[, n]) # create complex vector
      hes[n, m] <- Im(f(comp.x + h*pmat2[, m], ...) - f(comp.x - h*pmat2[, m], ...)) # second derivative calculation
    }
  }
  hes <- hes / (2*h*h)          
  hes                           # scale and return hessian
}

# testing
fpc_sec(x = c(1, 1), h = 1e-8, f = rosenbrock)
fpc_sec(x = c(5, 3, 6, 4), h = 1e-8, f = f)

### all methods decribed for complex-step differentiation adapted from:
### Ridout, M. S. (2009). Statistical Applications of the Complex-Step Method of 
### Numerical Differentiation. The American Statistician, 63(1), 66-74.

# this function uses the derivative functions created above in a quasi-newton optimization module
complex_newton <- function(objective, parm, hessian = F, tol = .Machine$double.eps, ...){
  k <- 1
  diff <- 1
  obj_old <- objective(parm, ...)
  grad <- fpc(x = parm, h = 1e-8, f = objective, ...)
  
  while(abs(diff) > tol){
    hess <- fpc_sec(x = parm, h = 1e-8, f = objective, ...)
    p    <- solve(hess) %*% -grad
    parm <- parm + p
    grad <- fpc(x = parm, h = 1e-8, f = objective, ...)
    obj  <- objective(parm, ...)
    diff <- obj - obj_old
    obj_old <- obj
    k <- k + 1
  }
  
  if(hessian){
    out <- list(parm = parm,
                value = objective(parm, ...),
                ierations = k,
                hessian = hess)
  } else {
    out <- list(parm = parm,
                value = objective(parm, ...),
                ierations = k) 
  }
  
  out
}


######### updated module -- estimates first as second derivatives simultaneously ##########
fpc <- function(var, h, f, ...){
  nvar <- length(var)
  jac  <- numeric(length = nvar)
  hes  <- matrix(0, nvar, nvar)
  pmat <- diag(1, nvar)
  for(m in 1:nvar){
    pmat.i <- h*pmat[, m]
    for(n in 1:m){
      comp.var <- complex(length.out = nvar, real = var, imaginary = h*pmat[, n])
      eval.a   <- f(comp.var + pmat.i, ...)
      eval.s   <- f(comp.var - pmat.i, ...)
      if(n == m){
        jac[m] <- Im(eval.a + eval.s)
      }
      hes[n, m] <- hes[m, n] <- Im(eval.a - eval.s)
    }
  }
  jac <- jac / (2*h)
  hes <- hes / (2*h*h)
  out <- list(value = Re(eval.a),
              first_deriv = jac,
              second_deriv = hes)
  out
}


## updated optimizer--reflects changes in derivative modules ##
complex_newton <- function(objective, parm, hessian = F, tol = sqrt(.Machine$double.eps), max.iter = Inf, ...){
  iter <- 1
  grad <- Inf
  
  while(sqrt(sum(grad^2)) > tol & iter < max.iter){
    eval <- fpc(var = parm, h = 1e-10, f = objective, ...)
    grad <- eval$first_deriv
    hess <- eval$second_deriv
    qr.h <- qr(hess)
    p    <- qr.solve(qr.h, -grad)
    parm <- parm + p
    iter <- iter + 1
  }
  
  if(hessian){
    out <- list(parm = parm,
                value = eval$value,
                ierations = iter,
                hessian = hess)
  } else {
    out <- list(parm = parm,
                value = eval$value,
                ierations = iter) 
  }
  
  out
}

# testing on rosebrock function -- gives good results
complex_newton(objective = rosenbrock, parm = c(1, -1.2))

## testing complex newton optimizer on logistic regression function
# Linear combination of predictors
v<-function(X, beta){
  v<-X%*%beta
  v
}

# Estimated Probability through sigmoid transformation
p_hat<-function(v){
  p<-1/(1+exp(-v))
  p
}

# Log-likelihood function
ll<-function(y, p){
  -sum(y*log(p)+(1-y)*log(1-p))
}


# logistic regression optimizer module
log_reg <- function(beta, X, y){
  v_i<-v(X = X, beta = beta)
  p<-p_hat(v = v_i)
  ll(y = y, p = p)
}

X<-log_sim[3:7]
y<-log_sim["decision"]

X<-cbind(1, X)
X<-as.matrix(X)
y<-as.matrix(y)

out <- complex_newton(objective = log_reg, parm = rep(0, 6), hessian = T, X = X, y = y)









