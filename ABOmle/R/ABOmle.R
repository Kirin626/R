ABOmle <- function(nO, nA, nB, nAB, p = 0.3, q = 0.3, acu = 1e-3, max_iter = 10, method = "N", special = FALSE){
  if(special){
    n <- nO + nA + nB +nAB
    r <- (nO/n)^(1/2); p <- (r^2+nA/n)^(1/2) - r; q <- 1 - r -p
  }
  # Input the nllf and calculate its Jacobi
  nllf <- expression((-2)*nO*log(1-p-q)-nA*log(p)-nA*log(2-p-2*q)
                     -nB*log(q)-nB*log(2-2*p-q)-nAB*(log(p)+log(q)))
  nllf.jacobi <- matrix(data = c(D(nllf, "p"), D(nllf, "q")), nrow = 2)
  nllf.jacobi <- as.expression(nllf.jacobi)
  if(method == "N"){
    # calculate its Hessian.
    nllf.hessian <- matrix(data = c(D(D(nllf, "p"),"p"),
                                    D(D(nllf, "p"),"q"),
                                    D(D(nllf, "q"),"p"),
                                    D(D(nllf, "q"),"q")),
                           nrow = 2)
    nllf.hessian <- as.expression(nllf.hessian)
    # Newton-Raphson Method
    p0 <- p + 1; q0 <- q + 1; iter <- 0
    while(abs(p-p0) + abs(q-q0) > acu & iter < max_iter){
      p0 <- p; q0 <- q; iter <- iter + 1
      temp.jacobi <- matrix(
        data = c(eval(nllf.jacobi[1, 1]), eval(nllf.jacobi[2, 1])),
        nrow = 2)
      temp.hessian <- matrix(
        data = c(
          eval(nllf.hessian[1, 1]),eval(nllf.hessian[2, 1]),
          eval(nllf.hessian[1, 2]),eval(nllf.hessian[2, 2])),
        nrow = 2)
      temp <- c(p, q) - solve(temp.hessian, temp.jacobi)
      p <- temp[1];q <- temp[2]
      if(!(p > 0 & p < 1 & q > 0 & q < 1)) {
        p <- runif(1, 0, 1)
        q <- runif(1, 0, 1 - p)
      }
    }
    if(abs(p-p0) + abs(q-q0) > 1e-3)
      cat("Cannot converge! Try a lager 'max_iter'.\n")
    else list(r = 1 - p - q, p = p, q = q, iter = iter)
  }
  else if(method == "F"){
    ## Input Fisher Information
    nllf.fisher <- matrix(
      data = c(
        expression(n*(2/p+p/(2-p-2*q)+4*q/(2-2*p-q))),
        expression(n*(2+2*p/(2-p-2*q)+2*q/(2-2*p-q))),
        expression(n*(2+2*p/(2-p-2*q)+2*q/(2-2*p-q))),
        expression(n*(2/q+q/(2-q-2*p)+4*p/(2-2*q-p)))
      ),
      nrow = 2
    )
    # Fisher Scoring
    n <- nO + nA + nB +nAB
    p0 <- p +1; q0 <- q +1; iter <- 0
    while(abs(p-p0) + abs(q-q0) > acu & iter < max_iter){
      p0 <- p; q0 <- q; iter <- iter + 1
      temp.jacobi <- matrix(
        data = c(eval(nllf.jacobi[1, 1]), eval(nllf.jacobi[2, 1])),
        nrow = 2)
      temp.fisher <- matrix(
        data = c(
          eval(nllf.fisher[1, 1]),eval(nllf.fisher[2, 1]),
          eval(nllf.fisher[1, 2]),eval(nllf.fisher[2, 2])),
        nrow = 2)
      temp <- c(p, q) - solve(temp.fisher, temp.jacobi)
      p <- temp[1];q <- temp[2]
      if(!(p > 0 & p < 1 & q > 0 & q < 1)) {
        p <- runif(1, 0, 1)
        q <- runif(1, 0, 1 - p)
      }
    }
    if(abs(p-p0) + abs(q-q0) > 1e-3 )
      cat("Cannot converge! Try a lager 'max_iter'.")
    else  list(r = 1 - p - q, p = p, q = q, iter = iter)
  }
}






