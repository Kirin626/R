generate_samples <- function(n = 1e+04, p = 0.3, q = 0.3){
  r <- 1 - p - q
  phenotype <- c("O", "A", "B", "AB")
  prob <- c(r^2, p^2+2*p*r, q^2+2*q*r, 2*p*q)
  return(table(sample(phenotype, n, replace = T, prob = prob)))
}
