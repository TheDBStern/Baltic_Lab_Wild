files <- dir('.', pattern = 'omega.out')

X = lapply(
  files,
  function(x) 
  {
  	print(x)
  	table <- as.matrix(read.table(x))
    table
  }
)
 
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

res <- apply(Y, c(1, 2), mean, na.rm = TRUE)