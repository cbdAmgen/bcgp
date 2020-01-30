create_x_and_x_test <- function(d, n, n_test, grid_test, randomX1D,
                                grid_test_size){

  if(d == 1){
    if(isTRUE(randomX1D)){
      x <- matrix(sort(runif(n, 0, 1)), ncol = 1)
    }else{
      x <- matrix(seq(0, 1, length.out = n), ncol = 1)
    }

    x_test <- matrix(seq(0, 1, length.out = n_test), ncol = 1)

  }else{
    x <- create_xmat(n, d)
    if(isTRUE(grid_test)){
      one_side <- seq(0, 1, length.out = grid_test_size)
      x_test <- as.matrix(expand.grid(rep(list(one_side), d)))
    }else{
      x_test <- create_xmat(n_test, d)
    }

  }
  return(list(x = x, x_test = x_test))
}

create_xmat <- function(n, d){
  if(requireNamespace("lhs", quietly = TRUE) && n <= 50 && d <= 5){
    x <- matrix(scaleX(lhs::optimumLHS(n, d)), ncol = d)
  }else{
    x <- matrix(scaleX(matrix(runif(n * d), nrow = n, ncol = d)), ncol = d)
  }
  return(x)
}
