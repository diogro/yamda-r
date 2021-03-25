
# https://stackoverflow.com/questions/12982528/how-to-create-an-r-function-programmatically
make_function <- function(args, body, env = parent.frame()) {
  f <- function() {}
  formals(f) <- args
  body(f) <- body
  environment(f) <- env
  f
}
make_alist <- function(args) {
  res <- replicate(length(args), substitute())
  names(res) <- args
  res
}