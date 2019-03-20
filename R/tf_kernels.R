# tensorflow implementations of common kernels

# bias (or constant) kernel
# k(x, y) = \sigma^2
tf_bias <- function(X, X_prime, variance) {
  
    # calculate dims from X and X_prime
    dims_out <- tf$stack(c(tf$shape(X)[0L], tf$shape(X)[1L], tf$shape(X_prime)[1L]))
  
    # create and return covariance matrix
    tf$fill(dims_out, tf$squeeze(variance))

}

# squared_exponential kernel (RBF)
tf_rbf <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_slice(X, active_dims)
  X_prime <- tf_slice(X_prime, active_dims)

  # calculate squared distances
  r2 <- scaled_square_dist(X, X_prime, lengthscales)
  
  # construct and return RBF kernel
  variance * tf$exp(-r2 / tf$constant(2.0, dtype = tf$float64))
  
}

# rational_quadratic kernel
# tf_rational_quadratic <- function(X, Y = NULL, sigma, alpha, lengthscales = NULL) {
#   
#   # precalculate variance
#   sigma2 <- sigma * sigma
#   
#   # calculate squared distances (scaled if needed)
#   r2 <- scaled_square_dist(X, Y, lengthscales)
#   
#   # construct and return rational quadratic kernel
#   sigma2 * (1 + r2 / (2 * alpha)) ^ -alpha
#   
# }

# linear kernel (base class)
#' tf_linear <- function(X, Y = NULL, variance) {
#'   
#'   # full kernel
#'   def K(self, X, X2=None, presliced=False):
#'     if not presliced:
#'     X, X2 = self._slice(X, X2)
#'     if X2 is None:
#'       return tf.matmul(X * self.variance, X, transpose_b=True)
#'     else:
#'       return tf.tensordot(X * self.variance, X2, [[-1], [-1]])
#'     
#'     # diag only
#'     params_as_tensors
#'     def Kdiag(self, X, presliced=False):
#'       if not presliced:
#'       X, _ = self._slice(X, None)
#'     return tf.reduce_sum(tf.square(X) * self.variance, -1)
#'     
#' }

# polynomial kernel (linear class)
#' tf_polynomial <- function(X, Y = NULL, variance, offset, degree) {
#'   
#'   params_as_tensors
#'   def K(self, X, X2=None, presliced=False):
#'     return (Linear.K(self, X, X2, presliced=presliced) + self.offset) ** self.degree
#'   
#'   params_as_tensors
#'   def Kdiag(self, X, presliced=False):
#'     return (Linear.Kdiag(self, X, presliced=presliced) + self.offset) ** self.degree
#'   
#' }

# rescale, calculate, and return clipped Euclidean distance
# scaled_dist <- function(X, Y, lengthscales = NULL) {
# 
#   stop("kernels based on absolute Euclidean distance are not implemented (yet)", call. = FALSE)
#     
# }

# # exponential kernel (stationary class)
# tf_exponential <- function(X, Y = NULL, sigma, lengthscales = NULL) {
#   
#   # precalculate variance
#   sigma2 <- sigma * sigma
#   
#   # calculate squared distances (scaled if needed)
#   r <- scaled_dist(X, Y, lengthscales)
#   
#   # construct and return exponential kernel
#   sigma2 * tf$exp(-0.5 * r)
# 
# }

# # Matern12 kernel (stationary class)
# tf_Matern12 <- function(X, Y = NULL, sigma, lengthscales = NULL) {
#   
#   # precalculate variance
#   sigma2 <- sigma * sigma
#   
#   # calculate squared distances (scaled if needed)
#   r <- scaled_dist(X, Y, lengthscales)
#   
#   # construct and return Matern12 kernel
#   sigma2 * tf$exp(-r)
#   
# }
# 
# # Matern32 kernel (stationary class)
# tf_Matern32 <- function(X, Y = NULL, sigma, lengthscales = NULL) {
#   
#   # precalculate variance
#   sigma2 <- sigma * sigma
#   
#   # calculate squared distances (scaled if needed)
#   r <- scaled_dist(X, Y, lengthscales)
# 
#   # precalculate root3
#   sqrt3 <- sqrt(3.0)
#     
#   # construct and return Matern32 kernel
#   sigma2 * (1.0 + sqrt3 * r) * tf$exp(-sqrt3 * r)
# 
# }
# 
# # Matern52 kernel (stationary class)
# tf_Matern52 <- function(X, Y = NULL, sigma, lengthscales = NULL) {
#   
#   # precalculate variance
#   sigma2 <- sigma * sigma
#   
#   # calculate squared distances (scaled if needed)
#   r <- scaled_dist(X, Y, lengthscales)
#   
#   # precalculate root5
#   sqrt5 <- sqrt(5.0)
#   
#   # construct and return Matern52 kernel
#   sigma2 * (1.0 + sqrt5 * r + 5.0 / 3.0 * tf$square(r)) * tf$exp(-sqrt5 * r)
#   
# }
# 
# # cosine kernel (stationary class)
# tf_cosine <- function(X, Y = NULL, sigma, lengthscales = NULL) {
#   
#   # precalculate variance
#   sigma2 <- sigma * sigma
#   
#   # calculate squared distances (scaled if needed)
#   r <- scaled_dist(X, Y, lengthscales)
#   
#   # construct and return cosine kernel  
#   sigma2 * tf$cos(r)
# 
# }
# 
# # calculate weighted products used in arc_cosine kernel
# weighted_product <- function(X, Y = NULL, weight_variance, bias_variance) {
#   
#   if (is.null(Y)) {
#     out <- tf$reduce_sum(weight_variance * tf$square(X), axis = -1L) + bias_variance
#   } else {
#     out <- tf$matmul(weight_variance * X), Y, transpose_b = TRUE) + bias_variance
#   }
# 
#   out
#   
# }
# 
# # calculate J (family of functions in reference paper for arc_cosine kernel)
# calc_J <- function(theta, order) {
#   
#   if (order == 0L) {
#     out <- pi - theta
#   }
#   if (order == 1L) {
#     out <- tf$sin(theta) + (pi - theta) * tf$cos(theta)
#   }
#   if (order == 2L) {
#     out <- 3.0 * tf$sin(theta) * tf$cos(theta) + (pi - theta) * (1.0 + 2.0 * tf$cos(theta) ^ 2)
#   }
#   
# }
# 
# # arc_cosine kernel (base kernel class)
# tf_arc_cosine <- function(X, Y = NULL, sigma, order, weight, bias, lengthscales = NULL) {
#   
#   stop("arc_cosine kernel is not implemented", call. = FALSE)
#   
#   # check parameters
#   # if (!(order %in% c(0L, 1L, 2L)))
#   #   stop("order must be one of 0, 1, or 2", call. = FALSE)
#   
#   # precalculate variances
#   # sigma2 <- sigma * sigma
#   # weight_sigma2 <- weight_sigma * weight_sigma
#   # bias_sigma2 <- bias_sigma * bias_sigma
#   
#   # calculate full kernel
#   # X_denominator <- tf$sqrt(weighted_product(X, NULL, weight_sigma2, bias_sigma2))
#   # if (is.null(Y)) {
#   #   Y <- X
#   #   Y_denominator <- X_denominator
#   # } else {
#   #   Y_denominator <- tf$sqrt(weighted_product(Y, NULL, weight_sigma2, bias_sigma2))
#   # }
#   # numerator <- weighted_product(X, Y, weight_sigma2, bias_sigma2)
#   # X_denominator <- tf$expand_dims(X_denominator, -1L)
#   # Y_denominator <- tf$matrix_transpose(tf$expand_dims(Y_denominator, -1L))
#   # cos_theta <- numerator / X_denominator / Y_denominator
#   # jitter <- 1e-15
#   # theta <- tf$acos(jitter + 1.0 - 2.0 * jitter) * cos_theta
#   # K <-  sigma2 * (1.0 / pi) * calc_J(theta) * X_denominator ^ order * Y_denominator ^ order
# 
#   # calculate diagonal of kernel  
#   # Kdiag <- sigma2 * (1.0 / pi) * J(theta) * weighted_product(X, weight_sigma2, bias_sigma2) ^ order
#   
#   NULL
#   
# }
# 
# # periodic kernel (base kernel class)
# tf_periodic <- function(X, Y = NULL, sigma, period, lengthscales = NULL) {
#   
#   stop("periodic kernel is not implemented", call. = FALSE)
#   
#   # precalculate variance
#   # sigma2 <- sigma * sigma
#   
#   # add dummy dimension so we can use broadcasting
#   # f <- tf$expand_dims(X, -2L)  #  ... x N x 1 x D
#   # f2 <- tf$expand_dims(Y, -3L) # ... x 1 x M x D
# 
#   # calculate distances  
#   # r <- pi * (f - f2) / period
#   # r <- tf$reduce_sum(tf$square(tf$sin(r) / lengthscales), -1L)
#   
#   # construct and return periodic kernel  
#   # sigma2 * tf$exp(-0.5 * r)
# 
#   NULL
#   
# }

tf_Prod <- function(kernel_a, kernel_b) {
  
  tf$multiply(kernel_a, kernel_b)
  
}

tf_Add <- function(kernel_a, kernel_b) {
  
  tf$add(kernel_a, kernel_b)

}

# rescale, calculate, and return clipped squared distance
scaled_square_dist <- function(X, X_prime, lengthscales) {
  
  X <- X / lengthscales
  X_prime <- X_prime / lengthscales
  
  Xs <- tf$reduce_sum(tf$square(X), axis = -1L)
  Xs_prime <- tf$reduce_sum(tf$square(X_prime), axis = -1L)
  
  dist <- tf$constant(-2.0, dtype = tf$float64) * tf$matmul(X, X_prime, transpose_b = TRUE)  
  dist <- dist + tf$transpose(Xs)
  dist <- dist + Xs_prime
  
  # return value clipped around single float precision
  tf$maximum(dist, 1e-40)
  
}

# helper function to pull out slices of tensors and add final dim if dropped
tf_slice <- function(X, dims) {
  
  X <- tf$gather(X, dims, axis = -1L)
  
  if (length(dims) == 1)
    X <- tf$expand_dims(X, axis = -1L)
  
  X
  
}
# 
# # combine as module for export via internals
# tf_kernels_module <- module(tf_static,
#                             tf_white,
#                             tf_constant,
#                             tf_bias,
#                             tf_squared_exponential,
#                             tf_rational_quadratic,
#                             tf_linear,
#                             tf_polynomial,
#                             tf_exponential,
#                             tf_Matern12,
#                             tf_Matern32,
#                             tf_Matern52,
#                             tf_cosine,
#                             tf_arc_cosine,
#                             tf_periodic)
