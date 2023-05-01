# bmi510.R

#' Rando: Sample dataframes and vectors.
#'
#' @description 
#' This function checks if the input x is an atomic vector or a dataframe-like
#' object, and applies the appropriate operation to obtain a sample of size n.
#'
#' @param x Either an atomic vector or dataframe-like object.
#' @param n Number of samples or rows to be obtained from the input ().
#' @param replace Logical value that indicates if sample is done with replacement.
#' 
#' @return An atomic vector or dataframe-like object of size n.
#' 
#' @export
#' 
#' @examples 
#' av = c(1,2,3,4,5,6,7,9)
#' rando(av, n = 5, replace = FALSE)
#' 
#' df = data.frame(x = c(1,2,3,4,5,6,7,8,9), y = c(9,8,7,6,5,4,3,2,1))
#' rando(df, n = 5, replace = TRUE)
rando = function(x,n=1,replace=T) {
  if (is.atomic(x)) {
    return (sample(x, size = n, replace = replace))
  }
  else if (is.data.frame(x)) {
    return (dplyr::sample_n(x, size = n, replace = replace))
  }
}


#' is_min: Minimum value in atomic vector.
#'
#' @description 
#' This function takes an atomic vector x as an input, finds the minimum value, and
#' returns a logical vector with TRUE where x equals its minimum value, and FALSE
#' for all other positions.
#'
#' @param x Atomic vector.
#' @param na.rm Logical value that indicates if NA values should be removed from x.
#' 
#' @return A logical vector with TRUE where x equals its minimum value.
#' 
#' @export
#' 
#' @examples 
#' av = c(1,2,3,4,5,6,7,9)
#' is_min(av, na.rm=TRUE)
is_min = function(x,na.rm=T) {
  if (is.atomic(x)) {
    return (x == min(x, na.rm = na.rm))
  }
}


#' is_max: Minimum value in atomic vector.
#'
#' @description 
#' This function takes an atomic vector x as an input, finds the maximum value, and
#' returns a logical vector with TRUE where x equals its maximum value, and FALSE
#' for all other positions.
#'
#' @param x Atomic vector.
#' @param na.rm Logical value that indicates if NA values should be removed from x.
#' 
#' @return A logical vector with TRUE where x equals its maximum value.
#' 
#' @export
#' 
#' @examples 
#' av = c(1,2,3,4,5,6,7,9)
#' is_max(av, na.rm=TRUE)
is_max = function(x,na.rm=T) {
  if (is.atomic(x)) {
    return (x == max(x, na.rm = na.rm))
  }
}


#' rep_mat: Replicates matrix in row and column dimension.
#'
#' @description 
#' This function takes a matrix or a dataframe x as an input, converts it to 
#' a matrix, and returns a matrix with M replicates of x in the row dimension
#' an dN replicates of x in the column dimension.
#'
#' @param x Matrix or dataframe.
#' @param M Integer that indicates the replicates of x in the row dimension
#' @param N Integer that indicates the replicates of x in the column dimension
#' 
#' @return A matrix with replicates of x in the row and the column dimension.
#' 
#' @export
#' 
#' @examples 
#' x = matrix(1:25,ncol=5)
#' rep_mat(x, M = 3, N = 2)
rep_mat = function(x, M=1, N=1) {
  if (is.matrix(x) || is.data.frame(x)) {
    x = as.matrix(x)
    x = matrix(rep(t(x), M), ncol = ncol(x), byrow = TRUE)
    x = matrix(rep(x, N), nrow = nrow(x), byrow = FALSE)
    return (x)
  }
}


#' classes: Classes of variables in a tibble.
#'
#' @description 
#' This function takes a tibble x as an input, finds the class of each variable, 
#' and returns a character vector with the corresponding classes.
#'
#' @param x Tibble.
#' 
#' @return A character vector with the identified classes.
#' 
#' @export
#' 
#' @examples 
#' tib = dplyr::tibble(x = c("a", "b"), y = c(1, 2))
#' classes(tib)
classes = function(x) {
  if (tibble::is_tibble(x)) {
    return (unlist(dplyr::summarize_all(x, dplyr::funs(class)), use.names = FALSE))
  }
}


#' df_scale: Scaled numerical variables.
#'
#' @description 
#' This function takes a tibble x as an input, finds the class of each variable, 
#' and returns a tibble in which the numeric variables have been scaled.
#'
#' @param x Tibble.
#' @param center Logical value that indicates if the columns should be centered.
#' @param scale Logical value that indicates if the columns should be scaled
#' 
#' @return Scaled tibble.
#' 
#' @export
#' 
#' @examples 
#' tib = dplyr::tibble(x = c("a", "b", "c"), y = c(1, 2, 3), z = c(10, 2, 7))
#' df_scale(tib, center = TRUE, scale = TRUE)
df_scale = function(x, center = T, scale = T) {
  if (tibble::is_tibble(x)) {
    classes_vector = classes(x)
    x[, classes_vector == "numeric"] = scale(x[, classes_vector == "numeric"], center = center, scale = scale)
    return (x)
  }
}


#' log_likelihood_norm: Log-likelihood normal distribution.
#'
#' @description 
#' This function takes a vector x as an input, and returns the Log-likelihood
#' under the normal distribution with the corresponding parametrizations.
#'
#' @param x Numeric vector (sample).
#' @param mean Mean of the distribution.
#' @param sd Standard deviation of the distribution.
#' 
#' @return Numeric value - Log-likelihood.
#' 
#' @export
#' 
#' @examples 
#' x = rnorm(100, mean = 0, sd = 1)
#' log_likelihood_norm(x = x, mean = 0, sd = 1)
log_likelihood_norm = function(x, mean, sd) {
  return (sum(log(dnorm(x = x, mean = mean, sd = sd))))
}


#' log_likelihood_unif: Log-likelihood uniform distribution.
#'
#' @description 
#' This function takes a vector x as an input, and returns the Log-likelihood
#' under the uniform distribution with the corresponding parametrizations.
#'
#' @param x Numeric vector (sample).
#' @param min Lower limit of the distribution.
#' @param max Upper limit of the distribution.
#' 
#' @return Numeric value - Log-likelihood.
#' 
#' @export
#' 
#' @examples 
#' x = runif(100, min = 0, max = 10)
#' log_likelihood_unif(x = x, min = 0, max = 10)
log_likelihood_unif = function(x, min, max) {
  return (sum(log(dunif(x = x, min = min, max = max))))
}


#' log_likelihood_chisq: Log-likelihood chi-square distribution.
#'
#' @description 
#' This function takes a vector x as an input, and returns the Log-likelihood
#' under the chi-squared distribution with the corresponding parametrizations.
#'
#' @param x Numeric vector (sample).
#' @param df Degrees of freedom.
#' 
#' @return Numeric value - Log-likelihood.
#' 
#' @export
#' 
#' @examples 
#' x = rchisq(100, df = 4)
#' log_likelihood_chisq(x = x, df = 4)
log_likelihood_chisq = function(x, df) {
  return (sum(log(dchisq(x = x, df = df))))
}


#' log_likelihood_f: Log-likelihood F distribution.
#'
#' @description 
#' This function takes a vector x as an input, and returns the Log-likelihood
#' under the F distribution with the corresponding parametrizations.
#'
#' @param x Numeric vector (sample).
#' @param df1 Degrees of freedom.
#' @param df2 Degrees of freedom.
#' 
#' @return Numeric value - Log-likelihood.
#' 
#' @export
#' 
#' @examples 
#' x = rf(100, df1 = 2, df2 = 5)
#' log_likelihood_f(x = x, df1 = 2, df2 = 5)
log_likelihood_f = function(x, df1, df2) {
  return (sum(log(df(x = x, df1 = df1, df2 = df2))))
}


#' log_likelihood_t: Log-likelihood t distribution.
#'
#' @description 
#' This function takes a vector x as an input, and returns the Log-likelihood
#' under the t distribution with the corresponding parametrizations.
#'
#' @param x Numeric vector (sample).
#' @param df Degrees of freedom.
#' 
#' @return Numeric value - Log-likelihood.
#' 
#' @export
#' 
#' @examples 
#' x = rt(100, df = 2)
#' log_likelihood_t(x = x, df = 2)
log_likelihood_t = function(x, df) {
  return (sum(log(dt(x = x, df = df))))
}


#' sensitivity: Sensitivity given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' sensitivity.
#'
#' @param pred Logical or numerical vector of binary predictions.
#' @param truth Logical or numerical vector of binary ground truth values
#' 
#' @return Numeric value - sensitivity.
#' 
#' @export
#' 
#' @examples 
#' sensitivity(c(1,1,0,1,0), c(1,0,0,1,0))
sensitivity = function(pred,truth) {
  if (length(pred) == length(truth)) {
    tp = sum(pred & truth)
    fn = sum(!pred & truth) 
    return(tp / (tp + fn))
  }
}


#' specificity: Specificity given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' specificity
#'
#' @param pred Logical or numerical vector of binary predictions.
#' @param truth Logical or numerical vector of binary ground truth values
#' 
#' @return Numeric value - specificity.
#' 
#' @export
#' 
#' @examples 
#' specificity(c(1,1,0,1,0), c(1,0,0,1,0))
specificity = function(pred,truth) {
  if (length(pred) == length(truth)) {
    tn = sum(!pred & !truth)
    fp = sum(pred & !truth)
    return(tn / (tn + fp))
  }
}


#' precision: Precision given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' precision
#'
#' @param pred Logical or numerical vector of binary predictions.
#' @param truth Logical or numerical vector of binary ground truth values
#' 
#' @return Numeric value - precision.
#' 
#' @export
#' 
#' @examples 
#' precision(c(1,1,0,1,0), c(1,0,0,1,0))
precision = function(pred,truth) {
  if (length(pred) == length(truth)) {
    tp = sum(pred & truth)
    fp = sum(pred & !truth) 
    return(tp / (tp + fp))
  }
}


#' recall: Recall given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' recall
#'
#' @param pred Logical or numerical vector of binary predictions.
#' @param truth Logical or numerical vector of binary ground truth values
#' 
#' @return Numeric value - recall
#' 
#' @export
#' 
#' @examples 
#' recall(c(1,1,0,1,0), c(1,0,0,1,0))
recall = function(pred,truth) {
  if (length(pred) == length(truth)){
    tp = sum(pred & truth)
    fn = sum(!pred & truth) 
    return(tp / (tp + fn))
  }
}


#' accuracy: Accuracy given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' accuracy
#'
#' @param pred Logical or numerical vector of binary predictions.
#' @param truth Logical or numerical vector of binary ground truth values.
#' 
#' @return Numeric value - accuracy.
#' 
#' @export
#' 
#' @examples 
#' accuracy(c(1,1,0,1,0), c(1,0,0,1,0))
accuracy = function(pred,truth) {
  if (length(pred) == length(truth)) {
    tp = sum(pred & truth)
    tn = sum(!pred & !truth)
    return((tp + tn) / length(pred))
  }
}


#' f1: F1-Score given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' F1-Score
#'
#' @param pred Logical or numerical vector of binary predictions.
#' @param truth Logical or numerical vector of binary ground truth values
#' 
#' @return Numeric value - F1-Score.
#' 
#' @export
#' 
#' @examples 
#' f1(c(1,1,0,1,0), c(1,0,0,1,0))
f1 = function(pred,truth) {
  if (length(pred) == length(truth)) {
    tp = sum(pred & truth)
    fn = sum(!pred & truth) 
    fp = sum(pred & !truth) 
    return(tp / (tp + (1/2) * (fp + fn)))
  }
}


#' minimum_n_per_group: Sample size.
#'
#' @description 
#' This function takes the desired delta, d, and the desired power, power,  and 
#' returns the number of observations per group using the power calculation for
#' a two sample t-test.
#'
#' @param d Numeric value of true difference in means.
#' @param power Numeric value of the power of the test.
#' 
#' @return Numeric value - rounded number of samples per group.
#' 
#' @export
#' 
#' @examples 
#' minimum_n_per_group(d = 0.5, power = 0.8)
minimum_n_per_group = function(d,power = 0.8) {
  test = power.t.test(n = NULL, delta = d, power = power, type = "two.sample")
  return(ceiling(test$n))
}


#' r2: R squared given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' R squared
#'
#' @param pred Continuous numerical vector of predictions.
#' @param truth Continuous numerical vector of ground truth values
#' 
#' @return Numeric value - R Squared.
#' 
#' @export
#' 
#' @examples 
#' r2(c(2.6,2.2,4.2,3.9,6.0), c(2.3,3.2,4.7,3.9,3.3))
r2 = function(pred,truth) {
  if (length(pred) == length(truth)) {
    return (cor(truth, pred) ** 2)
  }
}


#' adj_R2: Adjusted R squared given predictions and ground truth.
#'
#' @description 
#' This function takes two vectors, pred and truth, as inputs, and returns the 
#' adjusted R squared.
#'
#' @param pred Continuous numerical vector of predictions.
#' @param truth Continuous numerical vector of ground truth values
#' @param truth Numeric value that indicates the number of model parameters, excluding the intercept.
#' 
#' @return Numeric value - Adjusted R Squared.
#' 
#' @export
#' 
#' @examples 
#' adj_R2(c(2.6,2.2,4.2,3.9,6.0), c(2.3,3.2,4.7,3.9,3.3), n_p = 2)
adj_R2 = function(pred,truth,n_p) {
  if (length(pred) == length(truth)) {
    rsq = r2(pred, truth) 
    n = length(pred)
    return (1 - (((1 - rsq) * (n - 1)) / (n - n_p - 1)))
  }
}