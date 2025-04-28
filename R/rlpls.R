#' A main function of RLPLS
#'
#' @description
#' This is a main function of sparse partial least square regression with random LASSO (RLPLS).
#'
#' @param X Matrix of predictors.
#' @param Y Vector or matrix of responses.
#' @param maxk The number of hidden components, with default 5.
#' @param eta0 A thresholding parameter between 0 and 1, with default 0.4.
#' @param kppa0 A parameter between 0 and 0.5 managing the effect of the concavity of the objective function and the closeness of direction vectors and surrogate vectors. It affects a model for multivariate response variables. The default is 0.3.
#' @param Btimes The number of boostrapping for the first direction vector.
#' @param b_var_num A vector indicating the number of selected variables for bootstrapping. The length should be 2 and each element needs to be smaller than the number of predictor variables.
#' @param var.selection_alpha A significance level for variable selection.
#' @param model_fit A PLS algorithm for fitting a model. The alternatives are \code{"kernelpls"}, \code{"widekernelpls"}, \code{"simpls"}, or \code{"oscorespls"}. The default is \code{"simpls"}. The selection is independent from the selection of \code{variable_select}.
#' @param variable_select A PLS algorithm for variable selection. The alternatives are \code{"pls2"} or \code{"simpls"}. The default is \code{"simpls"}. The selection is independent from the selection of \code{model_fit}.
#' @param scale.x logical value whether to scale predictor variables.
#' @param scale.y logical value whether to scale response variables.
#' @param eps0 An effective zero, with default 1e-4.
#' @param maxiter The maximum number of iterations when fitting direction vectors, with default is 5000.
#' @param seed A seed number.
#'
#' @return a list of rlpls result.
#'
#' @examples
#' set.seed(111)
#' # Set the size of samples and true coefficients
#' n = 30
#' beta0 = c(runif(20,-10,-5),runif(20,5,10),rep(0,10))
#'
#' # Construct predictors & response
#' # the structure is similar with that from Chun and Keles (2010)
#' xxx = yyy = c(NULL)
#' H1 = c(rep(2,times=25),rep(4,times=n-25))
#' H2 = rep(10,n)
#' for(i in 1:length(beta0)){
#'   if(i <= 30){
#'     xxx = cbind(xxx, H1+rnorm(n,0,1))
#'   }else{
#'     xxx = cbind(xxx, H2+rnorm(n,0,1))
#'   }
#' };rm(i)# for i
#'
#' yyy = xxx %*% beta0 + rnorm(n,0,1.5^2)
#'
#' # Run RLPLS
#' library(rlpls)
#' result = rlpls(X = xxx, Y = yyy, b_var_num = c(25,20), Btimes = 5000 , variable_select = "pls2")
#' result$beta.rlpls
#'
#' @export
rlpls = function(X, Y, maxk = 5, eta0 = 0.4, kppa0 = 0.3,
                   Btimes = 1000, b_var_num = c(NA, NA), var.selection_alpha = 0.05,
                   model_fit = "simpls", variable_select = "simpls",
                   scale.x = TRUE, scale.y = TRUE,
                   maxiter = 5000, eps0 = 1e-4, seed = 2022){
  set.seed(seed)

  ## Step 0 : SET initial values ----
  p = ncol(X); q = ncol(Y); N = nrow(X)
  X = matrix(X, nrow=N, ncol=p)

  Y0 = Y
  X0 = X

  # scale X (variance)
  if(scale.x){
    vec1 = matrix(1, nrow = 1, ncol = N)
    sig.x = sqrt((vec1 %*% X0^2)/(N - 1))
    if (any(sig.x < .Machine$double.eps)) stop("Some columns of the matrix X have zero-variance.")
    X0 = scale(X0, center = FALSE, scale = sig.x)
  }else{
    sig.x = rep(1, p) # for 'return'
  } # scale X (variance)

  # scale Y
  if(scale.y) {
    vec1 = matrix(1, nrow = 1, ncol = N)
    sig.y = sqrt((vec1 %*% Y0 ^ 2) / (N - 1))
    if (any(sig.y < .Machine$double.eps))  stop("Some columns of the matrix Y have zero-variance.")
    Y0 = scale(Y0, center = FALSE, scale = sig.y)
  } else{
    sig.y = rep(1, q) # for 'return'
  } # scale.y


  beta.hat = matrix(0,p,q, dimnames = list(paste0("X",1:p), paste0("Y",1:q)))
  beta.mat = list()

  # initial dataset
  X1 = X0
  Y1 = Y0

  for(k in 1:maxk){

    ## Step 1 : FIND 1st vector ----
    ### rlasso/Step1 : Boostrap Samples ----------------------------------------
    b1_w1.list = matrix(0, nrow = Btimes, ncol = p) # zero for unselected variables
    b2_w1.list = matrix(0, nrow = Btimes, ncol = p) # zero for unselected variables
    for(bth in 1:Btimes){

      ## 1) Draw B(Btimes) bootstrap samples with size n by sampling with replacement
      b_id = sample(1:N, N, replace = TRUE)

      ## 2) For the b1-th bootstrap sample, randomly select q1(=b_var_num) candidate variables
      b_var = sort(sample(1:p, b_var_num[1], replace = FALSE))
      b_X1 = as.matrix(X1[b_id,b_var])
      b_Y1 = as.matrix(Y1[b_id,])

      ## 3) Find 1st vector
      b1_w1.list[bth,b_var] = sparse1st(Z = t(b_X1) %*% b_Y1,
                                        eta = eta0, kappa = kppa0, eps = eps0, maxstep = maxiter)

    } # for bth (bootstrap)

    ## 4) Compute the importance measure
    var_imp = abs(apply(b1_w1.list,2,mean))

    ### rlasso/Step2 : Selecting Variables ----------------------------------------

    for(bth in 1:Btimes){
      ## 1) Draw another set of B(Btimes) bootstrap samples with size n
      ##    by sampling with replacement from the original training data set.
      b2_id = sample(1:N, N, replace = TRUE)

      ## 2) For the b2th bootstrap sample, randomly select q2 candidate variables
      ##    with selection probability of xj proportional to its importance Ij
      ##    obtained in Step 1
      b2_var = sort(sample((1:p)[var_imp != 0], min(b_var_num[2], sum(var_imp != 0)), prob = (var_imp/sum(var_imp))[var_imp != 0]))
      b2_X1 = as.matrix(X1[b2_id,b2_var])
      b2_Y1 = as.matrix(Y1[b2_id,])

      ## 3) Adaptive lasso -> obtain the estimator b_beta
      b2_w1.list[bth,b2_var] = sparse1st(Z = t(b2_X1) %*% b2_Y1,
                                         eta = eta0, kappa = kppa0, eps = eps0, maxstep = maxiter)
    } # for bth (bootstrap)

    ## 4) Compute the final estimator beta from b_beta
    ## with parametric testing (Park et al.,2015)
    if(var.selection_alpha > 0){
      tmp.w = apply(b2_w1.list,2,mean)

      # the average of the selection ratio, var.selection_alpha = 0.05
      hat.pi = sum(b2_w1.list != 0)/(nrow(b2_w1.list)*ncol(b2_w1.list))
      C_j = apply(b2_w1.list,2,function(x){sum(x != 0)})
      hat.w = as.numeric(C_j >= qbinom(1-var.selection_alpha, Btimes, hat.pi))
    }else{
      hat.w = apply(b2_w1.list,2,mean)
    }
    setA  = unique((1:p)[hat.w != 0 | beta.hat[,1] != 0])

    # SELECT variables w.r.t hat.w
    X.a = X0[, setA, drop = FALSE]

    ## Step 2 : PLS_direction vector ----------------------------------------

    # RUN PLS with reduced.X
    fit.pls = pls::plsr(Y1 ~ X.a, ncomp = min(k, length(setA)), method = "simpls", scale = FALSE)

    beta.hat = matrix(0, nrow = p, ncol = q)
    beta.hat[setA, ] = matrix(fit.pls$coefficients[,,min(k, length(setA)), drop=FALSE], nrow = length(setA), ncol = q)
    beta.mat[[k]] = beta.hat

    # Step 3 : UPDATE X1 or Y1 depending on variable_selection method --------
    if(variable_select == "pls2"){
      Y1 = Y0 - X1 %*% beta.hat # X1 = X0
    }else if(variable_select == "simpls"){
      P.a = fit.pls$projection
      X1 = X0
      X1[,setA] = X.a - X.a%*%P.a%*%solve(t(P.a)%*%P.a)%*%t(P.a) # update raw data
    } # if :: variable_select

  } # for k

  return(list(beta.rlpls  = beta.hat, dataset = list(x = X, y = Y), pls.fit.method = model_fit, variable.select.method = variable_select))
}
