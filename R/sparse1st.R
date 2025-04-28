#'
#' A function to find the first direction vector
#'
sparse1st = function( Z, eta, kappa, eps, maxstep){

  p = nrow(Z)
  q = ncol(Z)
  Z = Z / median( abs(Z) ) #Znorm1=median( abs(Z) )

  # main iterations
  # if univariate response, then just soft thresholding
  if ( q == 1 ){
    w = ust( Z, eta )
  }else if ( q > 1 ){

    # if multivariate response

    M = Z %*% t(Z)
    dist = 10
    i = 1

    # main iteration: optimize c and w iteratively
    # use svd solution if kappa==0.5
    if ( kappa==0.5 ){

      # initial value for a & c (outside the unit circle)
      w = rep(10, p) #matrix( 10, p, 1 )
      w.old = w

      while ( dist>eps & i<=maxstep ){
        # optimize w for fixed c
        mcsvd = svd( M%*%w )
        tmp = mcsvd$u %*% t(mcsvd$v)

        # optimize c for fixed w
        # soft thresholding ( assuming lambda2 -> Inf )

        w = ust( M%*%tmp, eta )

        # calculate discrepancy between a & c

        dist = max( abs( w - w.old ) )
        w.old = w
        i = i + 1
      }

      # solve equation if 0<kappa<0.5
    }else if( kappa>0 & kappa<0.5 ){

      kappa2 = ( 1 - kappa ) / ( 1 - 2*kappa )

      # initial value for c (outside the unit circle)
      w = rep(10, p) #matrix( 10, p, 1 )
      w.old = w

      # define function for Lagrange part
      h = function(lambda){
        alpha = solve( M + lambda*diag(p) ) %*% M %*% w
        obj = t(alpha) %*% alpha - 1/kappa2^2             # 1/kappa2^2 is correct (typo in paper)
        return(obj)
      }

      # control size of M & w if too small
      if( h(eps) * h(1e+30) > 0 ){
        while( h(eps) <= 1e+5 ){M = 2*M; w = 2*w}
      }

      while( dist>eps & i<=maxstep ){

        # control size of M & w if too small
        if ( h(eps) * h(1e+30) > 0 ){
          while( h(eps) <= 1e+5 ) {
            M = 2*M; w = 2*w
          } # while
        } # if

        # optimize a for fixed c
        lambda.st = uniroot( h, c( eps, 1e+30 ) )$root
        tmp = kappa2 * solve( M + lambda.st * diag(p) ) %*% M %*% w

        # optimize c for fixed w
        # soft thresholding ( assuming lambda2 -> Inf )
        w = ust( M%*%tmp, eta )

        # calculate discrepancy between a & c
        dist = max( abs( w - w.old ) )
        w.old = w
        i = i + 1
      } # while
    } # if : kappa
  } # if : q

  return(w)
}
