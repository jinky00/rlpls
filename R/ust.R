#'
#' Univariate Soft Thresholding Estimator
#'
ust = function( b, eta ){
  b.ust = matrix( 0, length(b), 1 )
  if ( eta < 1 ){
    bb = abs(b) - eta * max( abs(b) )
    b.ust[ bb>=0 ] = bb[ bb>=0 ] * (sign(b))[ bb>=0 ]
  }
  return(b.ust)
}
