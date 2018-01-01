# this R scipt contains a function "REL" and a minimal example.



REL = function(X, tau){

  # INPUT:
  #   X: an (n*m) matrix, where n is the sample size and m is the number of moments
  #   tau: tuing parameter
  
  # OUTPUT:
  #   if the optimization is successful, the function returns
  #       $sol: the optimizer ( n-dimensional vector )
  #       $opt: the value of the criterion function ( \sum log pi )
  #   if the optimization fails, the function returns a warning message
  
  # DEPENDENCE:
  #   A solver called MOSEK. 
  #       download: https://www.mosek.com/downloads/
  #       academic license: https://www.mosek.com/products/academic-licenses/
  #   An R called package "Rmosek"
  #       installation instruction: https://docs.mosek.com/8.0/rmosek/install.html
  
  
  library(Rmosek)

  
  prob = list(sense = "max")

  prob$A = Matrix( rbind( 1,  t(X) ) ) # make a (m+1)*n matrix  
  prob$c = rep(0,n) 
  
  
  prob$bc = rbind( blc = c( 1, rep(-tau,m) ), buc = c( 1,  rep(tau,m) ) ) 
  prob$bx = rbind( blx = rep(0,n), bux = rep(1,n) )
  
 
  opro <- matrix ( list (), nrow =5, ncol = n )
  rownames ( opro ) <- c(" type ","j","f","g","h")
  for (i in 1:n){ opro[,i] = list("LOG", i, 1, 1, 0 ) }
  

  prob$scopt <- list ( opro = opro  )
    
  r = mosek ( prob, opts = list(soldetail = 2, verbose = 0) )

  opt = r$sol$itr$pobjval # primal objective value  
  sol = r$sol$itr$xx

#    if (r$response$msg == "MSK_RES_ERR_USER_NLO_FUNC: The user-defined nonlinear function reported an error."){
#      J = 1000}
  
  if ( abs( sum(sol) - 1 ) > 0.000001 ) { 
    return( warnings("solution fails!") )
  } else{
   return( list( sol = sol, opt = opt ) )
  }

}

#########################################
# a minimal example

n = 100;
m = 4;
tau = .01;
X = matrix( rnorm( n * m), nrow = n)

rel = REL(X, tau)

print(rel)

