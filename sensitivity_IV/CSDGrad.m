function [ vCSDGrad ] = CSDGrad( x0, dx )

nParams = length( x0 ) ;
vCSDGrad = zeros( nParams, 1 ) ;
xPlus = x0 + 1i * dx ;

for ix = 1 : nParams
    x1 = x0 ;
    x1( ix ) = xPlus( ix ) ;
    [ fval ] = MomentComponents( x1 ) ;
    vCSDGrad = imag( fval / dx ) ;
end