% in linear model, it is exactly the Bai and Ng.

function [IC] = IC(M, b)

[s, V] = score_linearIV(gg, sd_g, hh); % score function with
v_M = hessian(M, b); % Hessian function with a subset M

IC = s_M' * pinv( v_M ) * s_M;
end


