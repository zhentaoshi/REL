function [lambda, dLambda] = varyingFixedEffect(b1, T)

bw = 3;
t = 1:T;

u = 1/bw * ( b1 * T - t );
lambda = 1 + (t/T) .* exp( -u.^2 ); % the lambda function
dLambda =    2/3 * t .* exp( -u.^2 ) .* (b1*T - t); % derivative of lambda
end