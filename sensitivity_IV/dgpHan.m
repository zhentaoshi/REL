function [y] = dgpHan(b)
global n m 
b1 = b(1);
b2 = b(2);

lambda = varyingFixedEffect(b1, m);
lambda = [1, lambda]; % add the coefficient of the initial value

alpha = 1+ random('Normal', 0,1,n,1);
e = 1 * random('Normal', 0, 1, n, m+1);

y = repmat(lambda, [n 1]) .* repmat(alpha, [1 m+1]) + b2 + e;
% y includes the initial value
end
