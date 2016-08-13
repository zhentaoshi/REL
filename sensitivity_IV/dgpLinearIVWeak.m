function [y, X, Z] = dgpLinearIVWeak(b, rho, useful)

%rho: magnitude of endogeneity
global m n

K = length(b);

%% instruments and disturbances
Z1 = random('Normal', 0, 1, n, m );
E = random('Normal', 0, 1, n, 2);
Sig = [1, rho; rho, 1];
E = E * Sig;
e = E(:,1);
v = E(:,2);

%% structural equation

x1 = 1 +  Z1(:,1:useful) * 0.1/sqrt(useful) * ones(useful, 1) + 0.5 * v;
x2 = random('Normal', 0, 1, n, K - 1);
X = [x1, x2];
Z = [Z1, x2];

y = X * b + 0.5 * e;

end
