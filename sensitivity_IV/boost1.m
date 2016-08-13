clear; clc;

T = 100;
K = 50;
eff = 5;
R = 8;

[y, X, Z] = DGP(T, K, eff);
[yt, Xt, Zt] = DGP(T, K, eff);

%% initial boosting
Q_vec = zeros( K , 1);
for k = 1:K 
    [b , b_tild, Q_vec(k)] = IVsub(k, X, Z, y, K, T, 0);
end
[Q, S0] = min(Q_vec); % note I is the index in term of length(M)
Q; % print Q
S0;

%% 


Q = zeros(1, R);
Qt= zeros(1, R);
Q2 = zeros(1,R);
B = zeros(K,R);
B_tild = zeros(K,R);
S = B;

for r = 1:R
    [B(:,r), B_tild(:,r), S, Q2(r)] = boost0(S0, X, Z, y, K, T);
    disp('Q2')
    Q2(r);
    Qt(r) = ( 1/T * Zt'*(yt - Xt*B(:,r) ) )' * ( 1/T * Zt'*(yt - Xt*B(:,r)) );
    S0 = S
end
%%
Q2
plot(1:R, ([Q2', Qt']));
% %% compare Q2(b2) - Q3(b2)
% ID = 1:K;
% Qd = zeros(R-1, 2);
% for r = 2:(R-1)
%     Qa = QS( B(:,r),  X,Z,y,T);
%     Qb =  QS(  B_tild(:,r+1), X,Z,y,T) ;
%     Qd(r,1) = Qa - Qb;
%     
%     Qa = QS(  B_tild(:,r+1),  X,Z,y,T);
%     Qb =  QS(  B(:,r+1), X,Z,y,T); 
%     Qd(r,2) = Qa - Qb;
% end
% Qd
% Qdd = Qd(:,2) + Qd(:,1)
%% boostReg


%% lasso
% bl = IVlasso(50, X, Z, y, K )