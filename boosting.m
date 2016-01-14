clear; 

 seed = 190;% 188 was a good seed. 188;
 rng(seed)

useful = 5; 
useless = 500;
global  n y x z m tau

%%
n = 500; % sample size
Rep = 50; % # of replication

C1 = 0.5; % constant in the tuning parameter
beta0 = [1;1];
dx = size(beta0, 1);
m = useful + useless + dx - 1;
tau =  C1 * sqrt( log(m) / n );  % tuning parameter

rho = 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sup_s = @(b) max(abs(mean( g(b') ) ) ); 
EL_fval   = @(b) gamm_msk(b); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-5;
tolX = 1e-8;
optionsSearch = optimset('Display', 'notify', 'MaxIter', 1e+7, 'MaxFunEvals', 1e+7, 'TolFun', tol, 'TolX', tolX);

 optionsCon    = optimoptions('fmincon','Algorithm', 'trust-region-reflective', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'on');

 optionsKN    = optimoptions('ktrlink','Algorithm', 'interior-point', 'Display', 'off',...
     'FinDiffType', 'central', 'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,...
     'GradObj', 'on');


r = 1;
T = 10; % T is the number of maximum number of boosting iterations
M_hat_seq = zeros(Rep, T);
while r <= Rep    
  [y, x, z] = dgpLinearIV(beta0, rho, useful);
  
  % [b_PEL] = ktrlink( EL_fval, beta0',[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
  b_PEL = beta0';
  
  
  
  for t = 1:T
      if t == 1
          M0 = [];
          [M_hat] = boost0(M0, b_PEL);
          M0 = M_hat;
      else
          [M_hat] = boost0(M0, b_PEL);
          M0 = sort( [M0 M_hat] ); % augment the index
      end
      M_hat_seq(r, t) = M_hat;
  end
  
   r = r + 1;
end
%%


% title = ['BB_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', num2str(C1), '_seed_', num2str(seed) ];
% export(BB, 'file', [title, 'BB.csv'], 'Delimiter',',');

