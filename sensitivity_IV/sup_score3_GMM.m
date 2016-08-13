% clear;
% load('DGP_Han_BB_n_240_m_160_Rep_500_C1_0.5_seed_204.mat')
% rng(seed);

global  n y x z m tau DGP W_GMM

% DGP = 'Han';
tol = 1e-8;
tolX = 1e-10;
optionsSearch = optimset('Display', 'notify', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolFun', tol, 'TolX', tolX);
optionsCon    = optimset('Algorithm', 'active-set', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');
optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');


%%
tic
r = 1;
T = 10;
M_hat_seq = zeros(Rep, T);
B_GMM0 = zeros( size(B1) );
B_GMM1 = zeros(size(B1) );
% B31_std = B3_std;

while r <= Rep
    
    if strcmp(DGP, 'linearIV' )
	rng(seed_v_real(r))
        [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data
        C1 = 0.5; % constant in the tuning parameter
        tau =  C1 * sqrt( log(m) / n );  % tuning parameter
        lowerLimit = [-5; -5];
        upperLimit = [ 5;  5];

    elseif strcmp(DGP, 'Han' )
         y = dgpHan(beta0); % generate the data
         tau = Tau(r,1);
    end
    %%
    bGMMInit = B1(r,:)';
    b_GMM0 = ktrlink( 'GMM_Q0', bGMMInit,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsKN );
    [~ , gg_GMM0, ~ ] = MomentComponents(b_GMM0);
   
    W_GMM = cov(gg_GMM0); % optimal-weighting matrix 
    b_GMM  = ktrlink( 'GMM_Q', b_GMM0,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsKN );
    toc
    %%
    B_GMM0(r, :) = b_GMM0;
    B_GMM1(r, :) = b_GMM;
    
    [b_GMM0, b_GMM]    

    r = r + 1;
    r
end


title = ['GMM_DGP_', DGP, '_B_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1), '_seed_', num2str(seed) ];
save([title, '.mat']);

