clear;
load('B_corr_n_100_m_75_Rep_500_C1_0.5_seed_101.mat')
rng(seed);

global  n y x z m DGP dx tau
DGP = 'linearIV';
tau = 0;

B4 = zeros( size(B3) );
B4_std = zeros(size(B3) );

Bo = B4;
Bo_std = Bo;

EL_fval   = @(b) gamm_msk(b); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-5;
tolX = 1e-8;
optionsSearch = optimset('Display', 'notify', 'MaxIter', 1e+7, 'MaxFunEvals', 1e+7, 'TolFun', tol, 'TolX', tolX);

 optionsCon    = optimset('Algorithm', 'trust-region-reflective', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');

 optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');


%%
tic
r = 1;

while r <= Rep
    
    if strcmp(DGP, 'linearIV' )
        [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data
        mm = m;
        
        M_hat = 1:4;
        m = length(M_hat);
        z = z(:, M_hat); % select the choosen moments
        
        [b_EL] = fmincon( EL_fval, beta0,[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
        [~, gam] = gamm_msk(b_EL);
        
        
        [H, V, totalProb] = h_deriv( b_EL, gam );
        H = H';    
        
        D = H' * V * H;
        Var_b = inv(D);
        
        b_EL_std = sqrt(n) * (b_EL - beta0) ./sqrt( diag(Var_b) );
        
        B4(r, :) = b_EL;
        B4_std(r, :) = b_EL_std;
        
        
        %% oracle
        
        [b_oracle] = ktrlink( EL_fval, beta0,[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
        [~, gam] = gamm_msk(b_oracle);
        [b_EL, b_oracle]
        
        [H, V, totalProb] = h_deriv( b_oracle, gam );
        H = H';    
        
        D = H' * V * H;
        Var_b = inv(D);
        
        b_oracle_std = sqrt(n) * (b_oracle - beta0) ./sqrt( diag(Var_b) );
        
        Bo(r, :) = b_oracle;
        Bo_std(r, :) = b_oracle_std;
        
        r = r + 1;
    end
    r
    m = mm;
    toc
end
title = ['B_oracle_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_seed_', num2str(seed) ];
save([title, '.mat'], 'B3_std', 'B4', 'B4_std', 'Bo', 'Bo_std', 'seed', 'DGP', 'Rep','useful', 'rho', 'n', 'm', 'beta0', 'dx');
export( dataset(B3_std, B4_std, Bo_std), 'file', [title, '.csv'], 'Delimiter', ',');

