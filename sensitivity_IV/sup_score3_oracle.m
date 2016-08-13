
% load('D:\My Box Files\data_EL_sim\DGP_EKK_B_corr_n_1200_m_80_Rep_500_C1_1_seed_301.mat')
% load('D:\My Box Files\data_EL_sim\DGP_EKK_B_corr_n_1200_m_120_Rep_500_C1_1_seed_301.mat')
% load('D:\My Box Files\data_EL_sim\DGP_EKK_B_corr_n_1800_m_80_Rep_500_C1_1seed_301.mat')
load('D:\My Box Files\data_EL_sim\DGP_EKK_B_corr_n_1800_m_120_Rep_500_C1_1_seed_301.mat')


global  n y x z m DGP dx tau M_hat
tau = 0;

B4 = zeros( size(B3) );
B4_std = zeros(size(B3) );

Bo = B4;
Bo_std = Bo;

EL_fval   = @(b) gamm_msk(b); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-5;
tolX = 1e-8;
% optionsSearch = optimset('Display', 'notify', 'MaxIter', 1e+7, 'MaxFunEvals', 1e+7, 'TolFun', tol, 'TolX', tolX);

% optionsCon    = optimset('Algorithm', 'trust-region-reflective', 'Display', 'off',...
%    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'on');

 optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');

%%
tic
r = 1;
seedVar = seed;
while r <= Rep
    
    if strcmp(DGP, 'linearIV' )
        rng(seed);
        [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data
        mm = m;
        
        M_hat = M_hat_seq(r, 1:ceil(n^(1/3) ) );
        m = length(M_hat);
        z = z(:, M_hat); % select the choosen moments
        
        
        [b_EL] = ktrlink( EL_fval, beta0,[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
        [~, gam] = gamm_msk(b_EL);
        
        
        [H, V, totalProb] = h_deriv( b_EL, gam );
        H = H';    
        
        D = H' * V * H;
        Var_b = inv(D);
        
        b_EL_std = sqrt(n) * (b_EL - beta0) ./sqrt( diag(Var_b) );
        
        B4(r, :) = b_EL;
        B4_std(r, :) = b_EL_std;
        
        
        % oracle
        M_hat = 1:4;
        m = length(M_hat);
        z = z(:, M_hat); % select the choosen moments
        
        [b_oracle] = ktrlink( EL_fval, bete0,[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
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
    
        %%
    if strcmp(DGP, 'Han' )
        %% Han's seed has to be redo
        y = dgpHan(beta0); % generate the data
        mm = m;
        
        % post etimation
        M_hat = M_hat_seq(r, 1:ceil(n^(1/3) ) );
        
        [b_EL] = ktrlink( EL_fval, beta0,[],[],[],[],lowerLimit,upperLimit ,[], optionsKN );
        [~, gam] = gamm_msk(b_EL);
        
        
        [H, V, totalProb] = h_deriv( b_EL, gam );
        H = H';    
        
        D = H' * V * H;
        Var_b = inv(D);
        
        b_EL_std = sqrt(n) * (b_EL - beta0) ./sqrt( diag(Var_b) );
        
        B4(r, :) = b_EL;
        B4_std(r, :) = b_EL_std;
        
        
        % oracle
        beta01 = beta0(1);
        M_hat = (round( beta01 * m) - 2):(round( beta01 * m) +2 );
        
        [b_oracle] = ktrlink( EL_fval, beta0,[],[],[],[],lowerLimit,upperLimit ,[], optionsKN );
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
    %%
    if strcmp(DGP, 'EKK' )
        seedVar = seedVar + 1; % varying seed
        rng(seedVar)
        lowerLimit = lowerLimit(1);
        upperLimit = upperLimit(1);
        
        M_hat = []; % cannot touch the random data generation!!! generate the full data
        init.N_nF = [ n,... % home country
            ceil( n * (0.6 + 0.2 * random('Uniform', 0, 1, 1, 6)) ) ,...
            ceil( n * (0.3 + 0.1 * random('Uniform', 0, 1, 1, m - 6 )) ) ]; % +2 ,... % strong moments
        init.X_bar = [1,  random('Uniform', 0, 1, 1, m) ];
        
        type = 'real';
        Ran.v = random('Uniform', 0, 1, n, 1);
        Ran.a = random('Normal',  0, 1, n , m+1);
        Ran.h = random('Normal',  0, 1, n , m+1);
        
        y  = n * EKK(beta0); % the real data;
        
        type = 'art'; % change to the 'artificial status'
        Ran.v = random('Uniform', 0, 1, Sn, 1);
        Ran.a = random('Normal',  0, 1, Sn , m+1);
        Ran.h = random('Normal',  0, 1, Sn , m+1);

        mm = m; % teporally save m
        
        % post-estimation
        M_hat = M_hat_seq(r, 1:ceil(n^(1/5) ) );
        y_full  = y;
        y   = y_full(:, M_hat);
        
        [b_EL] = ktrlink( EL_fval, beta0,[],[],[],[],lowerLimit,upperLimit ,[], optionsKN );
        [~, gam] = gamm_msk(b_EL);
        
        
        [H, V, totalProb] = h_deriv_EKK( b_EL, gam );
        H = H';    
        
        D = H' * V * H;
        Var_b = inv(D);
        
        b_EL_std = sqrt(n) * (b_EL - beta0) ./sqrt( diag(Var_b) );
        
        B4(r, :) = b_EL;
        B4_std(r, :) = b_EL_std;
        
        
        % oracle
        beta01 = beta0(1);
        M_hat =  1:ceil(n^(1/5) );
        y = y_full(:, M_hat);
        
        [b_oracle] = ktrlink( EL_fval, beta0,[],[],[],[],lowerLimit,upperLimit ,[], optionsKN );
        [~, gam] = gamm_msk(b_oracle);
        [b_EL, b_oracle]
        
        [H, V, totalProb] = h_deriv_EKK( b_oracle, gam );
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
title = ['DGP_', DGP, '_B_oracle_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_seed_', num2str(seed) ];
save([title, '.mat']);
export( dataset(B3_std, B4_std, Bo_std), 'file', [title, '.csv'], 'Delimiter', ',');

