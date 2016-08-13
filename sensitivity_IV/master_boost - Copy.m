% clear;

%load('F:\My Box Files\data_EL_sim\BB_n_100_m_75_Rep_500_C1_0.5_seed_101.mat')
%load('F:\My Box Files\data_EL_sim\BB_n_100_m_150_Rep_500_C1_0.5_seed_102.mat')
%load('F:\My Box Files\data_EL_sim\BB_n_200_m_75_Rep_500_C1_0.5_seed_103.mat')
% load('F:\My Box Files\data_EL_sim\BB_n_200_m_150_Rep_500_C1_0.5_seed_104.mat')

% load('F:\My Box Files\data_EL_sim\DGP_Han_BB_n_120_m_80_Rep_500_C1_0.5.mat')
% load('F:\My Box Files\data_EL_sim\DGP_Han_BB_n_120_m_160_Rep_500_C1_0.5.mat')
% load('F:\My Box Files\data_EL_sim\DGP_Han_BB_n_240_m_80_Rep_500_C1_0.5.mat')
% load('F:\My Box Files\data_EL_sim\DGP_Han_BB_n_240_m_160_Rep_500_C1_0.5.mat')
% rng(seed);

global  n y x z m tau DGP dx

% DGP = 'Han';

%%
tic
r = 1;
T = 10;
M_hat_seq = zeros(Rep, T);
B3 = zeros( size(B2) );
B3_std = zeros(size(B3) );
% B31_std = B3_std;
ProbCheck = zeros(Rep,1);

for r = 1:Rep
    rng(seed_v_real(r));
    r
    if strcmp(DGP, 'linearIV' )
        rng(seed_v_real(r))
        [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data
        C1 = 0.5; % constant in the tuning parameter
        tau =  C1 * sqrt( log(m) / n );  % tuning parameter

    elseif strcmp(DGP, 'Han' )
         y = dgpHan(beta0); % generate the data
         tau = Tau(r,1);
    end
    
    b_PEL = B2(r,:)';
    gam   = GAM(r, :)';
    % [H, V, p, totalProb] = h_deriv( beta0, gam );
    [H, V, p, ProbCheck(r)] = h_deriv( b_PEL, gam );
    H = H';
    
    
    %% boosting
    
    for t = 1:T
        if t == 1
            M0 = [];
            [M_hat] = boost(M0, H, V);
            M0 = M_hat;
        else
            [M_hat] = boost(M0, H, V);
            M0 = sort( [M0 M_hat] ); % augment the index
        end
        M_hat_seq(r, t) = M_hat;
    end
    M_hat = M_hat_seq(r, 1:ceil( dx * n^(1/5))  );
%     if strcmp(DGP, 'Han')
%         M_hat = M_hat_seq(r, 1:ceil( n^(1/5) ) );
%     end
    
    %% bias correction
    [hh, gg, sd_g] = MomentComponents(b_PEL);
    h_bar_p = sum( bsxfun(@times, hh, p), 1)';
    
    sandwich = H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * H(M_hat, :);
    Var_b = inv(sandwich);
    b_corr = b_PEL - sandwich \ ( H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * h_bar_p(M_hat, :) );
    
    
    [b_PEL, b_corr];
    B3(r,:) = b_corr;
    B3_std(r, :) = sqrt(n) * (b_corr - beta0)./sqrt( diag( Var_b ) );
    
    %%
end

if strcmp(DGP, 'linearIV')
title = ['DGP_', DGP, '_rho_', num2str(rho), '_B_corr_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1), '_seed_', num2str(seed) ];
else
    title = ['DGP_', DGP, '_B_corr_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1), '_seed_', num2str(seed) ];
end
save([title, '.mat']);

