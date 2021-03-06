% clear;

% load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.5.mat')
% load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.5.mat')
% load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.5.mat')
% load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.5.mat')

% load('DGP_Han_BB_n_240_m_80_Rep_500_C1_0.5.mat')

% rng(seed);

global  n y x z m tau DGP dx

% DGP = 'Han';

%%
tic
r = 1;
T = 10;
M_hat_seq = zeros(Rep, T ,2);
B3 = zeros( size(B2,1), 10 );
B3_std = zeros(size(B3) );
% B31_std = B3_std;
ProbCheck = zeros(Rep,1);

for r = 1:Rep
    rng(seed_v_real(r));
    r
    if strcmp(DGP, 'linearIV' )
        rng(seed_v_real(r))
        [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data
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
    
    %% boosting
    for ii = 1:2
        for t = 1:T
            if t == 1
                M0 = [];
                [M_hat] = boost(M0, H(:,ii), V);
                M0 = M_hat;
            else
                [M_hat] = boost(M0, H(:,ii), V);
                M0 = sort( [M0 M_hat] ); % augment the index
            end
            
           M_hat_seq(r, t, ii) = M_hat;
        end
    end
    % M_hat = M_hat_seq(r, 1:ceil( dx * n^(1/5))  );
    
    %% bias correction
    for jj = 1:10
        M_hat = unique( M_hat_seq( r, 1:jj, :) );
        [hh, gg, sd_g] = MomentComponents(b_PEL);
        h_bar_p = sum( bsxfun(@times, hh, p), 1)';
        
        sandwich = H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * H(M_hat, :);
        Var_b = inv(sandwich);
        b_corr = b_PEL - sandwich \ ( H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * h_bar_p(M_hat, :) );
        
        
        [b_PEL, b_corr];
        B3(r,jj) = b_corr(1);
    end
    
    %%
end

if strcmp(DGP, 'linearIV')
    RMSE = sqrt( mean( (B3 - 1).^2 ) );
    BIAS = mean( B3 - 1);
    out = [RMSE; BIAS]
    title = ['DGP_', DGP, '_rho_', num2str(rho), '_B_corr_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
        num2str(C1), '_seed_', num2str(seed) ];
elseif strcmp(DGP, 'Han')
    RMSE = sqrt( mean( (B3 - .9).^2 ) );
    BIAS = mean( B3 - .9);
    out = [RMSE; BIAS]
    title = ['DGP_', DGP, '_B_corr_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
        num2str(C1), '_seed_', num2str(seed) ];
end
save([title, '.mat']);
export(dataset(out), 'file', [title, '.csv'], 'Delimiter', ',');

