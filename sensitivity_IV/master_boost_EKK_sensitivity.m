


global  n m y  DGP seed 
global  Sn type Ran NotnanIndex % particular for EKK
% Rep = 10;
B2 = B(:,4); % PEL estimator


% %% don't touch this set of parameter. This is good for a known beta0
%      beta0 = 5;
%      param.lamb  = 1;
%      param.sigA  = .3; % changed from 1.69. Gives much stable results.
%      param.sigE  = .3; 
%      param.rho   = -0.5;
%     

tic

T = 15;
M_hat_seq = zeros(Rep, T);
B3 = zeros( size(B2,1),6 );
B3_std = zeros(size(B3) );
ProbCheck = zeros(Rep,1);

seedVar = seed; % seed change in each loop. so we can repeat it in other scripts.

for r = 1:Rep
    r
    seedVar = seedVar + 1; % varying seed
    rng(seedVar)
    
    
    type = 'real';
    Ran.v = random('Uniform', 0, 1, n, 1);
    Ran.a = random('Normal',  0, 1, n , m+1);
    Ran.h = random('Normal',  0, 1, n , m+1);
    
    y  = n * EKK(beta0); % the real data;
    
    type = 'art'; % change to the 'artificial status'
    Ran.v = random('Uniform', 0, 1, Sn, 1);
    Ran.a = random('Normal',  0, 1, Sn , m+1);
    Ran.h = random('Normal',  0, 1, Sn , m+1);
    
    b_PEL = B2(r);
    gam = GAM(r,:)';
    
    b_PEL = beta0;    
    gam = 1/n * zeros(size(GAM(1,:)'));
    [H, V, p, ProbCheck(r)] = h_deriv_EKK( b_PEL, gam );
    H = H';


    %% boosting
    
    for t = 1:T
        if t == 1
            M0 = [];
            [M_hat] = boost_EKK(M0, H, V);
            M0 = M_hat;
        else
            [M_hat] = boost_EKK(M0, H, V);
            M0 = sort( [M0 M_hat] ); % augment the index
        end
        
        if ~isempty(M_hat)
            M_hat_seq(r, t) = M_hat;
        else
            M_hat_seq(r, t) = t;
        end
    end
    % originally, 
    % M_hat = M_hat_seq(r, 1:(round(  n^(1/5)) ) ); % M_hat has removed the NaN.
    
    %% bias correction
    for jj = 1:10
        M_hat = M_hat_seq(r, 1:jj);
        hh = MomentComponents(b_PEL);
        hh = hh(:, NotnanIndex);
        h_bar_p = sum( bsxfun(@times, hh, p), 1)';

        sandwich = H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * H(M_hat, :);
        Var_b = inv(sandwich);
        b_corr = b_PEL - sandwich \ ( H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * h_bar_p(M_hat, :) );


       % [b_PEL, b_corr];
        B3(r,jj) = b_corr;
        B3_std(r, jj) = sqrt(n) * (b_corr - beta0)./sqrt( diag( Var_b ) );

        sum( abs(B3_std)>1.96 ) /r;
    end

end

RMSE = sqrt( mean( (B3 - 5).^2 ) );
BIAS = mean( B3 - 5);
out = [RMSE; BIAS]

%%
title = ['DGP_sensitivity_', DGP, '_B_corr_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1), '_seed_', num2str(seed) ];
% save([title, '.mat']);
export(dataset(out), 'file', [title, '.csv'], 'Delimiter', ',');

