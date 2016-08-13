clear
global  n y m tau dx seed DGP dd %#ok<NUSED>
global  Sn Ran   

dx = 5;

load('trade_data.mat')

rng(seed)
y  = sell; % the real data of entry.

Ran.v = random('Uniform', 0, 1, Sn, 1);
Ran.a = random('Normal',  0, 1, Sn , m+1);
Ran.h = random('Normal',  0, 1, Sn , m+1);

C1 = 1;
tau = C1 * sqrt( log(m)/n );

%%
tic
r = 1;
T = 10;
M_hat_seq = zeros(1, T);
M_hat_total = zeros(dx, T);


    [~, gam ]   = gamm_msk(b_PEL);
    [H, V, p, totalProb] = h_deriv_EKK_5( b_PEL, gam );
    H = H';
 
    
    %% boosting
    
    for dd = 1:5    
        for t = 1:T
            if t == 1
                M0 = [];
                [M_hat] = boost_EKK(M0, H, V);
                M0 = M_hat;
            else
                [M_hat] = boost_EKK(M0, H, V);
                M0 = sort( [M0 M_hat] ); % augment the index
            end
            M_hat_seq(r, t) = M_hat;
        end
        M_hat_total(dd,:) = M_hat_seq(r, : );
    end
    
    M_hat = unique( M_hat_total(:, 1:ceil((n/log(m))^(0.2)) ) )

    
    %% bias correction
    [hh, gg, sd_g] = MomentComponents(b_PEL);
    h_bar_p = sum( bsxfun(@times, hh, p), 1)';
    
    sandwich = H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * H(M_hat, :);
    Var_b = inv(sandwich);
    b_corr = b_PEL - sandwich \ ( H(M_hat,:)' * pinv( V(M_hat, M_hat) ) * h_bar_p(M_hat, :) );
    
    variance = sqrt(diag(Var_b))/sqrt(n);
    
    b_corr
    variance



