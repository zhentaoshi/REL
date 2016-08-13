function [y] = EKK(b)

global n m type Sn % M_hat


%%
if strcmp(type, 'real')
    S = n;
elseif strcmp(type, 'art')
    S = Sn;
end

mm_hat = 1:m;

[ delta_nF, X_nF, u_bar ] = EKK_simulator(S, mm_hat, b);


N_nF  = sum(delta_nF);
% X_bar = sum( X_nF )./N_nF;
y = X_nF./repmat( N_nF, [S 1]);
y(:,1) = [] ; % remove the home country

%% notation
% theT = theta_tilde
% lamb   = lambda
% rho  = rho
% sigA = sigma_alpha
% sigE = sigma_eta

% k1 = kappa_1
% k2 = kappa_2

        function [ delta_nF, X_nF, u_bar ] = EKK_simulator(S, mm_hat, b)
        % Eaton, Kortum and Kramarz (2011). 7 steps.
            theT = b(1);
            % lamb = b(2);
            % sigA = b(3);
            % sigE = b(4);
            % rho  = b(5);
        

        global init param Ran

       sigA = param.sigA;
       sigE = param.sigE;
       rho  = param.rho;
       lamb = param.lamb;
        
        a    = Ran.a;
        h    = Ran.h;
        v    = Ran.v;

         X_bar = init.X_bar;
         N_nF  = init.N_nF;
         
        % a = a(:, [0, mm_hat] + 1 ); % +1 because of the home country.
        % h = h(:, [0, mm_hat] + 1 );
        
         
        % X_bar = X_bar(:, [0, mm_hat] + 1 );
        % N_nF  = N_nF (:, [0, mm_hat] + 1 );
        %% step 0
        % in the paper, S is number of firms. s is the firm index. n is the index of countrys
        % in this script, I replace S by n, replace n by m

        %%  step 1
        % Eqaution (34)
        k1 = (theT / (theT - 1) - theT / (theT + lamb - 1) ) *...
            exp( 0.5 * ( sigA + 2* rho * sigA * sigE * (theT -1) + sigE * (theT - 1)^2 )  ) ;
        % Equation (35)
        k2 = exp( 0.5 * (theT * sigE)^2 );

        %% step 2
        % X_bar is directly calculated from the real data.
        sig_E_nF = k2/k1 * X_bar;

        %% step 3
        alpha = exp( sigA * sqrt( 1 - rho^2) * a + sigA * rho *  h );
        alpha(alpha>100) = 100; alpha(alpha<0.01) = 0.01;
        
        eta   = exp( sigE * h );
        eta(eta>100) = 100; eta(eta<0.01) = 0.01;
        
        ratio = alpha./eta;

        %% step 4 construct the entries hurdles
        % N_nF is from the data
        component0 = repmat( N_nF, [S 1 ]) .* eta;
        u_n_bar = component0.^theT; % /k1 in the paper. Not needed at all.

        %% step 5
        % let the first country be the home country.
        u_bar_X = max( u_n_bar(:, 2:length(mm_hat)),[], 2 ); % S * 1 vector
        u_bar_F = u_n_bar(:, 1);
        u_bar   = min( [ u_bar_X, u_bar_F] , [], 2 );

        %% step 6
        u = v .* u_bar; % n * 1 vector

        %% step 7 empty

        %% step 8
        delta_nF = ( repmat(u, [1 length(mm_hat)+1] ) <= u_n_bar ) ; % n * m matrix
        
        component1 = repmat(u, [1 length(mm_hat)+1])  ./ u_n_bar ;
        component2 = repmat(sig_E_nF, [S, 1]);
        X_nF0 = ratio .* (1-component1).^ ( lamb / theT)  .*...
            component1.^(-1/theT)  .* component2;
        X_nF = delta_nF .* X_nF0;

        end

end
