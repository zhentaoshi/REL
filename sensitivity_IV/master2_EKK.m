% though it is indeed easy to do 5 parameters, that is too slow to afford.

global  n y m tau DGP dx tol tolX seed
global init Sn param type Ran W_GMM % particular for EKK

seed = 301;
% m = 80;
% n = 200;
Rep = 500; % of Monte Carlo replication

tol = 1e-4;
tolX = 1e-6;

DGP = 'EKK';

%% parameter
if strcmp(DGP, 'EKK')   
    beta0 = [2.46];
    param.sigA  = 0.69; % changed from 1.69. Gives much stable results.
    param.sigE  = 0.34; 
    param.rho   = -0.65;
    param.lamb  = 0.91;
    
    Sn = n * 5; % the number of artificial firms    
    % C1 = 1; % tuning paramter. The sample size is different from the other two examples
    tau =  C1 * sqrt( log(m) / n );  % tuning parameter
    lowerLimit = [  1.46; 0.5; 0.2; 0.2; -.9 ]; % enlarged the parameter space
    upperLimit = [  3.46; 3  ; 2.0; 2.0;  .9 ];
end
dx = length(beta0);

optionsSearch = optimset('Display', 'notify', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolFun', tol, 'TolX', tolX);
optionsCon    = optimset('Algorithm', 'active-set', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');
optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');

%% function handles
sup_s = @(b) max(abs(mean( MomentComponents(b) ) ) ); 
EL_fval   = @(b) gamm_msk(b); 

%% pre-define
B = zeros(Rep, 4);
GAM = zeros(Rep, m);
Prob = zeros(Rep, 3); 
%1st column: exitflag of the sup-score
%2nd column: exitflag of PEL
%3rd column: total prob of PEL
Tau = zeros(Rep, 3);

%%

r = 1;
tic
seedVar = seed; % seed change in each loop. so we can repeat it in other scripts.
for r = 1:Rep
    %% data generation
    seedVar = seedVar + 1;
    rng(seedVar)
    
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
    
    %% estimation
    
    % one-dimensional sup-score. no initial value needed.
    [b_sup] = fminbnd( sup_s, lowerLimit(1), upperLimit(1) );
    
    
    % GMM. 
    % GMM0 is the equally weighted GMM. use 'b_sup' as initial value
    % GMM is the 'optimally weighted' GMM. use 'GMM0' as initial value
    % use 'g' function, instead of 'h' function
    bGMMInit = b_sup;
    b_GMM0 = fmincon( 'GMM_Q0', bGMMInit,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsCon );
    [~ , gg_GMM0, ~ ] = MomentComponents(b_GMM0);
   
    W_GMM = cov(gg_GMM0); % optimal-weighting matrix 
    b_GMM  = fmincon( 'GMM_Q', b_GMM0,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsCon );
    toc
    
    % PEL
    % 'invalid initial value' almost never happens.
    bPELInit = b_sup;
    [b_PEL, L, PELExitFlag] = fmincon( EL_fval, bPELInit,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsCon );
    [~, gam] = gamm_msk(b_PEL);
    GAM(r, :) = gam';
    % display
    [ b_sup, b_GMM, b_GMM0, b_PEL]'
    r
    L
    
    B(r,:) = [ b_sup, b_GMM, b_GMM0, b_PEL ]';
    toc
end

title = ['DGP_', DGP, '_BB_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1) ];
save([title, '.mat']);
