% though it is indeed easy to do 5 parameters, that is too slow to afford.

global  n y m tau DGP dx tol tolX seed
global init Sn param type Ran W_GMM % particular for EKK
% load('market_condition.mat')
%% market condition
% don't touch. 1.5, 0.3 is very good parameter.
fixed.N_nF = [ 1, 1.5 * ones(1,6), 0.3 * ones(1,m-6)];
fixed.X_bar = ones(1,m+1); 
%%
seed = 301;

%% the market condition

    init.N_nF = ceil( n * fixed.N_nF(1: (m+1)) );
    init.X_bar = fixed.X_bar(1 : (m+1) );
    
% Rep = 10; % of Monte Carlo replication

tol = 1e-4;
tolX = 1e-6;

DGP = 'EKK';

%% parameter
if strcmp(DGP, 'EKK')   
    
     beta0 = 5;
     param.lamb  = 1;
     param.sigA  = .3; % siaA = 0.3 gave unbiased estimator. I am trying to give it some bias.
     param.sigE  = .3; 
     param.rho   = -0.5;

    
    Sn = n * 5; % the number of artificial firms    
    C1 = 0.5; % tuning paramter. The sample size is different from the other two examples
    tau =   C1 * sqrt( log(m) / n );  % tuning parameter
    lowerLimit =   beta0 -1 ;
    upperLimit =   beta0 +1 ;
end
dx = length(beta0);

optionsSearch = optimset('Display', 'notify', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolFun', tol, 'TolX', tolX);
optionsCon   = optimset('Algorithm', 'active-set', 'Display', 'off', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');
optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');

%% function handles

%%

r = 1;
tic
seedVar = seed; % seed change in each loop. so we can repeat it in other scripts.
for r = 1:Rep
    %% data generation
    seedVar = seedVar + 1;
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
    
    %% estimation

%     GMM. 
%     GMM0 is the equally weighted GMM. use 'b_sup' as initial value
%     GMM is the 'optimally weighted' GMM. use 'GMM0' as initial value
%     use 'g' function, instead of 'h' function
    b_sup = B(r,1);
    b_PEL = B(r,4);
    bGMMInit = b_sup;
    b_GMM0 = fmincon( 'GMM_Q0', bGMMInit,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsCon );
    [~ , gg_GMM0, ~ ] = MomentComponents(b_GMM0);
   
    W_GMM = cov(gg_GMM0); % optimal-weighting matrix 
    b_GMM  = fmincon( 'GMM_Q', b_GMM0,[],[],[],[],...
        lowerLimit(1:dx), upperLimit(1:dx),[], optionsCon );
    toc

    

    
    [ b_sup, b_GMM, b_GMM0, b_PEL]'
    r
    
    B(r,2:3) = [  b_GMM0, b_GMM ];
    toc
end

title = ['full_DGP_', DGP, '_BB_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1), '_seed_', num2str(seed) ];
save([title, '.mat']);
