
seed = 800;
% rng(seed)

global  n y x z m tau DGP dx boot

% m = 20;
% n = 100;
% Rep = 100;
MCRep = 120; % # of Monte Carlo replication

DGP = 'Han';

if strcmp( DGP, 'linearIV' )
    useful = 4; 
    rho = 0.6; % the endogeneity
    beta0 = [1;1]; % column vector. Good.
    beta_hypo = mul * beta0;
    C1 = 0.5; % constant in the tuning parameter
end

if strcmp( DGP, 'Han')
    beta0 = [.9; -2];
    beta_hypo =  mul * beta0;
    C1 = 0.5;
end

tau =  C1 * sqrt( log(m) / n );  % tuning parameter
dx = length(beta0);

sup_s = @(b) max(abs(mean( MomentComponents(b) ) ) ); 
EL_fval   = @(b) gamm_msk(b); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% optionsSearch = optimset('Display', 'notify', 'MaxIter', 1e+7, 'MaxFunEvals', 1e+7, 'TolFun', tol, 'TolX', tolX);
% 
%  optionsCon    = optimset('Algorithm', 'trust-region-reflective', 'Display', 'off',...
%     'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');
% 
%  optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off',...
%     'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');

GAMBoot = zeros(MCRep, m);
GG = zeros(m, MCRep);
L_all = zeros(MCRep, 1);

%% linear IV
tic

% [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data

decision = zeros(1, Rep);

for r = 1:Rep
    if strcmp(DGP, 'linearIV' )
        boot = 0;
        [d0.y, d0.x, d0.z] = dgpLinearIV(beta0, rho, useful); % generate the data
        % d0.* is the real data
        y = d0.y; x = d0.x; z = d0.z;
        testStat = LM(beta_hypo);
              
        % bootstrap
        boot = 1;
        GB = zeros(MCRep, m);
        gammB = zeros(MCRep, m);
        bootStat0 = zeros(MCRep, 1);
        for MCr = 1:MCRep           
            bootsample = randsample(1:n, n, 1);
            y = d0.y(bootsample,:);
            x = d0.x(bootsample,:);
            z = d0.z(bootsample,:);
            [bootStat0(MCr), GB(MCr, :), gammB(MCr, :) ] = LM(beta_hypo);
        end
        
        


    elseif strcmp(DGP, 'Han')
        boot = 0;
        d0.y = dgpHan(beta0); % generate the data
        y = d0.y;
        testStat = LM(beta_hypo);
              
        % bootstrap
        boot = 1;
        GB = zeros(MCRep, m);
        gammB = zeros(MCRep, m);
        bootStat0 = zeros(MCRep, 1);
        for MCr = 1:MCRep           
            bootsample = randsample(1:n, n, 1);
            y = d0.y(bootsample,:);
            [bootStat0(MCr), gammB(MCr, :), GB(MCr, :) ] = LM(beta_hypo);
        end
    end
    
%% back to the real world
%         bootStat1 = diag( ( GB - repmat( mean(GB), [MCRep, 1]) )' * gammB  );
        adj =  gammB * mean(GB)';
        bootStat1 = bootStat0 - adj;
        q95 = quantile( bootStat1, 0.95);
        if (testStat > q95)
            decision(r) = 1;
        end
    
    toc
    display( sprintf( 'testStat = %f', testStat ) )
    display( sprintf( 'q95 = %f', q95) )
    display( sprintf( 'iteration = %i', r) )
    display( sprintf( 'decision = %f', sum(decision)/r ) )
    display( sprintf( 'p_value = %f' , mean( testStat > bootStat1 ) ) )
    display( sprintf( ' ' ) )
        
end

toc



