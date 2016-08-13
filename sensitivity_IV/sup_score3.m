global  n y x z m tau DGP dx

Rep = 500; % # of Monte Carlo replication

tol = 1e-8;
tolX = 1e-10;

DGP = 'linearIV';

if strcmp( DGP, 'linearIV' )
    seed = 301;
    useful = 4; 
    rho = 0.6; % the endogeneity
    beta0 = [1;1]; % column vector. Good.
    C1 = 0.5; % constant in the tuning parameter
end

if strcmp( DGP, 'Han')
    seed = 201;
    beta0 = [0.9; -2];
    C1 = 0.5;
    lowerLimit = [ 0.2; -5];
    upperLimit = [1.01; 1 ];
end


tau =  C1 * sqrt( log(m) / n );  % tuning parameter
dx = length(beta0);

sup_s = @(b) max(abs(mean( MomentComponents(b) ) ) ); 
EL_fval   = @(b) gamm_msk(b); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optionsSearch = optimset('Display', 'notify', 'MaxIter', 1e+7, 'MaxFunEvals', 1e+7, 'TolFun', tol, 'TolX', tolX);
optionsCon    = optimset('Algorithm', 'trust-region-reflective', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');
 optionsKN    = optimset('Algorithm', 'interior-point', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');


B1 = zeros(Rep, dx);
B2 = zeros(Rep, dx);
GAM = zeros(Rep, m);
Prob = zeros(Rep, 3); 
%1st column: exitflag of the sup-score
%2nd column: exitflag of PEL
%3rd column: total prob of PEL
Tau = zeros(Rep, 3);
seed_v_real = zeros(Rep,1);

%%
tic
r = 1;
r_no_repeat = 0;
multiStart = 11;
while r <= Rep
    r_no_repeat = r_no_repeat + 1;
   %% linear IV 
    
   if strcmp(DGP, 'linearIV' )
       seed_v = seed + r_no_repeat;
       rng(seed_v);
       
    [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data

    [xbase, ybase] = meshgrid( linspace(-5, 5, multiStart^2) );
    bSupMulti = zeros(multiStart^2, dx);
    fSupMulti = zeros(multiStart^2, 1);
    
    % sup-score
    for i = 1:multiStart^2
        i;
      bInit = [xbase(i); ybase(i)];
     [bSupMulti(i, :), fSupMulti(i) ]  = fminsearch( sup_s, bInit, optionsSearch );
    end
    
    % toc
    [supMinValue, which] = min( fSupMulti );
    b_sup = bSupMulti(which,:)';
    
    if all( abs(b_sup) < 5 ) % to make sure b_sup is in the compact set
        seed_v_real(r) = seed_v;
        
        [b_PEL] = ktrlink( EL_fval, b_sup,[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
        [~, gam] = gamm_msk(b_PEL);
        GAM(r, :) = gam';
        [H, V, totalProb] = h_deriv( b_PEL, gam );
        
        B1(r, :) = b_sup;
        B2(r, :) = b_PEL;
        
        [b_sup b_PEL]        
        r = r + 1;
    end
   end
   
   
   %% Han's DGP
   
   if strcmp(DGP, 'Han')
       seed_v = seed + r_no_repeat;
       rng(seed_v);
       y = dgpHan(beta0); % generate the data
       [xbase, ybase] = meshgrid( linspace(lowerLimit(1), upperLimit(1), multiStart),...
           linspace(lowerLimit(2), upperLimit(2), multiStart) );
       bSupMulti = zeros(multiStart^2, dx);
       fSupMulti = zeros(multiStart^2, 1);
       supExitFlag = zeros(multiStart^2, 1);
       
       for i = 1:multiStart^2
           i;
           bInit = [xbase(i); ybase(i)];
           [bSupMulti(i, :), fSupMulti(i), supExitFlag(i) ]  = fminsearch( sup_s, bInit, optionsSearch );
       end
       
       [supMinValue, which] = min( fSupMulti );
       b_sup = bSupMulti(which,:)';
       Prob(r, 1) = supExitFlag(which);
       
       g_sort = sort( (abs(mean( MomentComponents(b_sup) ) ) ), 'descend');
       
       if b_sup(1) <= upperLimit(1) && b_sup(1) >= lowerLimit(1) &&...
               b_sup(2) <= upperLimit(2) && b_sup(2) >= lowerLimit(2) % to make sure b_sup is in the compact set
           %        tau = g_sort(C1);
           Tau(r, 1) = tau;
           bPELInit = b_sup;
           [L0, gam] = gamm_msk(bPELInit);
           
           if L0 < 1e+10 % if the initial value is infeasible
               seed_v_real(r) = seed_v;
               
               [b_PEL, ~, PELExitFlag] = ktrlink( EL_fval, bPELInit,[],[],[],[],...
                   lowerLimit, upperLimit,[], optionsKN );
               
               Prob(r, 2) = PELExitFlag;
               
               [~, gam] = gamm_msk(b_PEL);
               GAM(r, :) = gam';
               
               Prob(r, 3) = probSum(b_PEL, gam);
               
               [b_sup b_PEL]
               
               B1(r, :) = b_sup;
               B2(r, :) = b_PEL;
               
               r = r + 1;
           end
       end
   end
    r
    toc
end
%%
% B1
% B2
% 
% [norms( bsxfun(@minus, B1, beta0'), 2, 2), norms( bsxfun(@minus, B2, beta0'), 2, 2)]
ID = (1:Rep)';
BB = dataset( ID, Prob, B1, B2);


title = ['DGP_', DGP, '_BB_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1) ];
save([title, '.mat']);

% export(BB, 'file', [title, 'BB.csv'], 'Delimiter',',');

