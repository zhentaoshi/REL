% sup_score works much worse when useless = 20.

clear; 

 seed = 660;% 188 was a good seed. 188;
 rng(seed)

useful = 120; 
useless = 2;
global  n y x z m tau

%%
n = 100; % sample size
Rep = 50; % # of replication

C1 = 0.5; % constant in the tuning parameter
beta0 = [1;1];
dx = size(beta0, 1);
m = useful + useless + dx - 1;
tau =  C1 * sqrt( log(m) / n );  % tuning parameter


m = useful + useless + dx - 1;

rho = 0.6;

sup_s = @(b) max(abs(mean( g(b') ) ) ); 
EL_fval   = @(b) gamm_msk(b); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-5;
tolX = 1e-8;
optionsSearch = optimset('Display', 'notify', 'MaxIter', 1e+7, 'MaxFunEvals', 1e+4, 'TolFun', tol, 'TolX', tolX);

 optionsCon    = optimoptions('fmincon','Algorithm', 'trust-region-reflective', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'on');

 optionsKN    = optimoptions('ktrlink','Algorithm', 'interior-point', 'Display', 'off',...
     'FinDiffType', 'central', 'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,...
     'GradObj', 'on');


B1 = zeros(Rep, dx);
B2 = zeros(Rep, dx);
Prob = zeros(Rep, 1);

tic
r = 1;
while r <= Rep
    
  [y, x, z] = dgpLinearIVWeak(beta0, rho, useful);

    multiStart = 5;
    [xbase, ybase] = meshgrid( linspace(-3, 3, multiStart^2) );
    bSupMulti = zeros(multiStart^2, dx);
    fSupMulti = zeros(multiStart^2, 1);
    
    for i = 1:multiStart^2
        i;
      bInit = [xbase(i), ybase(i)];
     [bSupMulti(i, :), fSupMulti(i) ,exitflag ]  = fminsearch( sup_s, bInit, optionsSearch );
     if exitflag ~=1
         break
     end
    end
    
    toc

    %       [b_sup]  = fminsearch( sup_s, beta0', optionsSearch );
    if  exitflag == 1; % to make sure b_sup is in the compact set
        [~, which] = min( fSupMulti );
        b_sup = bSupMulti(which,:);
        
        if ~ any( abs(b_sup)  > 5 )
            b_init_EL = b_sup;
            
            [b_PEL] = ktrlink( EL_fval, b_init_EL,[],[],[],[],-5*ones(dx,1),5*ones(dx,1) ,[], optionsKN );
            [~, ~, ~, Prob(r)] = EL_fval(b_PEL);
            % b_PEL = b_sup;
            
            % exitFlag_PEL = 1;
            %% add the feature that if PEL fail, just let it be b_sup
            % if norm( b_sup - b_PEL, 2) > 5 || exitFlag_PEL ~= 1
            %     r = r - 1;
            % else
            B1(r, :) = b_sup;
            B2(r, :) = b_PEL;
            
            [b_sup b_PEL]
            r = r + 1;
        end
    end
    r
    toc
end
%%
valid =  all( abs(B1) < 5, 2);

B1 = B1(valid,:);
B2 = B2(valid, :);
Prob = Prob(valid);

[norms( bsxfun(@minus, B1, beta0'), 2, 2), norms( bsxfun(@minus, B2, beta0'), 2, 2)]


%%
ID = (1:size(B1, 1) )';
BB = dataset( ID, Prob, B1, B2);

title = ['BB_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', num2str(C1), '_seed_', num2str(seed) ];
export(BB, 'file', [title, 'BB.csv'], 'Delimiter',',');

