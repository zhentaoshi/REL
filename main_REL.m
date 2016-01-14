global  n y m tau dx seed 
global  Sn Ran

load('trade_data.mat')
seed = 10;
rng(seed);

y  = sell; % the real data of entry.

Ran.v = random('Uniform', 0, 1, Sn, 1);
Ran.a = random('Normal',  0, 1, Sn , m+1);
Ran.h = random('Normal',  0, 1, Sn , m+1);

C1 = 1; % the tuning parameter
tau = C1 * sqrt( log(m)/n );

%% multiple start to determine the initial value
R = 10;
% set the number of initial values. In practice should be a large number, say R = 1000

BInit = zeros(dx, R);
validInit = zeros(R,1);
L0Init = zeros(R, 1);

tic
for r = 1:R
    rng(r)
    br = random('Uniform', lowerLimit, upperLimit, dx, 1);
    BInit(:, r) = br;
    [L0, gam_r] = gamm_msk(br); % check if it is a valid starting point
    L0Init(r) = L0;
    if L0 < 1e+10 && ~all( gam_r == 0 )
        validInit(r) = 1;
    end
    [r, sum(validInit)];
    toc
end

[~, init_index] = min(L0Init);
bPELInit = BInit(:, init_index);

%% REL estimation
optionsCon    = optimset('Algorithm', 'active-set', 'Display', 'off',...
    'FunValCheck', 'on', 'MaxIter', 1000, 'TolFun', tol, 'TolX', tolX ,'GradObj', 'off');

[b_PEL, L, PELExitFlag] = fmincon( EL_fval, bPELInit,[],[],[],[],...
    lowerLimit(1:dx), upperLimit(1:dx),[], optionsCon );
[~, gam] = gamm_msk(b_PEL);
b_PEL;
