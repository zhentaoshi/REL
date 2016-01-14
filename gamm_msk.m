function [L_hat, gam] = gamm_msk( b )
% b: the trial value of the coefficient
% tau: the tuning parameter
% b
global tau 
% the interpretatino is still the same. 

[hh, ~, ~] = MomentComponents(b);
n = size(hh,1);
m0 = size(hh, 2);
mm = m0;

% if strcmp(DGP, 'EKK')    
%     hh_mean = mean(hh);
%     naIndex = isnan(hh_mean);
%     hh = hh(:, ~naIndex );
%     mm = size(hh, 2);
% else
%     mm = m0;
% end
%% 


% the decision variables are the probablity 'p'
prob.opr = repmat('log', [n 1]);
prob.opri = zeros(n,1);
prob.oprj = (1:n)';
prob.oprf = ones(n,1);
prob.oprg = ones(n,1);

prob.c = sparse( zeros(n, 1) );
prob.a = [ ones(1,n); hh'] ; % combine the sum of prob into the linear constraint
prob.blc = [ 1; -tau*ones(mm, 1) ];
prob.buc = [ 1;  tau*ones(mm, 1) ];
prob.blx = sparse( zeros(n, 1) ); % lower bound of decision variables
prob.bux =  ones(n, 1); %Inf * ones(n, 1);

[res] = mskscopt(prob.opr, ...
    prob.opri, prob.oprj, prob.oprf, prob.oprg,...
    prob.c, prob.a, prob.blc, prob.buc, prob.blx, prob.bux,[], 'maximize echo(0)' );

 
gam = sparse( res.sol.itr.y/n );   
gam(1) = []; % the gamma, L-multiplier. The first one is removed   
%% infeasible case
if ~strcmp( res.sol.itr.solsta, 'OPTIMAL') || all( gam == 0 )
    L_hat = 1e+20; % Inf; %1e+20;
    gam = ones(m0, 1);
else 
    L_hat = -res.sol.itr.pobjval; % the optimal value of the function
%     
%     if strcmp(DGP, 'EKK')
%         gamNA = zeros(m0, 1);
%         gamNA( ~naIndex ) = gam;
%         gam = gamNA;
%     end
end


end

