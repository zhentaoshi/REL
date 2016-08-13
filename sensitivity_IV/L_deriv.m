function [S, totalProb] = L_deriv( b, gam )

% b: the trial value of the coefficient
% gam: the estimated gamma under the given b in the inner loop

% this is the derivative for beta.
% in the current version (8/15/2013) I do not use it.
% I compared the case of 'linearIV' and 'Han', in which we have closed-form derivative.
% the difference is negligible whether we provide the solver with the analytic derivative.

global dx y x z n m DGP tau

[hh, gg, sd_g] = MomentComponents(b);

stdGamma = gam./sd_g;
p = 1./( n * ( 1 +  gg * stdGamma - tau * norm(gam,1) ) );
totalProb = sum(p) 
% always check if the prob adds up to 1. 
% If not, indicates serious problem, perhaps imported a wrong dataset.
if strcmp(DGP, 'linearIV')
   % standarization by the std error

    
    X_cell = mat2cell(x, ones(n, 1), dx);
    Z_cell = mat2cell(z, ones(n, 1), m);
    p_cell = num2cell(p);
    
    xz_cell = cellfun(@(x,z) x'*z, X_cell, Z_cell, 'UniformOutput', false);
    xzBar = cellfun(@(x) reshape(x, [1 dx * m]), xz_cell, 'UniformOutput', false);
    xzBar = mean( cell2mat(xzBar) );
    xzBar = reshape(xzBar, [dx m]);
    
    gg_cell = mat2cell(gg, ones(n, 1), m);
    ggBar = mean(gg); 
    
    
    
    sdDeriv_cell = cellfun(@(gg, xz) 0.5 * gg .* sd_g'.^(-3) * gam * (gg - ggBar) * ( -(xz -xzBar))',...
        gg_cell, xz_cell, 'UniformOutput', false);
    sdDeriv = mean( cell2mat(sdDeriv_cell ) );
%%    
    score_i = cellfun(@(p, x, z, gg) ( n * p * ( -(x' * z)) * stdGamma - n * p *  sdDeriv' )',...
        p_cell, X_cell, Z_cell, gg_cell, 'UniformOutput', false );
    score_i = cell2mat(score_i);
    
    S = - sum( score_i); % "-" because of negative loglikelihood
    
elseif strcmp(DGP, 'Han')
    S = 0;
end

end

