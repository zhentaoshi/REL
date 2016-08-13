function [H, V, p, totalProb] = h_deriv_EKK( b, gam )

% b: the trial value of the coefficient
% gam: the estimated gamma under the given b in the inner loop

% this is the derivative of L wrt beta, not moments wrt beta

global  n tau NotnanIndex

[hh, gg, sd_g ]= MomentComponents(b);
NotnanIndex = ~isnan(hh(1,:));
hh0 = zeros(size(hh));
hh0(:, NotnanIndex) = hh(:, NotnanIndex);
% in MOSEK, when NaN in hh, the lagrangian multiplier is 0.
% hh is creaed to match this rule in the estimation of p.

p = 1./( n * ( 1 +  hh0 * gam - tau * norm(gam,1) ) );
totalProb = sum(p);
% always check if the prob adds up to 1. 
% If not, indicates serious problem, perhaps imported a wrong dataset.
hh = hh(:, NotnanIndex); % the short version with nan removed.
mm = size(hh, 2);

% %%
p_cell = num2cell(p);
epsilon = 1e-2; % 1e-3 was very good. for the complex method
hhh(:,:,1) = MomentComponents(b+epsilon);
hhh(:,:,2) = MomentComponents(b-epsilon);
H0 = ( hhh(:,:,1) - hhh(:,:,2) )/ (2*epsilon);
% H0 = CSDGrad( b, epsilon);

H0 = H0(:, NotnanIndex); % the short version with nan removed.
H   = p' * H0;
% max(abs(H))


%%
% a chunk of old code.
%     %% this chuck is special for 'EKK'
%     % numerical difference
%     epsilon = 1e-6; % 1e-5 is a good value. the choice of epsilon partly depends on the magnitude of the function.
%     eval1 = sum( EKK(b + epsilon) ); % must be a n * m
%     eval2 = sum( EKK(b - epsilon) );
%     Dg = -( eval1 - eval2 )/(2*epsilon); 
% 	% negative or positive log-likelihood doesn't matter.
% 	% always calculate Dg, the Derivative of g, in the form how g is written.
% 	% in this case, it is negative because E(g) =  E(y) - theta(b)
%     
%     sd_g = sd_g';
%     % each cell is for one 'i'    
%     term1 = Dg ./sd_g; % the signal is positive. Correct. -1 * -1 = 1;
%     
%     term1 = repmat( term1, [n 1]);
%     COV = zeros(1,mm);
%     for j = 1:mm
%         COV(j) = mean( gg(:,j).* Dg(:,j) ) - mean(gg(:,j)) * mean(Dg(:,j));
%     end
%     
%     term2 = 0.5 *  gg .* repmat( COV./ (sd_g.^(-3)), [n 1] ); % derivative wrt b1 
%     H  = term1 - term2; 
%     H  = p'*H;

    %% this chuck is commmon for all DGP.
    hh_cell = mat2cell(hh, ones(n, 1), mm);
    VV_cell = cellfun(@(p, hh) p * hh' * hh,  p_cell, hh_cell, 'UniformOutput', false);
    VV1 = reshape( cell2mat(VV_cell), [mm n mm]);
    VV2 = permute(VV1, [1 3 2]);
    V = sum(VV2, 3);
end




