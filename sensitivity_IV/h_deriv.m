function [H, V, p, totalProb] = h_deriv( b, gam )

% b: the trial value of the coefficient
% gam: the estimated gamma under the given b in the inner loop

% this is the derivative of L wrt beta, not moments wrt beta
global dx y x z n m DGP tau 

[hh, gg, sd_g] = MomentComponents(b);
p = 1./( n * ( 1 +  gg * (gam./sd_g) - tau * norm(gam,1) ) );
totalProb = sum(p);
% always check if the prob adds up to 1. 
% If not, indicates serious problem, perhaps imported a wrong dataset.
mm = size(hh, 2);


%%
p_cell = num2cell(p);
gg_cell = mat2cell(gg, ones(n, 1), mm);

if strcmp(DGP, 'linearIV')

    X_cell = mat2cell(x, ones(n, 1), dx);
    Z_cell = mat2cell(z, ones(n, 1), m);    
    
    xz_cell = cellfun(@(x,z) x'*z, X_cell, Z_cell, 'UniformOutput', false);
    xzBar = cellfun(@(x) reshape(x, [1 dx * m]), xz_cell, 'UniformOutput', false);
    xzBar = mean( cell2mat(xzBar) );
    xzBar = reshape(xzBar, [dx m]);    
    
    ggBar = mean(gg); 
    
    % each cell is for one 'i'    
    term1_cell = cellfun(@(x, z)   - x' * (z ./sd_g' ), X_cell, Z_cell, 'UniformOutput', false );
    
    term3_cell = cellfun(@(gg, xz) bsxfun(@times, (gg - ggBar),  -(xz -xzBar)) ,...
        gg_cell, xz_cell, 'UniformOutput', false);
    term3 = reshape( cell2mat(term3_cell), [ dx n m]);
    term3 = permute(term3, [ 1 3 2]);
    term3  =  mean(term3, 3); % NEGTIVE? think again
        
    term2_cell = cellfun(@(gg, xz) bsxfun(@times, 0.5 * (sd_g.^(-3))' .*  gg, term3)  ,gg_cell, xz_cell, 'UniformOutput', false);

%%    
    H0 = cellfun(@(p, term1, term2)   p * ( term1 -  term2 ) ,...
        p_cell, term1_cell, term2_cell, 'UniformOutput', false );
    H1 = reshape( cell2mat(H0), [ dx n m]);
    H2 = permute(H1, [ 1 3 2]);
    H  = sum(H2, 3); % NEGTIVE? think again
end
    
if strcmp(DGP, 'Han')
    % numerical derivative is less error-pro.
    epsilon1 = 0.000001;
    epsilon2 = 0.01;
    
    H0 = zeros(n, m, dx);
    
    hh11 = MomentComponents([b(1) + epsilon1; b(2)]);
    hh12 = MomentComponents([b(1) - epsilon1; b(2)]);
    H0(:, :, 1) = (hh11 - hh12)/(2*epsilon1);
    hh21 = MomentComponents([b(1); b(2) + epsilon2 ]);
    hh22 = MomentComponents([b(1); b(2) - epsilon2 ]);
    H0(:, :, 2) = ( hh21 - hh22 )/(2*epsilon2);
    H(1,:) = p' * H0(:,:,1);
    H(2,:) = p' * H0(:,:,2);
    % H = permute(H, [2 3 1]);
    % H = H';
    
    
    
    
%     sd_g = sd_g';
%     b1 = b(1); b2 = b(2);
%     
%     y1 = y(:,1);
%     [lambda, dLambda ] = varyingFixedEffect(b1, m);
%     
%     if ~isempty(M_hat)
%         lambda = lambda(M_hat);
%         dLambda = dLambda(M_hat);
%     end
%     
%     Dg_1 =  repmat(y1 - b2, [1 mm] ).* repmat(dLambda, [n 1]);
%     Dg_2 =  repmat( lambda - 1, [n 1] );
%     
%     term_1 = zeros(dx * n, mm);
%     term_1( dx*(1:n)-1,  :) =   Dg_1./ repmat(sd_g, [n 1]); % derivative wrt b1
%     term_1( dx*(1:n)  ,  :) =   Dg_2./ repmat(sd_g, [n 1]); % derivative wrt b2
%     
%     COV = zeros(1,mm);
%     for j = 1:mm
%         COV(j) = mean( gg(:,j).* Dg_1(:,j) ) - mean(gg(:,j)) * mean(Dg_1(:,j));
%     end
%     
%     term_2 = zeros(dx * n, mm);
%     term_2( dx*(1:n) -1, :) = 0.5 *  gg .* repmat( COV./ (sd_g.^(-3)), [n 1] ); % derivative wrt b1
%     
%     COV = zeros(1,mm);
%     for j = 1:mm
%         COV(j) = mean( gg(:,j).* Dg_2(:,j) ) - mean(gg(:,j)) * mean(Dg_2(:,j));
%     end
%     term_2( dx*(1:n) , :) = 0.5 *  gg .* repmat( COV./ (sd_g.^(-3)), [n 1] );
%     
%     term1_cell = mat2cell(term_1, dx * ones(n, 1), mm );
%     term2_cell = mat2cell(term_2, dx * ones(n, 1), mm );
%     
%     H0 = cellfun(@(p, term1, term2)   p * ( term1 -  term2 ) ,...
%         p_cell, term1_cell, term2_cell, 'UniformOutput', false );
%     H1 = reshape( cell2mat(H0), [ dx n mm]);
%     H2 = permute(H1, [ 1 3 2]);
%     H  = sum(H2, 3); % NEGTIVE? think again   
end




    %%
    hh_cell = mat2cell(hh, ones(n, 1), mm);
    VV_cell = cellfun(@(p, hh) p * hh' * hh,  p_cell, hh_cell, 'UniformOutput', false);
    VV1 = reshape( cell2mat(VV_cell), [mm n mm]);
    VV2 = permute(VV1, [1 3 2]);
    V = sum(VV2, 3);
end

