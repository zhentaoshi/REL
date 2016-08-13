function [Q] = GMM_Q( beta1 )

global n DGP W_GMM
% the interpretatino is still the same. 

[~, gg, ~] = MomentComponents(beta1);

        % m0 = size(gg, 2);

        if strcmp(DGP, 'EKK')    
            gg_mean = mean(gg);
            naIndex = isnan(gg_mean);
            gg = gg(:, ~naIndex );
            W_no_na = W_GMM( ~naIndex, ~naIndex );
            if any(any(isnan(W_no_na)) )|| any(any(isinf(W_no_na)))
                V = eye(size(W_no_na));
            else
                V = pinv(W_no_na);
            end
            % mm = size(gg, 2);
        else
            V = pinv(W_GMM);
            % mm = size(gg, 2);
        end
        
        
        g_mean_no_na = mean(gg);
        Q = n * (g_mean_no_na * V * g_mean_no_na');
        
end


