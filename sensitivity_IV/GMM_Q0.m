function [Q] = GMM_Q0( beta1 )

global n DGP
% the interpretatino is still the same. 

[~, gg, ~] = MomentComponents(beta1);

        % m0 = size(gg, 2);

        if strcmp(DGP, 'EKK')    
            gg_mean = mean(gg);
            naIndex = isnan(gg_mean);
            gg = gg(:, ~naIndex );
            % mm = size(gg, 2);
        else
            % mm = m0;
        end
        
        g_mean_no_na = mean(gg);
        Q = n * (g_mean_no_na * g_mean_no_na');
        
end


