function [h, g, sd_g] = MomentComponents(b)

% gg is the raw moments evaluated at the give b
% sd_g is the standard deviation of gg
% hh = gg / sd_g

global x y z DGP m M_hat 
% global boot gg_adj

    n = size(y,1);
    
    if strcmp(DGP, 'linearIV')
        e =  y - x * b ;
        g =  bsxfun(@times, z,  e); % a row vector
    end
    
    %%

    if strcmp(DGP, 'Han')
        b1 = b(1);
        b2 = b(2);
        y1 = y(:,1); % DO NOT CHANGE Y. Y is a global variable.
        if isempty(M_hat)
            y2 = y(:,2:(m+1));
            lambda = varyingFixedEffect(b1, m);
            g = bsxfun(@minus, ...
                y2 - repmat(lambda, [n 1]) .* repmat(y1, [1 m]), ...
                (1-lambda)*b2 );
        else
            y2 = y(:, M_hat + 1);
            mm = length(M_hat);
            lambda = varyingFixedEffect(b1, m);
            lambda = lambda(M_hat);
            g = bsxfun(@minus, ...
                y2 - repmat(lambda, [n 1]) .* repmat(y1, [1 mm]), ...
                (1-lambda)*b2 );
        end
    end
    %%

    if strcmp(DGP, 'EKK')

        yArt = EKK(b); % this is the theoretical prediction
        yArt_moment = repmat( sum(yArt), [n 1]);
        g = y - yArt_moment; % true data is a matrix. Art data gives a moment.
        % the probability of course is on the true data points.
    end

%%
    sd_g = std(g, 1);
    h = g./repmat(sd_g, [n, 1]);
    sd_g = sd_g';
end