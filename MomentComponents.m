function [h, g, sd_g] = MomentComponents(b)

% gg is the raw moments evaluated at the give b
% sd_g is the standard deviation of gg
% hh = gg / sd_g

global y 

n = size(y,1);
 % only match the mean sell
    [~, yArt_sell] = EKK_reduced2(b); % this is the theoretical prediction
    g = bsxfun(@minus, y, yArt_sell ); % true data is a matrix. Art data gives a moment.

%%
sd_g = std(g, 1);
% g( abs( mean(g) ) > 1e+5, : ) = NaN;

h = g./repmat(sd_g, [n, 1]);
sd_g = sd_g';
end