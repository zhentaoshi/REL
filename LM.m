function [stat, gam, gb] = LM(b)

global n tau boot

[L0, gam] = gamm_msk(b);
hh = MomentComponents(b);
m  = size(hh, 2);
pen = tau * norm(gam,1);
pi = 1./(n * ( 1 + hh * gam - pen) );

if L0 < 1e+10 % no solution
    gb = pi' * hh;
    stat = pen; % pi' * hh * gam;
else
    if boot == 1;
        stat = 0;
        gb   = zeros(1,m);
        
    elseif boot == 0;
        stat = Inf;
        gb   = zeros(1,m);
    end
    % 
end

end