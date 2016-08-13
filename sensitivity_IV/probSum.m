function [totalProb, p] = probSum( b, gam )
global n tau

[hh, gg, sd_g] = MomentComponents(b);
p = 1./( n * ( 1 +  gg * (gam./sd_g) - tau * norm(gam,1) ) );
totalProb = sum(p);

end

