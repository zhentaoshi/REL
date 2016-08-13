%% boosting
function [M_hat] = boost0(M0, b)

global m y x z n
% M0: is the given set
% d: total number of structural parameters
% M: total number of moments
% this function add one more moment in each iteration

dx = size(x, 2);
M  = setdiff(1:m, M0); % M is the avaiable moments to explore
mm = length(M); % cardinality of the candidate set

% Q store the increment
Q = zeros(1,mm);

% in the loop each time we calculate the max eigenvalue as the increment.
S = 1/n * z' * x;

    e =  y - x * b' ;
    gg =  bsxfun(@times, z,  e); % a row vector
    sd_g = std(gg, 1);

    hh = gg./repmat(sd_g, [n, 1]);
    V = 1/n * ( hh' * hh );


% the selected coordinates

if isempty(M0)
    IC0 = zeros(dx, dx);
else
    S_M0 = S(M0, :);
    V_M0 = V(M0, M0);
    IC0 = S_M0' * pinv(V_M0) * S_M0;
end

for j = 1:mm
    Mj = [M0 M(j)]; % the candidate moments + the given moments from a previous run
    S_Mj = S(Mj, :);
    V_Mj = V(Mj, Mj);
    ICj = S_Mj' * pinv(V_Mj) * S_Mj;
  
	Q(j) = max( eig(ICj - IC0) );
end
[~, which] = max(Q); % find out the max increment
M_hat = M(which); % the index choose moment
end
