function [M_hat] = boost(M0, S, V)

global m dx
% M0: is the given set
% dx: total number of structural parameters
% M: total number of moments
% this function add one more moment in each iteration

M  = setdiff(1:m, M0); % M is the avaiable moments to explore
mm = length(M); % cardinality of the candidate set

% Q store the increment
Q = zeros(dx,mm);

% the selected coordinates

if isempty(M0)
    IC0 = zeros(dx, dx);
else
    S_M0 = S(M0, :);
    V_M0 = V(M0, M0);
    IC0 = S_M0' * pinv(V_M0) * S_M0;
end

for j = 1:mm
    Mj = [reshape( M0, 1, [] ), M(j)]; % the candidate moments + the given moments from a previous run
    S_Mj = S(Mj, :);
    V_Mj = V(Mj, Mj);
    
    % minimal eigenvalue to prevent collinearity
    minEig = min( abs( eig(V_Mj) ) );
    if minEig < 0.01;
        ICj = IC0;
    else
        ICj = S_Mj' * pinv(V_Mj) * S_Mj;
    end
    
    if  any( any( isnan( ICj-IC0 ) ) )  || any( any( isinf(ICj-IC0) ) )
        Q(:,j) = 0;
    else
        % Q(j) = max( eig(ICj - IC0) ); %I used to use the max eigenvalue as the criterion.
        % Q(j) = trace( ICj - IC0);
        Q(:,j) = diag( ICj - IC0);
    end    
end
[~, which1] = max(Q(1,:)); % find out the max increment
[~, which2] = max(Q(2,:)); % find out the max increment
M_hat =  [M(which1);M(which2)] ; % the index choose moment
end