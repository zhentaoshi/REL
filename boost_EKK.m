function [M_hat] = boost_EKK(M0, S, V)

global  dx NotnanIndex dd n

% M0: is the given set
% d: total number of structural parameters
% M: total number of moments
% this function add one more moment in each iteration

% M  = setdiff(1:sum(NotnanIndex), M0); % M is the avaiable moments to explore
mm = size(V,1); % cardinality of the candidate set
M = 1:mm;

% Q store the increment
Q = zeros(1,mm);

% the selected coordinates

if isempty(M0)
    IC0 = 0;
else
    S_M0 = S(M0, dd);
    V_M0 = V(M0, M0);
    IC0 = S_M0' * pinv(V_M0) * S_M0;
end
  
for j = 1:mm
    if ismember( M(j) , M0 )
        Q(j) = 0;
    else
        Mj = [M0 M(j)]; % the candidate moments + the given moments from a previous run
        S_Mj = S(Mj, dd);
        V_Mj = V(Mj, Mj);
        minEig = min( abs( eig(V_Mj) ) );
        if minEig < 0.05/(log(log(n)));    %|| % tuning parameter
            ICj = IC0;
        else
            ICj = S_Mj' * pinv(V_Mj) * S_Mj;
        end
        
        if  isnan( ICj-IC0 )  || isinf(ICj-IC0)
            Q(j) = 0;
        else
            Q(j) = ICj-IC0;
        end
        
    end

end
[~, which] = max(Q); % find out the max increment
M_hat = M(which);% THIS IS IMPORTANT!
end
