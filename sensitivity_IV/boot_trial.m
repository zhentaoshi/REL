

n = 100;
T = 499;
R = 400;

YB = zeros(T);
decision = zeros(R, 1);
yb_bar = zeros(T, 1);
for r = 1:R
    r;
    for t = 1:T


        y = random('Normal', 0, 1, n , 1);
        % bsample = ceil( n * random('Uniform', 0, 1, n, 1 ) ); 
        bsample = randsample(1:n, n, 1);
        yb = y(bsample);
        yb_bar(t) = mean(yb) - mean(y);
        

    end
    decision(r) = ( mean( mean(y) > yb_bar ) > 0.50 );
end

sort(decision)
mean(decision)