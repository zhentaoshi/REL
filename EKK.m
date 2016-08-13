function [entry_pred, X_bar_pred] = EKK(b)

global m Sn
[ delta_nF, X_nF, u_bar ] = EKK_simulator(Sn, m + 1, b);

N_nF_art  = sum(delta_nF);
entry_pred = mean(delta_nF);
entry_pred(1) = []; % remove the home country

X_bar_pred = sum(X_nF)./ N_nF_art;
X_bar_pred(1) = [] ; % remove the home country

%% notation
% theT = theta_tilde
% lamb   = lambda
% rho  = rho
% sigA = sigma_alpha
% sigE = sigma_eta

% k1 = kappa_1
% k2 = kappa_2

    function [ deltanf, xnF, ubar ] = EKK_simulator(S, mm, b)
        % Eaton, Kortum and Kramarz (2011). 7 steps.
        thetatil = para(1);
        lambda = para(2);
        sigma_alpha = para(3);
        sigma_eta = para(4);
        rho = para(5);

        global init Ran S dom

        a    = Ran.a;
        h    = Ran.h;
        v    = Ran.v;

        X_bar = init.X_bar;
        N_nF  = init.N_nF;
        %% step 0
        % in the paper, S is number of firms. s is the firm index. n is the index of countrys
        % in this script, I replace S by n, replace n by m

        %step1 (34) and (35) in EKK
        kappa1 = (thetatil/(thetatil-1) - thetatil/(thetatil+lambda-1)) * exp((sigma_alpha^2 + 2*sigma_alpha*sigma_eta*rho*(thetatil-1) + sigma_eta^2*(thetatil-1)^2)/2);
        kappa2 = exp((thetatil*sigma_eta)^2/2);

        %step2 (28)
        sigmaEnf = kappa2/kappa1 * Xbarnf;    %Xbarnf 1*S, 1*S

        %step3
        ln_alpha = sigma_alpha*sqrt(1-rho^2)*a + sigma_alpha*rho*h;
        ln_eta = sigma_eta*h;
        alpha = exp(ln_alpha);  %n*S
        eta = exp(ln_eta);  %n*S

        %step4

        Nnf = repmat(Nnf,[n 1]);    %Nnf is 1*S from data, then transform to n*S
        ubarn = Nnf ./ kappa2 .* (eta .^thetatil);  %n*S

        %step5
        ubarn_temp = ubarn;
        ubarn_temp(:,Sdom) = [];  %n*(S-1)
        ubarX = (max(ubarn_temp'))';  %n*1
        ubar = (min([ubarn(:,Sdom) ubarX]'))';  %n*1

        %step6
        u = v.*ubar; % n*1

        %step 7

        %step 8
        deltanf = (repmat(u,[1 S+1])<=ubarn); % 0-1 matrix n*S
        utemp = repmat(u,[1 S+1]) ./ ubarn;    %n*S
        sigmaEnftemp = repmat(sigmaEnf,[n 1]);  %n*S
        xnf = (alpha ./ eta) .* (1 - (utemp).^(lambda/thetatil)) .* (utemp).^(-1/thetatil) .* sigmaEnftemp;
        xnf = xnf .* deltanf;

    end

end
