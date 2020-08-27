function exp_damage = objective5(x, d1, d2, r1, r2, p, q, s1, s2)

% Changes:  Same as 4, but "race" is based on Poisson process

%%% Name control variables %%%
X1 = x(1);
X2 = x(2);
I1 = s1*x(3);
I2 = s2*x(4);

%%% Infection probabilities in t1 %%%
Pr_S1_t1 = r1*(1 - X1)*(1 - I1);  % probability patch 1 infected in t1
Pr_S2_t1 = r2*(1 - X2)*(1 - I2);  % probability patch 2 infected in t1

%%% Infection probabilities in t2, given conditions for t2 spread %%%
lambda_q = log(1/(1-q));
lambda_p = log(1/(1-p));

Pr_S1_t2 = (1 - X1)*(1 - I1) ...  % probability direct patch controls don't work
           *(1 - exp(-lambda_p - I2*lambda_q)) ...  % probability more than zero events (p and q) occur
           *(lambda_p/(I2*lambda_q + lambda_p));  %  probability first event is physical spread
Pr_S2_t2 = (1 - X2)*(1 - I2) ...  % probability direct patch controls don't work
           *(1 - exp(-lambda_p - I1*lambda_q)) ...  % probability more than zero events (p and q) occur
           *(lambda_p/(I1*lambda_q + lambda_p));  %  probability first event is physical spread

%%% Expected infection in t2 = infection in t1 OR t2 %%%
exp_S1 = Pr_S1_t1 + ... % prob t1 infection
         (1 - Pr_S1_t1)*Pr_S2_t1*Pr_S1_t2; % prob of conditions for t2 infection and t2 infection
exp_S2 = Pr_S2_t1 + ... % prob t1 infection
         (1 - Pr_S2_t1)*Pr_S1_t1*Pr_S2_t2; % prob of conditions for t2 infection and t2 infection

exp_damage = d1*exp_S1 + d2*exp_S2;
