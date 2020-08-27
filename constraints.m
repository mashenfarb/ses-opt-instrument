
function [c, ceq] = constraints(x, cx1, cx2, ci1, ci2, b)

%%% Name control variables %%%
if length(x) == 2
    X1 = x(1);
    X2 = x(1);
    I1 = x(2);
    I2 = x(2);
elseif length(x) == 4
    X1 = x(1);
    X2 = x(2);
    I1 = x(3);
    I2 = x(4); 
end 

%%% Budget constraint for total cost of interventions %%%
budget_c = (cx1*X1^2 + cx2*X2^2 + ci1*I1^2 + ci2*I2^2) - b;     %% Budget constraint

c = budget_c;
ceq = [];
