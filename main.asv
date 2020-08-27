clear all
close all

%%%% Set working directory to main folder %%%%
[dir,fn,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(dir)

%%%%%%%%%%%%%%%%%%%
%%   Parameters  %%
%%%%%%%%%%%%%%%%%%%
%%%% Economic %%%%
cx1 = 1;  % Public intervention, patch 1 [0, inf)
cx2 = 1;  % Public intervention, patch 2 [0, inf)
ci1 = 1;  % Information intervention, patch 1 [0, inf)
ci2 = 1;  % Information intervention, patch 2 [0, inf)

d1 = 1;  % Nominal damages incurred per period in patch 1
d2 = 1;  % Nominal damages incurred per period in patch 2

b = 1;  % Budget for intervention

dis = .05;  % discount rate (not used)

%%%% Biological %%%%
r1 = .99;  % Initial rate of introduction, patch 1
r2 = .99;  % Initial rate of introduction, patch 2
p = .5;  % probability of patch-to-patch spread without control

%%%% Social %%%%
q = .2;  % Information spread rate
q1 = .5;  % Patch 1 to patch 2 influence
q2 = .5;  % Patch 2 to patch 1 influence

s1 = 1;  % Informational intervention efficacy, patch 1
s2 = 1;  % Informational intervention efficacy, patch 2



%%%%%%%%%
%%   Optimization Test - Homogenous Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set initial guess, bounds, and options %%%
x_0 = [0;1];  % [b/2;b/2];  % x = [X; I]
x_lb = [0;0];  % lower bound for control variables
x_ub = [1,1];  % upper bound for control variables

%%% Find optimal X and I with homogenous patch controls %%%
options = optimset('TolCon', 1e-8, 'TolFun', 1e-8, 'TolX', 1e-8, 'MaxFunEvals', 100000, 'MaxIter', 100000);
[opt_XI_homog,fval_homog,eflag,output] = fmincon(@(x)objective5_homogControls(x, d1, d2, r1, r2, p, q, s1, s2),...
                                                 x_0, [],[],[],[], x_lb, x_ub,...
                                                 @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
                                                 options);

%%% Print results %%%
sprintf('X = %.4f', opt_XI(1))
sprintf('I = %.4f', opt_XI(2))
sprintf('Expected Sum of Patch Damages = %.3f', fval_homog)
sprintf('Fmincon Exitflag: %i', eflag)


%%%%%%%%%
%%   Optimization Test - Patch-specific Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set initial guess, bounds, and options %%%
x_0 = [0; 0; 1; 1];  % x = [X1, X2, I1, I2]
x_lb = [0;0;0;0];    % lower bound for control variables
x_ub = [1,1,1,1];    % upper bound for control variables

%%% fmincon optimization call %%%
options = optimset('TolCon', 1e-12, 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 100000, 'MaxIter', 100000);
[opt_XI,fval,eflag,output] = fmincon(@(x)objective5(x, d1, d2, r1, r2, p, q, s1, s2),...
                                     x_0, [],[],[],[], x_lb, x_ub,...
                                     @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
                                     options);     

%%% Print results %%%
sprintf('X1 = %.4f', opt_XI(1))
sprintf('X2 = %.4f', opt_XI(2))
sprintf('I1 = %.4f', opt_XI(3))
sprintf('I2 = %.4f', opt_XI(4))
sprintf('Expected Sum of Patch Damages = %.3f', fval)
sprintf('Fmincon Exitflag: %i', eflag)


%%%%%%%%%%%
%%   Figure - Homogenous actors, QvsP Contour Plot, Percent of Budget to I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Set efficacy and budget parameters  %%%%
s1 = 1;
s2 = 1;
b = 1;

%%%%  Set Axes Limits  %%%%
p_lo = .00001;
p_hi = .99999;

q_lo = .00001;
q_hi = .99999;

%%%% Set Axes Intervals %%%%
n_p = 30;
n_q = 30;

%%%%  Run simulation, Iterating over parameter space (n_p x n_q)  %%%%
p_i_vec = linspace(p_lo, p_hi, n_p);  %% values of p to calculate (ecological connectivity)
q_j_vec = linspace(q_lo, q_hi, n_q);  %% values of q to calculate (social connectivity)
out_I_budgetshare = ones(length(p_i_vec),length(q_j_vec));
options = optimset('TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 100000, 'MaxIter', 100000, 'Display', 'None');
for i = 1:length(p_i_vec)
    for j = 1:length(q_j_vec)
        p_i = p_i_vec(i);
        q_j = q_j_vec(j);

        %%% Homogenous control case %%%
        x_lb = [0;0];  % lower bound for control variables
        x_ub = [1;1];  % upper bound for control variables
        x_0 = [.5; .5];
        [opt_XI_homog,fval_homog,eflag,output] = fmincon(@(x)objective5_homogControls(x, d1, d2, r1, r2, p_i, q_j, s1, s2),...
                                             x_0, [],[],[],[], x_lb, x_ub,...
                                             @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
                                             options);

%         %%%  Patch-specific Controls Case with Homog as Initial Guess  %%%
%         x_0 = [opt_XI_homog(1); opt_XI_homog(1); opt_XI_homog(2); opt_XI_homog(2)];
%         x_lb = [0;0;0;0];  % lower bound for control variables
%         x_ub = [1,1,1,1];  % upper bound for control variables
%         [opt_XI,fval,eflag,output] = fmincon(@(x)objective5(x, d1, d2, r1, r2, p_i, q_j, s1, s2),...
%                                              x_0, [],[],[],[], x_lb, x_ub,...
%                                              @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
%                                              options);
        %%% Save and Print Print results %%%
        total_I = ci1*opt_XI_homog(2)^2 + ci2*opt_XI_homog(2)^2;  %% total I costs (hard-coding costs from function)
        out_I_budgetshare(i,j) = total_I/b;     %% save proportion of budget on I
        disp(opt_XI_homog)
        disp([i;j])
        disp(eflag)
        disp(fval)
    end
end
%%  Make Plot %%
% s = pcolor(q_lo:((q_hi-q_lo)/(n_q - 1)):q_hi, p_lo:((p_hi-p_lo)/(n_p - 1)):p_hi, out_I_budgetshare);  % X = 1:n, Y = 1:m, where [m,n] = size(C)
% set(s,'Edgecolor','None')
s = contour(q_lo:((q_hi-q_lo)/(n_q - 1)):q_hi, p_lo:((p_hi-p_lo)/(n_p - 1)):p_hi, out_I_budgetshare, 15);
set(gca,'layer','top');
title({'Share of Budget to I Depends on Invasion Spread and Information Spread'})
ylabel('p - patch to patch infection rate')
xlabel('q - private land intervention spread rate')
colormap default
colorbar



%%%%%%%%%%
%%   Figures 2 and 3 -
%%   (2) Line Graph - Homogenous actors, Share of Budget to I as Budget Increases
%%   (3) Line Graph - Homogenous actors, Objective Values as Budget Increases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set Parameters Space for Budget, P, Q, and Efficacy %%%%
b_i_vec = linspace(0, 2, 101);
p_j_vec = [.1; .5; .6; .7; .8; .9;];
q_j_vec = [.1; .5; .6; .7; .8; .9;];
% p_j_vec = [.8;.5;.2;.8];
% q_j_vec = [.2;.5;.8;.8];

s1 = 1;  % information efficacy patch 1
s2 = 1;  % information efficacy patch 2

%%%%  Run simulation for parameter spaces  %%%%
options = optimset('TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 100000, 'MaxIter', 100000, 'Disp', 'None');
out_I_budgetshare = ones(length(p_j_vec), length(b_i_vec));
out_obj_value = ones(length(p_j_vec), length(b_i_vec));
for j = 1:length(p_j_vec)
    p_j = p_j_vec(j);
    q_j = q_j_vec(j);
    for i = 1:length(b_i_vec)
        b_i = b_i_vec(i);
        
        %%%  Set initial guess and bounds  %%%
        if b_i == b_i_vec(1)  % Initial guess after first run is previous results
            x_0 = [0;1];
        else
            x_0 = opt_XI_hom;
        end
        x_lb = [0;0];  % lower bound for control variables
        x_ub = [1;1];  % upper bound for control variables
        
        %%%  Find optimal X and I for homogenous patch controls  %%%
        [opt_XI_hom,fval,eflag,output] = fmincon(@(x)objective5_homogControls(x, d1, d2, r1, r2, p_j, q_j, s1, s2),...
                                             x_0, [],[],[],[], x_lb, x_ub,...
                                             @(x)constraints(x, cx1, cx2, ci1, ci2, b_i),...
                                             options);
         
        %%%  Final optimal X and I for patch-specific controls (only using homog for now)  %%%
%         if  b_i < 1.3 || b_i > 1.3
%             x_0 = [opt_XI_hom(1); opt_XI_hom(1); opt_XI_hom(2); opt_XI_hom(2)];
%         else
%             %%% Use previous sim as initial guess for high values %%%
%             x_0 = opt_XI;
%         end
%         x_lb = [0;0;0;0];  % lower bound for control variables
%         x_ub = [1,1,1,1];  % upper bound for control variables
%         [opt_XI,fval,eflag,output] = fmincon(@(x)objective5(x, d1, d2, r1, r2, p_j, q_j, s1, s2),...
%                                      x_0, [],[],[],[], x_lb, x_ub,...
%                                      @(x)constraints(x, cx1, cx2, ci1, ci2, b_i),...
%                                      options);
        %%% Print and Save results %%%
        X1_opt = opt_XI_hom(1);
        X2_opt = opt_XI_hom(1);
        I1_opt = opt_XI_hom(2);
        I2_opt = opt_XI_hom(2);
        total_I = ci1*I1_opt^2 + ci2*I2_opt^2;  %% total I costs (hard-coding costs from function)
        out_I_budgetshare(j,i) = total_I/b_i;   %% save proportion of budget on I
        out_obj_value(j,i) = fval;
        disp(opt_XI)
        disp([i;j])
        disp(eflag)
        disp(fval)
    end
end
%%  Figure 2 - Line Graph - Homogenous actors, Share of Budget to I as Budget Increases
legend_vec = string(zeros(length(p_j_vec), 1));
for j = 1:length(p_j_vec)
    plot(b_i_vec, out_I_budgetshare(j,:))
    hold on
    legend_vec(j) = sprintf('p=%g, q=%g', p_j_vec(j), q_j_vec(j));
end
hold off
title('Share of Optimal Intervention to I Depends on Budget, Invasion Spread, and Information Spread')
xlabel('Budget')
ylabel('Percent of Budget to Information')
legend(legend_vec)

%%  Figure 3 - Line Graph - Homogenous actors, Objective Values as Budget Increases
legend_vec = string(zeros(length(p_j_vec), 1));
for j = 1:length(p_j_vec)
    plot(b_i_vec, out_obj_value(j,:))
    hold on
    legend_vec(j) = sprintf('p=%g, q=%g', p_j_vec(j), q_j_vec(j));
end
hold off
title('Expected Damages and Budget')
xlabel('Budget')
ylabel('Expected Damages')
legend(legend_vec)




%%%%%%%%%
%%   Figure 4 - Test different efficacy parameters (work in progress)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = 1;
s_i_vec = linspace(0, 1, 100);
p_j_vec = [.1; .5; .9; .1; .9;];
q_j_vec = [.1; .5; .9; .9; .1;];

options = optimset('TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 100000, 'MaxIter', 100000, 'Disp', 'None');

out_I_budgetshare = ones(length(p_j_vec), length(s_i_vec));
out_obj_value = ones(length(p_j_vec), length(s_i_vec));

for j = 1:length(p_j_vec)
    p_j = p_j_vec(j);
    q_j = q_j_vec(j);
    for i = 1:length(s_i_vec)
        s1_i = s_i_vec(i);
        s2_i = s_i_vec(i);
        
        %%% Run homogenous controls for initial guess %%%
        x_0 = [0;1];
        x_lb = [0;0];  % lower bound for control variables
        x_ub = [1;1];  % upper bound for control variables
        [opt_XI_hom,fval,eflag,output] = fmincon(@(x)objective5_homogControls(x, d1, d2, r1, r2, p_j, q_j, s1_i, s2_i),...
                                             x_0, [],[],[],[], x_lb, x_ub,...
                                             @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
                                             options);

        %%% Print results %%%
        X1_opt = opt_XI_hom(1);
        X2_opt = opt_XI_hom(1);
        I1_opt = opt_XI_hom(2);
        I2_opt = opt_XI_hom(2);
        total_I = ci1*I1_opt^2 + ci2*I2_opt^2;
        out_I_budgetshare(j,i) = total_I/b;
        out_obj_value(j,i) = fval;
        disp(opt_XI_hom)
        disp([i;j])
        disp(eflag)
        disp(fval)
    end
end

%%  Figure 4 %%
legend_vec = string(zeros(length(p_j_vec), 1));
for j = 1:length(p_j_vec)
    plot(s_i_vec, out_I_budgetshare(j,:))
    hold on
    legend_vec(j) = sprintf('p=%g, q=%g', p_j_vec(j), q_j_vec(j));
end
hold off
title('Share of Optimal Intervention to I Depends on Information Efficacy, Invasion Spread, and Information Spread')
xlabel('Efficacy (s)')
ylabel('Percent of Budget to Information')
legend(legend_vec)

%%  Figure 5 %%
legend_vec = string(zeros(length(p_j_vec), 1));
for j = 1:length(p_j_vec)
    plot(s_i_vec, out_obj_value(j,:))
    hold on
    legend_vec(j) = sprintf('p=%g, q=%g', p_j_vec(j), q_j_vec(j));
end
hold off
title('Expected Damages and Efficacy')
xlabel('Efficacy')
ylabel('Expected Damages')
legend(legend_vec)
   