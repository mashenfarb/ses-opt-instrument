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
r1 = .5;  % Initial rate of introduction, patch 1
r2 = .5;  % Initial rate of introduction, patch 2
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
[opt_XI_homog,fval_homog,eflag,output] = fmincon(@(x)objective_homogControls(x, d1, d2, r1, r2, p, q, s1, s2),...
                                                 x_0, [],[],[],[], x_lb, x_ub,...
                                                 @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
                                                 options);

%%% Print results %%%
sprintf('X = %.4f', opt_XI_homog(1))
sprintf('I = %.4f', opt_XI_homog(2))
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
[opt_XI,fval,eflag,output] = fmincon(@(x)objective(x, d1, d2, r1, r2, p, q, s1, s2),...
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
s1 = .5;
s2 = .5;
b = 1;

%%%%  Set Axes Limits  %%%%
p_lo = .00001;
p_hi = .99999;

q_lo = .00001;
q_hi = .99999;

%%%% Set Axes Intervals %%%%
n_p = 50;
n_q = 50;

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
        [opt_XI_homog,fval_homog,eflag,output] = fmincon(@(x)objective_homogControls(x, d1, d2, r1, r2, p_i, q_j, s1, s2),...
                                             x_0, [],[],[],[], x_lb, x_ub,...
                                             @(x)constraints(x, cx1, cx2, ci1, ci2, b),...
                                             options);

%         %%%  Patch-specific Controls Case with Homog as Initial Guess  %%%
%         x_0 = [opt_XI_homog(1); opt_XI_homog(1); opt_XI_homog(2); opt_XI_homog(2)];
%         x_lb = [0;0;0;0];  % lower bound for control variables
%         x_ub = [1,1,1,1];  % upper bound for control variables
%         [opt_XI,fval,eflag,output] = fmincon(@(x)objective(x, d1, d2, r1, r2, p_i, q_j, s1, s2),...
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
title(sprintf('Share of Budget to I, a1 = %g, a2 = %g', s1, s2))
ylabel('p - patch to patch infection rate')
xlabel('q - private land intervention spread rate')
colormap default
colorbar
% caxis([0,1])



%%%%%%%%%%
%%   Figures 2 and 3 -
%%   (2) Line Graph - Homogenous actors, Share of Budget to I as Budget Increases
%%   (3) Line Graph - Homogenous actors, Objective Values as Budget Increases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set Parameters Space for Budget, P, Q, and Efficacy %%%%
b_i_vec = linspace(0, 2, 101);
p_j_vec = [.99; .8; .5; .2; .01];
q_j_vec = [.99; .8; .5; .2; .01];
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
        [opt_XI_hom,fval,eflag,output] = fmincon(@(x)objective_homogControls(x, d1, d2, r1, r2, p_j, q_j, s1, s2),...
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
%         [opt_XI,fval,eflag,output] = fmincon(@(x)objective(x, d1, d2, r1, r2, p_j, q_j, s1, s2),...
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
xticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'})
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
xticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'})
ylabel('Expected Damages')
legend(legend_vec)


%%%%%%%%%%
%%   Figures 6 and 7 -  X and I as a function of budget, vary efficacy, homogenous patches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set Parameters Space for Budget, P, Q, and Efficacy %%%%
b_i_vec = linspace(0, 2, 101);
s_j_vec = [1; .8; .5; .2; 0];
% s1_j_vec = s_j_vec;
% s2_j_vec = s_j_vec;
s1_j_vec = [1; .8; .5; .2; 0];
s2_j_vec = [1; .8; .5; .2; 0];

%%%%  Run simulation for parameter spaces  %%%%
options = optimset('TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 100000, 'MaxIter', 100000, 'Disp', 'None');
out_I_budgetshare = ones(length(s1_j_vec), length(b_i_vec));
out_X_budgetshare = ones(length(s1_j_vec), length(b_i_vec));
out_obj_value = ones(length(s1_j_vec), length(b_i_vec));
for j = 1:length(s_j_vec)
    s1_j = s1_j_vec(j);
    s2_j = s2_j_vec(j);
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
        [opt_XI_hom,fval_hom,eflag_hom,output_hom] = fmincon(@(x)objective_homogControls(x, d1, d2, r1, r2, p, q, s1_j, s2_j),...
                                             x_0, [],[],[],[], x_lb, x_ub,...
                                             @(x)constraints(x, cx1, cx2, ci1, ci2, b_i),...
                                             options);
        
        %%%  Final optimal X and I for patch-specific controls %%%
        if b_i == b_i_vec(1)  % Initial guess after first run is previous results
            x_0 = [opt_XI_hom(1); opt_XI_hom(1); opt_XI_hom(2); opt_XI_hom(2)];
        else
            x_0 = opt_XI;
        end
        x_lb = [0;0;0;0];  % lower bound for control variables
        x_ub = [1,1,1,1];  % upper bound for control variables
        [opt_XI,fval,eflag,output] = fmincon(@(x)objective(x, d1, d2, r1, r2, p, q, s1_j, s2_j),...
                                     x_0, [],[],[],[], x_lb, x_ub,...
                                     @(x)constraints(x, cx1, cx2, ci1, ci2, b_i),...
                                     options);
        
        %%% Print and Save results %%%
        if s1_j == s2_j
            X1_opt = opt_XI_hom(1);
            X2_opt = opt_XI_hom(1);
            I1_opt = opt_XI_hom(2);
            I2_opt = opt_XI_hom(2);
            out_obj_value(j,i) = fval_hom;
        else
            X1_opt = opt_XI(1);
            X2_opt = opt_XI(2);
            I1_opt = opt_XI(3);
            I2_opt = opt_XI(4);
            out_obj_value(j,i) = fval;
        end
        total_I = ci1*I1_opt^2 + ci2*I2_opt^2;  %% total I costs (hard-coding costs from function)
        total_X = cx1*X1_opt^2 + cx2*X2_opt^2;  %% X costs (hard-coded cost function)
        out_I_budgetshare(j,i) = total_I/b_i;   %% save proportion of budget on I
        out_X_budgetshare(j,i) = total_X/b_i;   %% proportion of budget to X
        
        disp(opt_XI)
        disp([i;j])
        disp(eflag)
        disp(fval)
    end
end
%%  Figure 6 - Line Graph - Shares of X and I as function of budget, varying efficacy, homog actors
legend_vec = string(zeros(length(s1_j_vec)*2, 1));
color_vec = [linspace(.1,1,length(s1_j_vec)); ...
             linspace(.1,.1,length(s1_j_vec)); ...
             linspace(1,.1,length(s1_j_vec))];
for j = 1:length(s1_j_vec)
    j_legend = j*2-1;
    plot(b_i_vec, out_I_budgetshare(j,:), 'LineStyle', '-', 'Color', color_vec(:,j))
    hold on
    legend_vec(j_legend,1) = sprintf('I Budget Share, s1=%g, s2=%g', s1_j_vec(j), s2_j_vec(j));
    
    plot(b_i_vec, out_X_budgetshare(j,:), 'LineStyle', '-.', 'Color', color_vec(:,j))
    hold on
    legend_vec(j_legend+1,1) = sprintf('X Budget Share, s1=%g, s2=%g', s1_j_vec(j), s2_j_vec(j));
end
hold off
title('Effect of Efficacy on Optimal Intervention Mix as Function of Budget')
xlabel('Budget')
xticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'})
ylabel('Share of Budget')
legend(legend_vec)

%%  Figure 7 - Line Graph - Expected damages as a function of budget, varying efficacy, homog actors
legend_vec = string(zeros(length(s_j_vec), 1));
color_vec = [linspace(.1,1,length(s_j_vec)); ...
             linspace(.1,.1,length(s_j_vec)); ...
             linspace(1,.1,length(s_j_vec))];
for j = 1:length(s_j_vec)
    plot(b_i_vec, out_obj_value(j,:), 'LineStyle', '-', 'Color', color_vec(:,j))
    hold on
    legend_vec(j,1) = sprintf('Expected Damages, s1=%g s2=%g', s1_j_vec(j), s2_j_vec(j));
end
hold off
title('Effect of Efficacy on Expected Damages as Function of Budget')
xlabel('Budget')
xticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'})
ylabel('Expected Damges')
legend(legend_vec)



%%%%%%%%%%
%%   Figures 8 -  Exploring patch heterogeneity and efficacy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set Parameters Space for Budget, P, Q, and Efficacy %%%%
b_i_vec = linspace(0, 2, 101);
s1_j_vec = [.8; .2; .8; .5];
s2_j_vec = [.8; .2; .2; .2];

q = .5;
p = .5;

%%%%  Run simulation for parameter spaces  %%%%
options = optimset('TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 100000, 'MaxIter', 100000, 'Disp', 'None');
out_I_budgetshare = ones(length(s1_j_vec), length(b_i_vec));
out_X_budgetshare = ones(length(s1_j_vec), length(b_i_vec));
out_obj_value = ones(length(s1_j_vec), length(b_i_vec));
out_X1_vec = ones(length(s1_j_vec), length(b_i_vec));
out_X2_vec = ones(length(s1_j_vec), length(b_i_vec));
out_I1_vec = ones(length(s1_j_vec), length(b_i_vec));
out_I2_vec = ones(length(s1_j_vec), length(b_i_vec));
for j = 1:length(s1_j_vec)
    s1_j = s1_j_vec(j);
    s2_j = s2_j_vec(j);
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
        [opt_XI_hom,fval_hom,eflag_hom,output_hom] = fmincon(@(x)objective_homogControls(x, d1, d2, r1, r2, p, q, s1_j, s2_j),...
                                             x_0, [],[],[],[], x_lb, x_ub,...
                                             @(x)constraints(x, cx1, cx2, ci1, ci2, b_i),...
                                             options);
        
        %%%  Final optimal X and I for patch-specific controls %%%
        if b_i == b_i_vec(1)  % Initial guess after first run is previous $results
            x_0 = [opt_XI_hom(1); opt_XI_hom(1); opt_XI_hom(2); opt_XI_hom(2)];
        else
            x_0 = opt_XI;
        end
        x_lb = [0;0;0;0];  % lower bound for control variables
        x_ub = [1,1,1,1];  % upper bound for control variables
        [opt_XI,fval,eflag,output] = fmincon(@(x)objective(x, d1, d2, r1, r2, p, q, s1_j, s2_j),...
                                     x_0, [],[],[],[], x_lb, x_ub,...
                                     @(x)constraints(x, cx1, cx2, ci1, ci2, b_i),...
                                     options);
        
        %%% Print and Save results %%%
        if s1_j == s2_j
            X1_opt = opt_XI_hom(1);
            X2_opt = opt_XI_hom(1);
            I1_opt = opt_XI_hom(2);
            I2_opt = opt_XI_hom(2);
            out_obj_value(j,i) = fval_hom;
        else
            X1_opt = opt_XI(1);
            X2_opt = opt_XI(2);
            I1_opt = opt_XI(3);
            I2_opt = opt_XI(4);
            out_obj_value(j,i) = fval;
        end
        total_I = ci1*I1_opt^2 + ci2*I2_opt^2;  %% total I costs (hard-coding costs from function)
        total_X = cx1*X1_opt^2 + cx2*X2_opt^2;  %% X costs (hard-coded cost function)
        out_I_budgetshare(j,i) = total_I/b_i;   %% save proportion of budget on I
        out_X_budgetshare(j,i) = total_X/b_i;   %% proportion of budget to X
        out_X1_vec(j,i) = X1_opt;
        out_X2_vec(j,i) = X2_opt;
        out_I1_vec(j,i) = I1_opt;
        out_I2_vec(j,i) = I2_opt;
        
        disp(opt_XI)
        disp([i;j])
        disp(eflag)
        disp(fval)
    end
end
%%  Figure 8 - Line Graph - Plot X1, X2, I1, I2. Four separate graphs for efficacy combinations

legend_vec = string(zeros(length(s_j_vec), 1));
color_vec = [linspace(.1,1,length(s_j_vec)); ...
             linspace(.1,.1,length(s_j_vec)); ...
             linspace(1,.1,length(s_j_vec))];
figure
for j = 1:length(s1_j_vec)
    ax = subplot(2,2,j);
    plot(b_i_vec, out_X1_vec(j,:), 'LineStyle', '-', 'Color', 'r') % color_vec(:,j))
    hold on
    plot(b_i_vec, out_X2_vec(j,:), 'LineStyle', ':', 'Linewidth', 2, 'Color', 'b') % color_vec(:,j))
    hold on
    plot(b_i_vec, out_I1_vec(j,:), 'LineStyle', '-', 'Color', 'g') % color_vec(:,j))
    hold on
    plot(b_i_vec, out_I2_vec(j,:), 'LineStyle', ':', 'Linewidth', 2, 'Color', 'k') % color_vec(:,j))
    hold off
    title(ax, sprintf('a1=%g a2=%g', s1_j_vec(j), s2_j_vec(j)))
    legend(['X1'; 'X2'; 'I1'; 'I2'])
    xlabel('Budget')
    xticklabels({'0', '.25', '0.5', '0.75', '1.0'})
    ylabel('Intervention Level')
end
