function MinimumTimeBrachistochrone_LinearConstraint
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     MinimumTimeBrachistochrone_LinearConstraint.m
%    Compiler:      MATLAB R2022b
%    Date:          12 June, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve minimum-time Brachistochrone problem
%                   subject to a linear path constraint x + y >= -1.5
%    References:    Ch 3. Applied Optimal Control, 1975, A.E. Bryson. Jr, Yu-Chi Ho
%                   Antony, T., & Grant, M. J. (2018). Path constraint regularization in optimal control problems using saturation functions. 2018 AIAA Atmospheric Flight Mechanics Conference. https://doi.org/10.2514/6.2018-0018                 

clear; close all; clc;

% Parameters
t0 = 0;
x0 = [0; 0; 0]; % [x, y, v] (m and m/s)
xf = [10; -10]; % m
g = -9.81; % m/s^2

% Linear Constraint
xLine = linspace(0,10,1000);
yLine = -1.5-xLine;

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-14,'TolX',1e-16,...
    'UseParallel',true); % fsolve

%% Numerical Solution 
lam0_guess = [-0.0853743837971961; 0.163600488622047; -0.101936799182887; 0.0637722104987520; -0.0320927581757651; 0.0640095166701569; 0.926218792177670; 0.510173298116292; 7.24525873769311e-12; -0.149146594295948; 0.195693246797812];

iter = 1; err = 1; errTol = 1e-10; iterMax = 100;
    while err > errTol && iter < iterMax
        [lam0,~] = fsolve(@costFunctionConstrained,lam0_guess,options,x0,xf,opts_ode,g,t0);
        err = norm(costFunctionConstrained(lam0_guess,x0,xf,opts_ode,g,t0));
        lam0_guess = lam0;
        iter = iter+1; 
    end
 
% Trajectory
lam0_01 = lam0(1:3); lam0_12 = lam0(4:6);
t1 = lam0(7); t2 = lam0(8); tf = lam0(9);
[~,X_01] = ode89(@unconstrainedBrachistochrone,[t0 t1],[x0; lam0_01],opts_ode,g);
[~,X_12] = ode89(@constrainedBrachistochrone,[t1 t2],[X_01(end,1:3)'; lam0_12],opts_ode,g);
[~,X_2f] = ode89(@unconstrainedBrachistochrone,[t2 tf],X_12(end,1:6)',opts_ode,g);

% Combined Trajectories
X_minT = [X_01; X_12; X_2f];

% Plots
PlotSolution(X_minT,x0,xf,xLine,yLine); % plots

end


%% Functions
function PlotSolution(X_minT,x0,xf,xLine,yLine)

figure; grid on; hold on; xlim([0 10]);
plot(xLine,yLine,'--r','LineWidth',1,'DisplayName','Constraint (x + y >= -1.5)');
plot(X_minT(:,1),X_minT(:,2),'-k','LineWidth',1,'DisplayName','Trajectory');
xlabel('Distance (m)'); ylabel('Distance (m)'); 
title('Constrained Trajectory of Brachistochrone Problem',['x0 = ',num2str(x0(1:2)'),'  |  ','xf = ',num2str(xf(1:2)')]); 
legend('show','Location','best')

end 


function Xdot = unconstrainedBrachistochrone(t,X,g)

x = X(1:3);
lambda = X(4:6);

[theta,~] = hamiltonianUnconstrained(t,X',g);

% State Dynamics
xDot = [x(3)*cos(theta);...
        x(3)*sin(theta);...
        g*sin(theta)];

% Costate Dynamics
lambdaDot = [0;...
             0;...
             -lambda(1)*cos(theta) - lambda(2)*sin(theta)];

Xdot = [xDot; lambdaDot];

end


function Xdot = constrainedBrachistochrone(t,X,g)

x = X(1:3);
lambda = X(4:6);

[theta,~] = hamiltonianConstrained(t,X',g);

mu1 = -1/(2*x(3)) * (g*lambda(3) + x(3)*(lambda(1)+lambda(2)));

% State Dynamics
xDot = [x(3)*cos(theta);...
        x(3)*sin(theta);...
        g*sin(theta)];

% Costate Dynamics
lambdaDot = [0;...
             0;...
             -lambda(1)*cos(theta) - lambda(2)*sin(theta) - mu1*(sin(theta)+cos(theta))];

Xdot = [xDot; lambdaDot];

end


function [theta,Huc] = hamiltonianUnconstrained(~,X,g)

    for ii = 1:length(X)

        if size(X,1) == 1
            x = X(1:3);
            lambda = X(4:6);
        else
            x = X(ii,1:3);
            lambda = X(ii,4:6); 
        end

        theta = [-2*atan((lambda(1)*x(3) + sqrt(g^2*lambda(3)^2 + 2*g*lambda(3)*lambda(2)*x(3) + lambda(1)^2*x(3)^2 + lambda(2)^2*x(3)^2))/(g*lambda(3) + lambda(2)*x(3)));...
                 -2*atan((lambda(1)*x(3) - sqrt(g^2*lambda(3)^2 + 2*g*lambda(3)*lambda(2)*x(3) + lambda(1)^2*x(3)^2 + lambda(2)^2*x(3)^2))/(g*lambda(3) + lambda(2)*x(3)))];
        
        for jj = 1:length(theta)
        
            % State Dynamics
            xDot = [x(3)*cos(theta(jj));...
                    x(3)*sin(theta(jj));...
                    g*sin(theta(jj))];
            
            % Hamiltonian
            H(jj) = (lambda * xDot) + 1;

        end

        [Huc(ii),I] = min(H);
        theta = theta(I);
        
    end

end


function [theta,Hc] = hamiltonianConstrained(~,X,g)

    for ii = 1:length(X)

        theta = [pi/4 + pi/2; pi/4 - pi/2];

        if size(X,1) == 1
            x = X(1:3);
            lambda = X(4:6);
        else
            x = X(ii,1:3);
            lambda = X(ii,4:6); 
        end

        mu1 = -1/(2*x(3)) * (g*lambda(3) + x(3)*(lambda(1)+lambda(2)));
        
        for jj = 1:length(theta)
            
            % State Dynamics
            xDot = [x(3)*cos(theta(jj));...
                    x(3)*sin(theta(jj));...
                    g*sin(theta(jj))];
            
            % Hamiltonian
            H(jj) = (lambda * xDot) + 1 + mu1*(xDot(1)+xDot(2));

        end

        [Hc(ii),I] = min(H);
        theta = theta(I);
            
    end

end


function err = costFunctionConstrained(lam0_guess,x0,xf,opts_ode,g,t0)

% Control Variables
lam0_01 = lam0_guess(1:3); lam0_12 = lam0_guess(4:6);
t1 = lam0_guess(7); t2 = lam0_guess(8); tf = lam0_guess(9);
piVec = lam0_guess(10:11);

% Trajectories (unconstrained and constrained arcs)
[t_01,X_01] = ode89(@unconstrainedBrachistochrone,[t0 t1],[x0; lam0_01],opts_ode,g);
[t_12,X_12] = ode89(@constrainedBrachistochrone,[t1 t2],[X_01(end,1:3)'; lam0_12],opts_ode,g);
[t_2f,X_2f] = ode89(@unconstrainedBrachistochrone,[t2 tf],X_12(end,1:6)',opts_ode,g);

% Hamiltonian
[~,H_01] = hamiltonianUnconstrained(t_01,X_01,g);
[~,H_12] = hamiltonianConstrained(t_12,X_12,g);
[~,H_2f] = hamiltonianUnconstrained(t_2f,X_2f,g); 

err = [X_12(1,1)+X_12(1,2)+1.5, X_01(end,4:5)-X_12(1,4:5)-piVec', X_01(end,6)-X_12(1,6), H_01(end)-H_12(1),...
       H_12(end)-H_2f(1), H_2f(end), X_2f(end,6), X_2f(end,1:2)-xf'];

end