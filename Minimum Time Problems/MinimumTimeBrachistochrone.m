function MinimumTimeBrachistochrone
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     MinimumTime_Brachistochrone.m
%    Compiler:      MATLAB R2022b
%    Date:          10 June, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve minimum-time Brachistochrone problem
%    References:    Ch 3. Applied Optimal Control, 1975, A.E. Bryson. Jr, Yu-Chi Ho

clear; close all; clc;

% Parameters
t0 = 0;
x0 = [0; 0; 0]; % [x, y, v] (m and m/s)
xf = [10; -10]; % m
g = -9.81; % m/s^2

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-14,'TolX',1e-16,...
    'UseParallel',false); % fsolve

%% Numerical Solution 
lam0_guess = [-0.066694440231030; 0.025469424836357; -0.101936799184506; 1.843277301347743];

iter = 1; err = 1; errTol = 1e-10; iterMax = 100;
    while err > errTol && iter < iterMax
        [lam0,~] = fsolve(@costFunctionUnconstrained,lam0_guess,options,x0,xf,opts_ode,g,t0);
        err = norm(costFunctionUnconstrained(lam0_guess,x0,xf,opts_ode,g,t0));
        lam0_guess = lam0;
        iter = iter+1; 
    end

% Trajectory
tf = lam0(4);
[~,X_minT] = ode89(@unconstrainedBrachistochrone,[t0 tf],[x0; lam0(1:3)],opts_ode,g);

% Plots
PlotSolution(X_minT,x0,xf); % plots

end


%% Functions
function PlotSolution(X_minT,x0,xf)

figure; grid on; hold on; xlim([0 10]);
plot(X_minT(:,1),X_minT(:,2),'-k','LineWidth',1,'DisplayName','Trajectory');
xlabel('Distance (m)'); ylabel('Distance (m)'); 
title('Trajectory of Brachistochrone Curve',['x0 = ',num2str(x0(1:2)'),'  |  ','xf = ',num2str(xf(1:2)')]); 
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


function err = costFunctionUnconstrained(lam0_guess,x0,xf,opts_ode,g,t0)

% Control Variables
lam0 = lam0_guess(1:3); 
tf = lam0_guess(4);

% Trajectory
[t,X] = ode89(@unconstrainedBrachistochrone,[t0 tf],[x0; lam0],opts_ode,g);

% Hamiltonian
[~,H] = hamiltonianUnconstrained(t,X,g);

err = [X(end,1:2)-xf', X(end,6), H(end)];

end
