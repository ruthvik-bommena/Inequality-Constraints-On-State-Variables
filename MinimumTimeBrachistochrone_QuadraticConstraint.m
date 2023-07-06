function MinimumTimeBrachistochrone_QuadraticConstraint
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     MinimumTimeBrachistochrone_QuadraticConstraint.m
%    Compiler:      MATLAB R2022b
%    Date:          12 June, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve minimum-time Brachistochrone problem
%                   subject to a quadratic path constraint x^2 - 14.6x - 5y ≤ -7.29

clear; close all; clc;

% Parameters
t0 = 0;
x0 = [0; 0; 0.1]; % [x, y, v] (m and m/s) (1 or 0.5)
xf = [10; -10]; % m
g = -9.81; % m/s^2

% Constraint
xPoly = linspace(0,10,1000);
yPoly = (xPoly.^2 - 14.6*xPoly + 7.29)/5;

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-14,'TolX',1e-16,...
    'UseParallel',true); % fsolve

%% Numerical Solution 
lam0_guess = [-0.0447275949571937; 0.189895744132303; -0.00170160736160922; -0.224743937728958; -0.0956247621545226; 0.00301582619216887; -0.0703851052613969; 0.0119380547842369; 0.0198920486791042; 1.28058867989786; -0.274577088050707; 0.00580899463559300; 0.180016342771764; 0.285520506286826; 0.0357916989139358; -0.107562816938759];

iter = 1; err = 1; errTol = 1e-10; iterMax = 100;
    while err > errTol && iter < iterMax
        [lam0,~] = fsolve(@costFunctionConstrained,lam0_guess,options,x0,xf,opts_ode,g,t0);
        err = norm(costFunctionConstrained(lam0_guess,x0,xf,opts_ode,g,t0));
        lam0_guess = lam0;
        iter = iter+1; 
    end
        
% Trajectory
lam0_01 = lam0(1:3); lam0_12 = lam0(4:6); lam0_2f = lam0(7:9);
t1 = lam0(10); t2 = lam0(11); tf = lam0(12);
[~,X_01] = ode89(@unconstrainedBrachistochrone,[t0 t1],[x0; lam0_01],opts_ode,g);
[~,X_12] = ode89(@constrainedBrachistochrone,[t1 t2],[X_01(end,1:3)'; lam0_12],opts_ode,g);
[~,X_2f] = ode89(@unconstrainedBrachistochrone,[t2 tf],[X_12(end,1:3)'; lam0_2f],opts_ode,g);

% Combined Trajectories
X_minT = [X_01; X_12; X_2f];

% Plots
PlotSolution(X_minT,x0,xf,xPoly,yPoly); % plots

end


%% Functions
function PlotSolution(X_minT,x0,xf,xPoly,yPoly)

figure; grid on; hold on; xlim([0 10]);
plot(xPoly,yPoly,'--r','LineWidth',1,'DisplayName','Constraint (x^2 - 14.6x - 5y ≤ -7.29)');
plot(X_minT(:,1),X_minT(:,2),'-k','LineWidth',1,'DisplayName','Trajectory');
xlabel('Distance (m)'); ylabel('Distance (m)'); 
title('Constrained Trajectory of Brachistochrone Problem',['x0 = ',num2str(x0(1:2)'),'  |  ','xf = ',num2str(xf(1:2)')]); 
legend('show','Location','best')

end 


function theta = thetaConstrained(X)

% Point on the curve
xPoint = X(1);
yPoint = X(2);

% Tangent line slope and equation
dydx = (2*xPoint - 14.6)/5;
yLine = 0;
xLine = ((yLine-yPoint)/dydx) + xPoint;

% Tangent vector
tanVec = [xPoint,yPoint]-[xLine,yLine];

% Angle (rad)
theta = atan(abs(tanVec(2)/tanVec(1)));

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


function Xdot = constrainedBrachistochrone(~,X,g)

x = X(1:3);
lambda = X(4:6);

% Control
theta = -thetaConstrained(X);

mu1 = (lambda(1)*x(3)*sin(theta) - lambda(2)*x(3)*cos(theta) - lambda(3)*g*cos(theta))/(-2*x(1)*x(3)*sin(theta) + 14.6*x(3)*sin(theta) - 5*x(3)*cos(theta));

% State Dynamics
xDot = [x(3)*cos(theta);...
        x(3)*sin(theta);...
        g*sin(theta)];

% Costate Dynamics
lambdaDot = [-2*mu1*x(3)*cos(theta);...
             0;...
             -lambda(1)*cos(theta) - lambda(2)*sin(theta) - mu1*(2*x(1)*cos(theta) - 14.6*cos(theta) - 5*sin(theta))];

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


function Hc = hamiltonianConstrained(~,X,g)

    for ii = 1:length(X)

        theta = thetaConstrained(X);

        if size(X,1) == 1
            x = X(1:3);
            lambda = X(4:6);
        else
            x = X(ii,1:3);
            lambda = X(ii,4:6); 
        end
       
        mu1 = (lambda(1)*x(3)*sin(theta) - lambda(2)*x(3)*cos(theta) - lambda(3)*g*cos(theta))/(-2*x(1)*x(3)*sin(theta) + 14.6*x(3)*sin(theta) - 5*x(3)*cos(theta));

        % State Dynamics
        xDot = [x(3)*cos(theta);...
                x(3)*sin(theta);...
                g*sin(theta)];
        
        % Hamiltonian
        Hc(ii) = (lambda * xDot) + 1 + mu1*(2*x(1)*xDot(1) - 14.6*xDot(1) - 5*xDot(2));

    end

end


function err = costFunctionConstrained(lam0_guess,x0,xf,opts_ode,g,t0)

% Control Variables
lam0_01 = lam0_guess(1:3); lam0_12 = lam0_guess(4:6); lam0_2f = lam0_guess(7:9);
t1 = lam0_guess(10); t2 = lam0_guess(11); tf = lam0_guess(12);
piVec1 = lam0_guess(13:14); piVec2 = lam0_guess(15:16);

% Trajectories (unconstrained and constrained arcs)
[t_01,X_01] = ode89(@unconstrainedBrachistochrone,[t0 t1],[x0; lam0_01],opts_ode,g);
[t_12,X_12] = ode89(@constrainedBrachistochrone,[t1 t2],[X_01(end,1:3)'; lam0_12],opts_ode,g);
[t_2f,X_2f] = ode89(@unconstrainedBrachistochrone,[t2 tf],[X_12(end,1:3)'; lam0_2f],opts_ode,g);

% Hamiltonian
[~,H_01] = hamiltonianUnconstrained(t_01,X_01,g);
H_12 = hamiltonianConstrained(t_12,X_12,g);
[~,H_2f] = hamiltonianUnconstrained(t_2f,X_2f,g); 

err = [X_12(1,1)^2 - 14.6*X_12(1,1) - 5*X_12(1,2) + 7.29, X_01(end,4:5)-X_12(1,4:5)-piVec1', X_01(end,6)-X_12(1,6), H_01(end)-H_12(1),...
       X_12(end,4:5)-X_2f(1,4:5)-piVec2', H_12(end)-H_2f(1), X_12(end,6)-X_2f(1,6), X_12(end,1:2)-[6.3,-9],...
       H_2f(end), X_2f(end,6), X_2f(end,1:2)-xf'];

end
