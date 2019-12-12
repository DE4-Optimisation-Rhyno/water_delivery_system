%% Water Transmission Subsystem Optimisation using non-gradient methods
clc
clear
tic
% x(1) is pipe radius r
% x(2) is pipe angle a
% x(3) is pump flow rate F
% pipe length and pump pressure head are substituted from simplification

%% Formulation

% Parameters
g = 9.81;
mu = 0.00089;
rho_water = 1000;
H = 12.5; %Average set value set via parametric analysis
ep = 0.0001525;

% Initial Values
x0 = [0.02,70,900]; %taking the average value informed by lower and upper bounds

% Lower and Upper Variable Bounds from Constraints
lb = [0.0762,10,500];
ub = [0.35,80,1000];

% Linear Constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear Constraints
nonlincon = @wt_nlcon;

% Objective Function in negative null form

objective = @(x) -sin(2.*x(2)).*-mu.*(H-20.*x(1).*sin(x(2))-((x(3).*g)./(2.*pi.*x(1)))-((1./(-1.8.*log10((ep./(7.4.*x(1))).^1.1 + ((241.*x(1))./(25000.*x(3)))))).^2.*mu.*20.*x(1).*x(3).^2)./(pi.^2.*x(1).*3.*g));

%% Optimisation using genetic algorithm
nvars = 3;
x_ga = ga(objective,nvars,A,b,Aeq,beq,lb,ub,nonlincon)

% Display solution
disp('Genetic Algorithm Solution')
disp(['r = ' num2str(x_ga(1))])
disp(['a = ' num2str(x_ga(2))])
disp(['F = ' num2str(x_ga(3))])
toc
%% Non-linear Constraints 

function [c,ceq] = wt_nlcon(x)
    c = [Constraint_g1(x),Constraint_g6(x),Constraint_g11(x)];
    ceq = 0;
end

% Constraint g1: Reynold's number cannot exceed laminar flow
function g1 = Constraint_g1(x)
    % Paramters
    mu = 0.00089;
    rho_water = 1000;
    
    % Variables
    r = x(1);
    F = x(3);
    
    % Constraint equation
    g1 = ((2.*rho_water.*F.*mu.^2)./(pi.*r))-4000;

end

% Constraint g6: Jet reaction force cannot exceed 809N
function g6 = Constraint_g6(x)
    % Paramters
    mu = 0.00089;
    
    % Variables
    r = x(1);
    a = x(2);
    F = x(3);
    
    % Constraint equation
    g6 = ((F.*mu.*cos(a))/(pi.*r.^2))-809;

end

% Constraint g11: Keeping full pipe via ratio between radius and flow rate
function g11 = Constraint_g11(x)
    % Paramters
    rho_water = 1000;

    % Variables
    r = x(1);
    F = x(3);
    
    % Constraint equation
    g11 = (F/r.^0.25)-1.5.*rho_water;

end

