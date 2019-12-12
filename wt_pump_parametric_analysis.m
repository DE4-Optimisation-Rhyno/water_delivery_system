%% Initial parametric analysis of water pumps
clc
clear

% Parameters
g = 9.81;
mu = 0.00089;
rho_water = 1000;
ep = 0.0001525;

% Variables (temporarily set to constants, based on average values from
% lower and upper bounds of constraints
r = 0.02;
l = 0.6;
a = 45;

% Testing values

F_upper = 920;
F_lower = 550;

H_upper = 20;
H_lower = 5;

% Averages

F = 735;
H = 12.5;

%% Plotting pump flow rate against spray travel T

F_values = linspace(F_lower,F_upper,100);

X = (ep./(7.4.*r)).^1.1 + ((241.*r)./(25000.*F_values));
Hf =((1./(-1.8.*log10(X))).^2.*mu.*l.*F_values.^2)./(pi.^2.*r.*3.*g);
T = sin(2.*a).*-mu.*(H-l.*sin(a)-((F_values.*g)./(2.*pi.*r))-Hf);

plot(F_values,T)
title('Effect of pump flow rate on spray travel')
xlabel('Flow rate (L/min)')
ylabel('Spray travel (m)')
set(gcf,'color','w');

%% Plotting pump pressure head against spray travel T
    
H_values = linspace(H_lower,H_upper,100);

X = (ep./(7.4.*r)).^1.1 + ((241.*r)./(25000.*F));
Hf =((1./(-1.8.*log10(X))).^2.*mu.*l.*F.^2)./(pi.^2.*r.*3.*g);
T = sin(2.*a).*-mu.*(H_values-l.*sin(a)-((F.*g)./(2.*pi.*r))-Hf);

figure();
plot(H_values,T)
title('Effect of pump pressure head on spray travel')
xlabel('Pressure head (m)')
ylabel('Spray travel (m)')
set(gcf,'color','w');