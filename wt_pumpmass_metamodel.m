%% Metamodel for correlating pump mass and flow rate
clc
clear

% Import discrete data for dependent variables from water pump specs
pump_data = csvread('pump_values2.csv'); % pump mass, flow rate, pump head

mp = pump_data(:,1);
Fp = pump_data(:,2);

% Scaling Fp to fit with m values
Fp_s = Fp/100;

p = polyfit(Fp_s,mp,3); % approximates relationship to n degrees

Fp_s_lb = min(Fp_s);
Fp_s_ub = max(Fp_s);
Fp1 = linspace(Fp_s_lb, Fp_s_ub, 100);
mp1 = polyval(p, Fp1);

plot(Fp_s,mp,'o'); % plots scatter graph for flow rate and mass
hold on
grid on
plot(Fp1,mp1); % plots approximated relationship
approx = sprintf('y = (%.3f) x^3 + (%.3f) x^2 + (%.3f) x + (%.1f)',p(1),p(2),p(3),p(4));
text(1,45,approx);
title('Relationship between flow rate and mass of a pump')
xlabel('Scaled Flow rate (0.01L/min)')
ylabel('Spray travel (m)')
hold off

%% Calculating pump mass, using flow rate as input rate

% Variables
F = 750;
F_s = 750/100;
% Objective Function

PumpMass = p(1).*F_s.^3 + p(2).*F_s.^2 + p(3).*F_s + p(4)