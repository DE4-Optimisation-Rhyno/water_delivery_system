%% Multi-objective Optimisation with Spray Travel and Mass
clc
clear

% Import discrete data for dependent variables from water pump specs
pump_shortlist = csvread('pump_values_shortlist2.csv'); % pump mass, flow rate, pump head
[row,col] = size(pump_shortlist);

% Parameters
g = 9.81;
mu = 0.00089;
rho = 1000;
ep = 0.0001525;

% Testing Variables
r = 0.0762;
l = r*20;
a = 44.7677;

% Declare H and f from the shortlist
Fsl = pump_shortlist(:,2);
Hsl = pump_shortlist(:,3);

% Create arrays for storing T and m values
T_calc = zeros(row,1);
m_calc = zeros(row,1);

% For loop to sub in different pump variables to find T and m
for px = 1:row
    T_calc(px) = max_T(ep,mu,g,r,l,a,Fsl(px),Hsl(px));
    m_calc(px) = min_m(Fsl(px));
end

%% Finding the pareto front for T against m

plot(-m_calc,T_calc,'o')
hold on
title('Trade-off between spray travel and mass using different pumps')
xlabel('Mass (kg)')
ylabel('Spray travel (m)')
set(gcf,'color','w');
hold off

%% Functions for the two objectives

%Calculating T - vars r,l,a,F,H
function T = max_T(ep,mu,g,r,l,a,F,H)
X = (ep./(7.4.*r)).^1.1 + ((241.*r)./(25000.*F));
Hf =((1./(-1.8.*log10(X))).^2.*mu.*l.*F.^2)./(pi.^2.*r.*3.*g);
T = sin(2.*a).*-mu.*(H-l.*sin(a)-((F.*g)./(2.*pi.*r))-Hf);
end

%Calculating m - var F
function m = min_m(F)
F_s = F/100;
m = 0.0087.*F_s.^3 -0.1663.*F_s.^2 + 1.9751.*F_s + 24.908;
end