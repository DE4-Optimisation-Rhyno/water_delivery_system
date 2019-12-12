%% Calculating spray travel
clc
clear

% Parameters
g = 9.81;
mu = 0.00089;
rho_water = 1000;
ep = 0.0001525;
%{
% Testing Variables
r = 0.0762;
l = r*20;
a = 44.7677; %44.7677/41.6261/57.3341
F = 788.0978;
H = 12.5;
%}

r = 0.0762;
l = r*20;
a = 44.7677; %44.7677/41.6261/57.3341
F = 970;
H = 4;

% Objective Function for Spray Travel
X = (ep./(7.4.*r)).^1.1 + ((241.*r)./(25000.*F));
Hf =((1./(-1.8.*log10(X))).^2.*mu.*l.*F.^2)./(pi.^2.*r.*3.*g);
T = sin(2.*a).*-mu.*(H-l.*sin(a)-((F.*g)./(2.*pi.*r))-Hf);

SprayTravel = T