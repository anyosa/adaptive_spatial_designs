function [sig,nu,ksi,mubeta,Sigmabeta]=sethyperparLIZARD;

% Prior Gaussian process settings
% Would depend on spatial smoothness/ heterogeneity
%
sig=3;  % standard deviation
nu=20;  % Range in m
ksi=0.001; % nugget

% Bayesian regression inputs
% Would depend on covariates
%
mubeta=[3;0.6;0.03;0.002];
Sigmabeta=[0.207 0.028   -0.0007   -0.0002;...
   0.0283    0.0055    0.0001    0.0000;...
   -0.0007    0.0001    0.0001    0.0000;...
   -0.0002    0.0000    0.0000    0.0001];...

