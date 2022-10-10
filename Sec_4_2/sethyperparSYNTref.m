function [sig,nu,ksi,mubeta,Sigmabeta,gamma]=sethyperparSYNTref;

% Prior Gaussian process settings
% Would depend on spatial smoothness/ heterogeneity
%
% Note, could also read the reference values here from file: 
% ReplicREF1.csv,...

% weight of variance
%gamma=0.1;
%gamma=0.5;
%gamma=0.3;
%gamma=0.9;
%gamma=2;
%gamma=3;
gamma=0.05;

% spat uncertainty
gamma1ref=1;
gamma1=gamma1ref;%-0.9+1.8*rand;
% regr uncertainty
gamma2ref=1;
gamma2=gamma2ref;%-0.9+1.8*rand;
sig0=1;  % standard deviation

sig=sig0*gamma1;

nuref=1500;  % Range in m
nu0=nuref;%-500+1000*rand;
nu=nu0;
%nu=nuref;

ksi=0.01; % nugget
 
% Bayesian regression inputs
% Would depend on covariates
mubeta=[-2;4];
SB12ref=-0.7;
SB11=1; SB22=1; SB12=SB12ref;%+0.2-0.4*rand;
Sigmabeta=gamma2^2*[SB11^2 SB12*SB11*SB22;SB12*SB11*SB22 SB22^2];


