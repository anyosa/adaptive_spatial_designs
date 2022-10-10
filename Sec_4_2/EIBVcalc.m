function IBV=EIBVcalc(mus,sigs,chis);

%
% Calculation of the expected integrated Bernoulli variance
% The file relies on a bi-variate cumulate distribution function
% The mean and variance required is computed from the expected
% values (mus), variance (sigs) and variance reduction (chis) terms in the inputs. 
% 

alfa=0.58; % Used to fit logit to probit link (Demidenko book)

r=alfa/sqrt(1+alfa^2*sigs^2);
cpar=-r^2*chis^2/(1+r^2*chis^2);
mu1=r*mus/sqrt(1+r^2*chis^2);
mu2=-r*mus/sqrt(1+r^2*chis^2);

IBV=mvncdf([mu1;mu2],[0;0],[1 cpar;cpar 1]);


