function [muhat,Sigmahat]=approxGauss(mu,C,H,y);

%
%
% Evaluation of the approximate Normal, 
%
% Inputs are prior mean mu and covariance matrix C and bernoulli data y
%
% 
%
% Outputs are the parameters in the conditional x|y
% ** muhat,Sigmahat
% 



k=size(y,1);
n=k;
A=[eye(n) H];
x0=mu;
Amu=A*mu;
niter=10;   % NUMBER OF ITERATIONS FOR MODE FITTING
for i=1:niter, 
    etaold=A*x0;
    % Transformed measurement z and covariance matrix C
    z=zeros(n,1);
    bder=(exp(etaold))./(ones(n,1)+exp(etaold));
    bder2=(exp(etaold))./((ones(n,1)+exp(etaold)).^2);
    z=(y-bder+(etaold.*bder2))./bder2;
    P=diag(1./bder2);  
    if (i<niter)
        AC=A*C;
        xl=mu+AC'*((AC*A'+P)\(z-Amu)); % Newton step
        dista=sum(real((xl-x0).*(xl-x0))); % Distance
        if (i==niter-1)
          disp(sprintf('Number of iterations is %d \n',i));
          disp(sprintf('Change is %0.5e \n',dista));
        end;
        x0=xl;
    end;
end;
% Return the posterior mode of the approximation
muhat=x0;

AC=A*C;
Sigmahat=C-AC'*((AC*A'+P)\AC);

