function ParameterSpecification;

%
%
%% Parameters
%
% bethat =
% 
%     4.7152
%   -14.3388
%    -0.0309
% 
% 
% Sigbet =
% 
%     9.5828   -3.1267   -0.1469
%    -2.4468    2.8973    0.0347
%    -0.1478    0.0454    0.0023
% 
%
% exp(thetahat)=(sigma0,ff,tau0)
%     0.7125
%     0.0556
%     0.0001
%


syntdata=1;

if (syntdata==1)
    [y,Hfull,sites,ddist]=gpsim;
    n=size(y,1);
    H=Hfull(:,2:end);
    A=ones(size(H));
        
end

figure(1); 
clf; 
hold;
for i=1:size(sites,1),
    sp=plot(sites(i,2),sites(i,1),'ko');
    set(sp,'MarkerSize',0.01+50*abs(H(i,1)));    
end
title('Circle size indicate covariate');
xlabel('longitude'); ylabel('latitude');

figure(2); 
clf; 
hold;
for i=1:size(sites,1),
    sp=plot(sites(i,2),sites(i,1),'ko');
    set(sp,'MarkerSize',0.01+0.001*A(i,1));    
end
title('Circle size indicate area');
%[d1km d2km]=lldistkm(latlon1,latlon2)
xlabel('longitude'); ylabel('latitude');

figure(3); 
clf; 
hold;
for i=1:size(y,1),
    if (y(i,1)==1)
        plot(sites(i,2),sites(i,1),'gx');
    else
        plot(sites(i,2),sites(i,1),'ro');       
    end
end
title('blue is 1 and red is 0');
xlabel('longitude'); ylabel('latitude');


%
% Simple Logit Reg
%

[b,dev,stats] = glmfit(H,y,'binomial','link','logit');

SS=stats.covb;

for ii=1:size(SS,1),
    disp(sprintf('Conf interval for variable %d',ii))
    [b(ii,1)-1.96*sqrt(SS(ii,ii)) b(ii,1)+1.96*sqrt(SS(ii,ii))] 
end

%
% initial value for covariance and regression
% set from simple plots
%

figure(7); 
clf;
imagesc(ddist); colorbar;

ssm=mean(diag([ones(size(H,1),1) H]*SS*[ones(size(H,1),1) H]'));
sigma0=sqrt(0.5*ssm);
ff=20;
muB0=b;
SigmaB0=SS;
tau=0.001;

theta0=[log(sigma0);log(ff);log(tau)];

%
% maximizing the approximate likelihood by Laplace approximation
%
Hr=[ones(size(y,1),1) H];
a=[Hr y];
funf=@(x)likeval(x,a,ddist,muB0,SigmaB0);

ff00=feval(funf,theta0);
options=optimset('Display','iter','MaxIter',500,'TolFun',0.0001);

[x,fval]=fminsearch(funf,theta0,options);%

disp('estimates:');
thetahat=x;
exp(thetahat)

disp('starting values:');
exp(theta0)

%% Posterior of beta 

n=size(y,1);
% mean
mu=[zeros(n,1);muB0]; 

% Covariance
%spat
sigma0=exp(thetahat(1,1));
ff=exp(thetahat(2,1));
tau=exp(thetahat(3,1));
SigmaS=sigma0^2*(1+ff*ddist).*exp(-ff*ddist)+tau^2*eye(n);
p=size(SigmaB0,1);
Sigma=[SigmaS zeros(n,p);zeros(p,n) SigmaB0];

% Approx Gauss
[mhat,Sigmahat]=approxGauss(mu,Sigma,Hr,y);

bethat=mhat(end-2:end,1)
Sigbet=Sigmahat(end-2:end,end-2:end)

sthat=sqrt(diag(Sigmahat(1:n,1:n)));

figure(11);
clf;
plot(mhat(1:n,1)); hold;
plot(mhat(1:n,1)+1.64*sthat(1:n,1),'r--');
plot(mhat(1:n,1)-1.64*sthat(1:n,1),'r--');
title('posterior latent process'); 

