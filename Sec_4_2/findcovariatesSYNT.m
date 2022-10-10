function [H,ids]=findcovariatesSYNT(domain,nv);

% assign covariates at all nodes
n1=nv(1,1); n2=nv(2,1);
delta=round(domain(1,2)-domain(1,1))/n1;

range=(domain(:,2)-domain(:,1));
mrange=mean(range);
delta=range./nv;
gE=domain(1,1)+delta(1,1)*[1:n1]';
gN=domain(2,1)+delta(2,1)*[1:n2]';
gM=[mean(gE);mean(gN)];
hE=(gE-mean(gE)).^2;
hN=(gN-mean(gN)).^2;
disst=hE*ones(1,n2)+ones(n1,1)*hN';
Hm=exp(-disst/(2*(0.23*mrange)^2));
Gm=((gE-min(gE))./(max(gE)-min(gE)))*ones(1,n2);
%H=[ones(n1*n2,1) mat2vec(Hm)];
H=[ones(n1*n2,1) mat2vec(Gm)];
ids=1;
