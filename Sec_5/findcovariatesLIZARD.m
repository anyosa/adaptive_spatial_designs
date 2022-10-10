function [H,ids]=findcovariatesLIZARD(domain,nv);

% assign covariates at all nodes
n1=nv(1,1); n2=nv(2,1);
delta=round(domain(1,2)-domain(1,1))/n1;

DATAMAT=csvread('merged_binary_features.csv',1,0);
sites=DATAMAT(:,1:2);
ids12=zeros(size(sites,1),2);
ids=zeros(size(sites,1),1);
for ii=1:size(sites,1),
    ids12(ii,1)=round((sites(ii,1)-domain(1,1))/delta);
    ids12(ii,2)=round((sites(ii,2)-domain(2,1))/delta);
    ids(ii,1)=ids12(ii,1)+(ids12(ii,2)-1)*n1;
end;
H=ones(n1*n2,4);
H(:,2:4)=-10*ones(n1*n2,3);
H(ids,2:4)=DATAMAT(:,4:6);


