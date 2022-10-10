function [EIBV]=corals_comp_SEQdesNEW(H,Sigmabeta,cmf,cm,mOLD,srsOLD,RiOLD,nodesNEW,nodes1NEW,nodes2NEW,nodesOLD,nodes1OLD,nodes2OLD);
%%%%%%%%%%%%%%%%%%%%
%
% BERNOULLI data, gaussian latent intensity field
% 1) Generation of synthetic data and print data to file
% 2) Estimation of parameters and prediction
%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE EIBV OF DESIGNS
%%%%%%%%%%%%

n=size(H,1);
mmean=zeros(n,1); % Allocate mean
srs=zeros(n,1); % Allocate variance

%toc;
%disp('Building posterior covariance matrices');

% Gaussian approximation (see separate file)
[srs,chivec]=BevalDesSEQNEW(mOLD,H,Sigmabeta,cmf,cm,srsOLD,RiOLD,nodesNEW,nodes1NEW,nodes2NEW,nodesOLD,nodes1OLD,nodes2OLD);

% Arrange standard deviation on the spatial grid
sigs=sqrt(srs);
chis=sqrt(chivec);

%%%%%%%%%%%%%%%%%%
% Compute EIBV
%%%%%%%
%toc;
%disp('Computing EIBV');

EIBVs=zeros(n,1);
for i=1:n,
    ff=ismember(i,[nodesOLD;nodesNEW]);
    if (ff==1)
        EBV=0;
    else
        mus=mOLD(i,1);
        ssi=sigs(i,1);
        chii=chis(i,1);
        EBV=EIBVcalc(mus,ssi,chii);
    end
    EIBVs(i,1)=EBV;
end;
EIBV=sum(EIBVs);
disp(sprintf('EIBV: %0.3g',EIBV));

