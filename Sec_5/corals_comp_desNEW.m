function [EIBV]=corals_comp_desNEW(H,Sigmabeta,cmf,cm,mubeta,nodes,nodes1,nodes2);

%%%%%%%%%%%%%%%%%%%%
%
% BERNOULLI data, gaussian latent intensity field
% 1) Generation of synthetic data and print data to file
% 2) Estimation of parameters and prediction
%
%%%%%%%%%%%%%%%%%%%%

%figure(6); clf; imagesc(vec2mat(H(:,2),n1)'); colorbar;
% Form H_{x,D}
HSigma=H*Sigmabeta;

muvec=H*mubeta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE EIBV OF DESIGNS
%%%%%%%%%%%%

n=size(H,1);
mmean=zeros(n,1); % Allocate mean
srs=zeros(n,1); % Allocate variance

%toc;
%disp('Building posterior covariance matrices');

% Gaussian approximation (see separate file below)
[srs,chivec]=BevalDes(H,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2);

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
    ff=ismember(i,nodes);
    if (ff==1)
        EBV=0;
    else
        mus=muvec(i,1);
        ssi=sigs(i,1);
        chii=chis(i,1);
        EBV=EIBVcalc(mus,ssi,chii);
    end
    EIBVs(i,1)=EBV;
end;
EIBV=sum(EIBVs);
disp(sprintf('EIBV: %3.2f',EIBV));
