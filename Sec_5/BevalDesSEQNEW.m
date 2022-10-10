function [srs,chivec]=BevalDesSEQNEW(mOLD,H,Sigmabeta,cmf,cm,srsOLD,Risub,nodesNEW,nodes1NEW,nodes2NEW,nodesOLD,nodes1OLD,nodes2OLD);

%
% File with variance reductions in SEQUENTIAL 
% design of the linearized likelihood model.
%
%
% INPUTS:
%
% H = matrix of covariates
% Sigmabeta = matrix of regression parameter covariance
% cmf = Fourier domain spatial covariance
% cm = Toeplitz matrix holding spatial covariance
% mOLD = old mean vector
% srsOLD = old variance vector
% RiOLD = old inverse measurement covariance
% nodesOLD = sampled grid cells so far
% nodes NEW = design grid cells for next float
%
%
%
% OUTPUTS:
%
% srs = approximate marginal posterior variance for all grid locations
% chivec = approximate marginal reduction in variance for all grid locations
%

nNEW=size(nodesNEW,1);
nOLD=size(nodesOLD,1);

nodes=[nodesOLD;nodesNEW];
nodes1=[nodes1OLD;nodes1NEW];
nodes2=[nodes2OLD;nodes2NEW];

[n1,n2]=size(cm);
n=n1*n2;

%x0v=H*mubeta;
%mb0Z=x0v(nodes,1);
%mb0=vec2mat(x0v,n1);

P22=zeros(size(nodesNEW,1),size(nodesNEW,1)); % Diagonal likelihood matrix
B_sp22=zeros(size(nodesNEW,1),size(nodesNEW,1)); % B=A * Sigma * A', using design and prior
% Set the B matrix from correlation structure and measurement nodes
for i=1:size(nodesNEW,1),
  for j=i:size(nodesNEW,1),
    B_sp22(i,j)=cm(1+mod2(nodes1NEW(j,1),nodes1NEW(i,1),n1),1+mod2(nodes2NEW(j,1),nodes2NEW(i,1),n2));
    if (j>i)
      B_sp22(j,i)=B_sp22(i,j);
    end;
  end;
end;
B_sp22=B_sp22+H(nodesNEW,:)*Sigmabeta*H(nodesNEW,:)';

B_sp12=zeros(size(nodesOLD,1),size(nodesNEW,1)); % B=A * Sigma * A', using design and prior
% Set the B matrix from correlation structure and measurement nodes
for i=1:size(nodesOLD,1),
  for j=1:size(nodesNEW,1),
    B_sp12(i,j)=cm(1+mod2(nodes1NEW(j,1),nodes1OLD(i,1),n1),1+mod2(nodes2NEW(j,1),nodes2OLD(i,1),n2));
   end;
end;
R12=B_sp12+H(nodesOLD,:)*Sigmabeta*H(nodesNEW,:)';

HSD=H*Sigmabeta*H(nodes,:)';
  
% linearization point to form Gaussian approximation
x0vn=mOLD(nodesNEW,1);
% derivatives of link
bder=exp(x0vn)./(1+exp(x0vn));
bder2=exp(x0vn)./((1+exp(x0vn)).^2);
P22=diag(1./bder2);  
R22=B_sp22+P22; 

Qi=inv(R22-R12'*Risub*R12);
Ri=[Risub+Risub*R12*Qi*R12'*Risub -Risub*R12*Qi;-Qi*R12'*Risub  Qi];

eigri=min(eig(Ri));
if (eigri<0.00000001)
    Ri=diag(max(0.001,diag(Ri)));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% srs part
%%
Lit=chol(Ri)';
VRmat=zeros(n,size(nodes,1));
% Spatial part, Fourier domain, 
% for each row of data Cholesky factor
for kk=1:size(nodes,1),
    uu=Lit(:,kk);
    uw=zeros(n,1);
    uw(nodes,1)=uu;
    uwkk=vec2mat(uw,n1);
    WW=real(fft2(cmf.*ifft2(uwkk))); 
    VRmat(:,kk)=mat2vec(WW)+HSD*uu;  
end
% Regression part
LSigma=chol(Sigmabeta)';
HVV=H*LSigma;

% Variances of x for these parameters sigma and nu
srs0=cm(1,1)*ones(n,1)+sum(HVV.^2,2)-sum(VRmat.^2,2);
srs=max(0.0001,srs0);

chivec0=srsOLD-srs;
chivec=max(0.0001,chivec0);

