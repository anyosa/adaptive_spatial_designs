function [srs,chivec]=BevalDes(H,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2);

%
%
% File with variance reductions in 1 step 
% design of the linearized likelihood model.
%
%

[n1,n2]=size(cm);
n=n1*n2;
x0v=H*mubeta;
HD=H(nodes,:);

mb0Z=x0v(nodes,1);
mb0=vec2mat(x0v,n1);
P=zeros(size(nodes,1),size(nodes,1)); % Diagonal likelihood matrix
B_sp=zeros(size(nodes,1),size(nodes,1)); % B=A * Sigma * A', using design and prior
% Set the B matrix from correlation structure and measurement nodes
for i=1:size(nodes,1),
  for j=i:size(nodes,1),
    B_sp(i,j)=cm(1+mod2(nodes1(j,1),nodes1(i,1),n1),1+mod2(nodes2(j,1),nodes2(i,1),n2));
    if (j>i)
      B_sp(j,i)=B_sp(i,j);
    end;
  end;
end;

B=B_sp+HD*Sigmabeta*HD';

HSD=H*Sigmabeta*HD';
  
% linearization point 
x0vn=x0v(nodes,1);
% derivatives of link
bder=exp(x0vn)./(1+exp(x0vn));
bder2=exp(x0vn)./((1+exp(x0vn)).^2);
P=diag(1./bder2);  
% Marginal covariance, with latent state x integrated out
R=B+P; 
Ri=inv(R);
Lit=chol(Ri)';
VRmat=zeros(n,size(nodes,1));
for kk=1:size(nodes,1),
    % Spatial part, Fourier domain
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
chivec=sum(VRmat.^2,2);
srs0=cm(1,1)*ones(n,1)+sum(HVV.^2,2)-chivec;

srs=max(0.0001,srs0);

