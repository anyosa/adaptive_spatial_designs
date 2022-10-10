function [P,z,u,Ri,mvec,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);

%
% Iterative linearization fitting Gaussian approximative posterior. And for covariance
% matrix of the linearized likelihood model.
%

[n1,n2]=size(cm);
n=n1*n2;
x0v=H*mubeta;
HD=H(nodes,:);
HSD=HSigma*HD';

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

convI=0;
cntI=0;
while (convI==0),
  
  %disp('Vector of latent variables'); 
  x0vn=x0v(nodes,1);
  
  % Transformed measurement z and covariance matrix C
  z=zeros(size(nodes,1),1);
  bder=exp(x0vn)./(1+exp(x0vn));
  bder2=exp(x0vn)./((1+exp(x0vn)).^2);
  z=(y-bder+(x0vn.*bder2))./bder2;
  P=diag(1./bder2);  
  % Marginal covariance, with latent state x integrated out
  R=B+P; 
  % Compute (A*Sigma A'+P)\(z-H_D*beta)
  ufo=R\(z-mb0Z);
  % Allocate and calculate u which is Ri*(z-A*mu), data impact via z
  uv=zeros(n,1);
  uv(nodes,1)=ufo;
  u=zeros(n1,n2);
  u=vec2mat(uv,n1);
  % Reset the point of linearization
  % Newton step,
  xl=mb0+real(fft2(cmf.*ifft2(u))); % Spatial part, Fourier domain
  xlv=mat2vec(xl);
  xlv=xlv+HSD*ufo; % Regression part
  dista=sum((xlv-x0v).*(xlv-x0v)); % Distance
  x0v=xlv;
  %mean(mean(xl))
  %pause(2);
  if (dista<0.0000001)
      convI=1;
      disp(sprintf('Number of iterations is %d \n',cntI));
      disp(sprintf('Change is %0.5e \n',dista));
  end
  cntI=cntI+1;
  if (cntI>20)
      convI=1;
  end;
end;
% Return the posterior mode of the approximation
mvec=x0v;

% calculating variances
Ri=inv(R);
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

