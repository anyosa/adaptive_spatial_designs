function centralindeces=getCentralInd(n1,n2,n);

%
% Cutting parts at all ends of grid
%

n1g=reshape([1:n1]'*ones(1,n2),n,1);
n2g=reshape(ones(n1,1)*[1:n2],n,1);
h1=ceil(n1*0.1);
h2=ceil(n2*0.1);
centralindeces=find((n1g>h1)&(n1g<n1-h1+0.1)&(n2g>h2)&(n2g<n2-h2+0.1));
