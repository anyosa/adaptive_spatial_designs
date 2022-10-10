function ll=logdet(f);


%L=chol(f);
%ll=2*sum(log(diag(L)));
[d,l]=eig(f);
ll=sum(log(diag(l)));

