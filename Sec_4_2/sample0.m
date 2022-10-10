
function x=sample0(cmf);

[n1,n2]=size(cmf);
xx=real(fft2(sqrt(cmf).*ifft2(randn(n1,n2))));
x=reshape(xx,n1*n2,1);

