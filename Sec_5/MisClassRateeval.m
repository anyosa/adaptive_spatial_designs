function MCrate=MisClassRateeval(mu,sigs2,nodes,yTRUE);

alfa=0.58;

r=alfa*mu./sqrt(1+alfa^2*sigs2);
PValV=normcdf(r);
MCRatevalV=zeros(size(mu,1),1);
MCrate=0;
for i=1:size(mu,1),
    ff=ismember(i,nodes);
    if (ff==1)
        MCrate=MCrate;
    else
        % misclassification rates
        MCRatevalV(i,1)=PValV(i,1)*(yTRUE(i,1)==0)+(1-PValV(i,1))*(yTRUE(i,1)==1);
        MCrate=MCrate+MCRatevalV(i,1);
    end
end

%figure(42); 
%clf;
%plot(MCRatevalV);
