function WCProb=WrongClassProbeval(mu,sigs2,nodes,yTRUE);

alfa=0.58;

r=alfa*mu./sqrt(1+alfa^2*sigs2);
PValV=normcdf(r);
WCProbvalV=zeros(size(mu,1),1);
WCProb=0;
for i=1:size(mu,1),
    ff=ismember(i,nodes);
    if (ff==1)
        WCProb=WCProb;
    else
        WCProbvalV(i,1)=PValV(i,1)*(yTRUE(i,1)==0)+(1-PValV(i,1))*(yTRUE(i,1)==1);
        WCProb=WCProb+WCProbvalV(i,1);
    end
end

%figure(22); 
%clf;
%plot(WCProbvalV);