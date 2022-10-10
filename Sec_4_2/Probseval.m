function probspred=Probseval(mu,sigs2,nodes,yTRUE);

alfa=0.58;

r=alfa*mu./sqrt(1+alfa^2*sigs2);
probs=normcdf(r);
for i=1:size(mu,1),
    ff=ismember(i,nodes);
    if (ff==1)
        probspred(i,1)=yTRUE(i,1);
    else
        probspred(i,1)=probs(i,1);
    end
end
