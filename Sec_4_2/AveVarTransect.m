function avtransect=AveVarTransect(mu,sigs2,nodes);

alfa=0.58;

r=alfa*mu(nodes,1)./sqrt(1+alfa^2*sigs2(nodes,1));
avtransect=mean(normcdf(r).*normcdf(-r));
