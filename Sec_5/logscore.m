function  LS=logscore(mu,sigs2,nodes,yGEN);

alfa=0.58;

r=alfa*mu./sqrt(1+alfa^2*sigs2);
PValV=normcdf(r);
LS=0;
for i=1:size(mu,1),
    ff=ismember(i,nodes);
    if (ff==1)
        gg=1;
    else
        LSv=-log(PValV(i,1))*yGEN(i,1)-log(1-PValV(i,1))*(1-yGEN(i,1));
        LS=LS+LSv;
    end
end
