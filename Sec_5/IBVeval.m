function IBVval=IBVeval(mu,sigs2,nodes);

alfa=0.58;

r=alfa*mu./sqrt(1+alfa^2*sigs2);
IBVvalV=normcdf(r).*normcdf(-r);
IBVval=0;
for i=1:size(mu,1),
    ff=ismember(i,nodes);
    if (ff==1)
        IBVval=IBVval;
    else
        IBVval=IBVval+IBVvalV(i,1);
    end
end

figure(12); 
clf;
plot(IBVvalV);


