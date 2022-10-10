function tstatdiff(score1,score2,score3);


B=size(score1,1);
diff12=score1-score2;
mdiff12=mean(diff12,1);
sdiff12=std(diff12,1);
tstat12=sqrt(B)*mdiff12./sdiff12;
diff13=score1-score3;
mdiff13=mean(diff13,1);
sdiff13=std(diff13,1);
tstat13=sqrt(B)*mdiff13./sdiff13;
diff23=score2-score3;
mdiff23=mean(diff23,1);
sdiff23=std(diff23,1);
tstat23=sqrt(B)*mdiff23./sdiff23;

disp('EIBV t-stat, pvalue');
disp(sprintf('ADAPT EIBV - ADAPT P: %0.3g (%0.3g)',tstat12(1),normcdf(-abs(tstat12(1))))); 
disp(sprintf('ADAPT EIBV - SPAT BAL: %0.3g (%0.3g)',tstat13(1),normcdf(-abs(tstat13(1))))); 
disp(sprintf('ADAPT P - SPAT BAL: %0.3g (%0.3g)',tstat23(1),normcdf(-abs(tstat23(1)))));

disp('Int Misclass prob : t-stat');
disp(sprintf('ADAPT EIBV - ADAPT P: %0.3g (%0.3g)',tstat12(2),normcdf(-abs(tstat12(2))))); 
disp(sprintf('ADAPT EIBV - SPAT BAL: %0.3g (%0.3g)',tstat13(2),normcdf(-abs(tstat13(2))))); 
disp(sprintf('ADAPT P - SPAT BAL: %0.3g (%0.3g)',tstat23(2),normcdf(-abs(tstat23(2)))));

disp('Negative log score - tstat');
disp(sprintf('ADAPT EIBV - ADAPT P: %0.3g (%0.3g)',tstat12(3),normcdf(-abs(tstat12(3))))); 
disp(sprintf('ADAPT EIBV - SPAT BAL: %0.3g (%0.3g)',tstat13(3),normcdf(-abs(tstat13(3))))); 
disp(sprintf('ADAPT P - SPAT BAL: %0.3g (%0.3g)',tstat23(3),normcdf(-abs(tstat23(3)))));


%plotksdiff(diff12,diff13,diff23);
