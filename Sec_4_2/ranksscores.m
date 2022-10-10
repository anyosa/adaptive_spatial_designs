function ranksscores(score1,score2,score3);

%
% rank scores
%

[rep,crit]=size(score1);

p12=zeros(1,3);
for ii=1:crit,
 [p,h,stats]=ranksum(score1(:,ii),score2(:,ii));
 p12(ii)=p;
 wstat12(ii)=stats.ranksum;
  [p,h,stats]=ranksum(score1(:,ii),score3(:,ii));
   p13(ii)=p;
 wstat13(ii)=stats.ranksum;

   [p,h,stats]=ranksum(score2(:,ii),score3(:,ii));
 p23(ii)=p;
 wstat23(ii)=stats.ranksum;

end
disp('EIBV : ranksum, pvalue');
disp(sprintf('ADAPT EIBV - ADAPT P: %2.0f (%0.2g)',wstat12(1),p12(1))); 
disp(sprintf('ADAPT EIBV - SPAT BAL: %2.0f (%0.2g)',wstat13(1),p13(1))); 
disp(sprintf('ADAPT P - SPAT BAL: %2.0f (%0.2g)',wstat23(1),p23(1)));

disp('Int Misclass prob : ranksum, pvalue');
disp(sprintf('ADAPT EIBV - ADAPT P: %2.0f (%0.2g)',wstat12(2),p12(2))); 
disp(sprintf('ADAPT EIBV - SPAT BAL: %2.0f (%0.2g)',wstat13(2),p13(2))); 
disp(sprintf('ADAPT P - SPAT BAL: %2.0f (%0.2g)',wstat23(2),p23(2)));

disp('Negative log score : ranksum, pvalue');
disp(sprintf('ADAPT EIBV - ADAPT P: %2.0f (%0.2g)',wstat12(3),p12(3))); 
disp(sprintf('ADAPT EIBV - SPAT BAL: %2.0f (%0.2g)',wstat13(3),p13(3))); 
disp(sprintf('ADAPT P - SPAT BAL: %2.0f (%0.2g)',wstat23(3),p23(3)));
