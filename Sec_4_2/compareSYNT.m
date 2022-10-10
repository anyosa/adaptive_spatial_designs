function compareSYNT;


rng(7);

% true model
[scorevecADAPT,scorevecP,scorevecSPAT]=ReplicateSynt;

% synt
[scorevecADAPTref,scorevecPref,scorevecSPATref]=ReplicateSyntref;

% t-tests
disp('--------------------');

tstatdiff(scorevecADAPT,scorevecP,scorevecSPAT)

tstatdiff(scorevecADAPTref,scorevecPref,scorevecSPATref)

B=100;
%B=2;
diff12=scorevecADAPT-scorevecADAPTref;
mdiff12=mean(diff12,1);
sdiff12=std(diff12,1);
tstat12=sqrt(B)*mdiff12./sdiff12;

disp('EIBV t-stat, pvalue');
disp(sprintf('ADAPT EIBV - ADAPT EIBV ref: %0.3g (%0.3g)',tstat12(1),normcdf(-abs(tstat12(1))))); 

disp('Int Misclass prob : t-stat');
disp(sprintf('ADAPT EIBV - ADAPT EIBV ref: %0.3g (%0.3g)',tstat12(2),normcdf(-abs(tstat12(2))))); 

disp('Negative log score - tstat');
disp(sprintf('ADAPT EIBV - ADAPT EIBV ref: %0.3g (%0.3g)',tstat12(3),normcdf(-abs(tstat12(3))))); 

gg=1;
