function corals_compare_several_designs;

figure(1); 
clf;
hold;

opstartVec=[2200 7300;2300 6700;2300 7000;2400 6900;2800 6700;...
    3000 7200;2900 6900;3100 6700;3200 6900;...
    3500 6800;3500 7200]';
opstartTime=[500 1000 1000 500 1000 ...
    500 500 1000 500 ...
    500 500];
EIBVvec=zeros(1,length(opstartTime));
for i=1:length(opstartTime),
    [EIBV]=corals_comp_des(opstartVec(:,i),opstartTime(i));
    EIBVvec(i)=EIBV;
end;
