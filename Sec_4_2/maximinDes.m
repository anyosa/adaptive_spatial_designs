function DesInd=maximinDes(opstartVec);

% Answer is DesInd=[2 3 8 9 10]

Itt=0;
ddmaxSOFAR=0;
while (Itt<10000)
    kk=randsample(1:10,5,false);
    dd1=(opstartVec(1,kk)'*ones(1,5)-ones(5,1)*opstartVec(1,kk));
    dd2=(opstartVec(2,kk)'*ones(1,5)-ones(5,1)*opstartVec(2,kk));
    dd=sqrt(dd1.^2+dd2.^2);
    dd=10000000*eye(5)+dd;
    ddmm=min(min(dd))
    if (ddmm>ddmaxSOFAR)
        kkBEST=kk;
        ddmaxSOFAR=ddmm;
    end;
    Itt=Itt+1
    
end

DesInd=kkBEST
ddmaxSOFAR

