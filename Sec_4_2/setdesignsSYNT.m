function [opstartVec,opstartTime]=setdesigns;

%
% Set the possible design drop locations (on the domain)
% and operation times (sec)
%


%opstartVec=[3300 7200;4000 8600;4400 7140;4400 8000;4200 8600;4200 8100;3600 7200;3700 7600;3800 7900;...
%    4300 7800;3900 7200;3400 6900;3400 8200;...
%    4000 7000;4300 7840;4000 7200]';
%opstartVec=[3900 7640;4300 7020;3340 6905;4190 7610;...
%    3798 7244;3501 7580;3852 8005;3380 8390;4350 8390]';%...
%    4150 7720;3900 7200;3400 6900;3400 8200;...
opstartVec=6850*ones(2,13);
opstartVec(1,1:13)=[3200:100:4400];

opstartTime=1500*ones(1,size(opstartVec,2));
%opstartTime=[300 1200 500 1000 1000 500 300 300 500 300 300 ...
%    200 200 200 200 200];
