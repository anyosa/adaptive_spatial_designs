function droploc_seq_designsLIZARD;

%
% Calculation of optimal sequential design
% for drop locations.
% At a deployment, real data from Lizard Island are read.
%
%

%
%
% Set domain and grid settings
[delta,min_east,max_east,min_north,max_north,domain,n1,n2,n,nv,vx,vy,vel]=setparLIZARD;

% Set list of candidate designs
% Set design drop location and operation time
[opstartVec,opstartTime]=setdesignsLIZARD;
 

%% Prior model
[sig,nu,ksi,mubeta,Sigmabeta]=sethyperparLIZARD;

%% Covariate inputs
[H,ids]=findcovariatesLIZARD(domain,nv);
figure(6); clf; imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(H(:,2),n1)'); colorbar;
axis xy; axis image;
% Form H_{x,D}
HSigma=H*Sigmabeta;

% Spatial covariance matrix
% Embed this on the torus for FFT calculations
% Construct toeplitz (see separate file below)
[cm,cmf]=toep_circ(sig,nu,delta,nv,ksi);

%% REAL DATA
DATAMAT=csvread('merged_binary_features.csv',1,0);
yREAL=zeros(n,1);
y0=DATAMAT(:,3);
yREAL(ids,1)=y0;
ys1ind=find(yREAL==1);

figure(6);
hold;
plot(min_east+ys1ind-n1*floor(ys1ind/n1),min_north+ceil(ys1ind/n1),'ko');
xlabel('Easting (m)'); ylabel('Northing (m)');

% Variances of x for these parameters sigma and nu
LSigma=chol(Sigmabeta)';
HVV=H*LSigma;

% Variances of x for these parameters sigma and nu
muvec=H*mubeta;
srs=cm(1,1)*ones(n,1)+sum(HVV.^2,2);
IBVval=IBVeval(muvec,srs,0);
disp(sprintf('Current IBV is %4.2f',IBVval));
IBVres=zeros(length(opstartTime)+1,1);
IBVres(1,1)=IBVval;
figure(8);
clf;
plot(1,IBVres(1,1),'ko');
hold;
ylabel('IBV');
xlabel('Sampling stage');

figure(19); 
clf; 
subplot(1,3,1), plot(0,IBVval,'ko'); ylabel('IBV');
xlabel('Stage'); hold;
MCrate=MisClassRateeval(muvec,srs,1,yREAL);
figure(19); 
subplot(1,3,2), plot(0,MCrate,'ko'); ylabel('Integrated mis-classification probability');
xlabel('Stage'); hold;
logsc=logscore(muvec,srs,1,yREAL);
subplot(1,3,3), plot(0,logsc,'ko'); ylabel('Negative log score');
xlabel('Stage'); hold;


probspred=Probseval(muvec,srs,1,yREAL);
iiiind=min(find(probspred==max(probspred)));
probspred(iiiind,1)=1;
figure(15)
clf;
imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(probspred,n1)');
title('Marginal class prediction');
colorbar;
axis xy;
axis image;
xlabel('East (m)'); ylabel('North (m)');

figure(25)
clf;
imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(probspred,n1)');
title('Predicted presence probabilities. (Initially.)');
colorbar;
colormap(jet);
axis xy;
axis image;
xlabel('Easting (m)'); ylabel('Northing (m)');

figure(1); 
clf;
hold;

%% STAGE 1 OF DESIGN
disp('STAGE 1');
tic;
EIBVvec=zeros(length(opstartTime),1);
for i=1:length(opstartTime),
    [nodesi,nodes1i,nodes2i]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
    [EIBV]=corals_comp_desNEW(H,Sigmabeta,cmf,cm,mubeta,nodesi,nodes1i,nodes2i);
    EIBVvec(i,1)=EIBV;
    toc;
end;

% Gather data according to selected design 
iimin=find(EIBVvec<min(EIBVvec)+0.01);
[nodes,nodes1,nodes2]=findnodes(opstartVec(:,iimin),opstartTime(iimin),vel,domain,nv);
Sel_line(1)=iimin;

% Get class labels of selected deployment
y=yREAL(nodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian approximation (see separate file below)
[P,z,u,Ri,margmean,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
mOLD=margmean;
srsOLD=srs;
RiOLD=Ri;
% Arrange mean value on the spatial grid
marge=vec2mat(margmean,n1);

IBVval=IBVeval(margmean,srs,nodes);
disp(sprintf('Current IBV is %4.2f',IBVval));
IBVres(2,1)=IBVval;
figure(8);
plot([1:2]',IBVres(1:2,1));
ylabel('IBV');
xlabel('Sampling stage');

figure(19); 
subplot(1,3,1), plot(1,IBVval,'ko'); ylabel('IBV');
xlabel('Stage');
MCrate=MisClassRateeval(margmean,srs,nodes,yREAL);
figure(19); 
subplot(1,3,2), plot(1,MCrate,'ko'); ylabel('Integrated mis-classification probability');
xlabel('Stage');
logsc=logscore(margmean,srs,nodes,yREAL);
figure(19); 
subplot(1,3,3), plot(1,logsc,'ko'); ylabel('Negative log score');
xlabel('Stage');


% Arrange standard deviation on the spatial grid
figure(11); 
clf;
plot(srs);
margstd=sqrt(vec2mat(srs,n1));

%%%%%%%%%%%%%%%%%%
% Plotting of process
%%%%%%%
pl=1;
if (pl==1)

    % Plot the predicted latent state
    figure(3)
    clf;
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],marge');
    title('Marginal predictions');
    colorbar;
    colormap(gray)
    axis xy;
    axis image;
    xlabel('East (m)'); ylabel('North (m)');
    
    figure(4)
    clf;
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],margstd');
    title('Marginal predicted standard deviations');
    colorbar;
    colormap(gray);
    axis xy;
    axis image;
    xlabel('East (m)'); ylabel('North (m)');
        
    probspred=Probseval(margmean,srs,nodes,yREAL);
    figure(15)
    clf;
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(probspred,n1)');
    title('Marginal class prediction (end stage)');
    colorbar;
    colormap(gray);
    axis xy;
    axis image;
    xlabel('East (m)'); ylabel('North (m)');

    figure(5)
    clf;
    plot(min_east+delta*nodes1,min_north+delta*nodes2,'.k');
    title('Measurement locations');
    axis([min_east max_east min_north max_north]);
    hold;
end;

desSAMPLED=zeros(length(opstartTime),1);
desSAMPLED(iimin,1)=1;

figure(1); 
clf;
hold;

%GG=7;
GG=14;
%% STAGE 2, 3,...
for ss=1:GG-1,
    
    disp(sprintf('STAGE %d',ss+1));

    indxx=find(desSAMPLED==1);
    opstartVecOLD=opstartVec(:,indxx);
    opstartTimeOLD=opstartTime(indxx);
    nodesOLD=nodes;
    nodes1OLD=nodes1;
    nodes2OLD=nodes2;

    EIBVvec=zeros(length(opstartTime),1);
    for i=1:length(opstartTime),
       if (desSAMPLED(i,1)==1)
            EIBV=10000000*100000;
       else
            [nodesNEWi,nodes1NEWi,nodes2NEWi]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
            [EIBVi]=corals_comp_SEQdesNEW(H,Sigmabeta,cmf,cm,mOLD,srsOLD,RiOLD,nodesNEWi,nodes1NEWi,nodes2NEWi,nodesOLD,nodes1OLD,nodes2OLD);
            EIBV=min(EIBVi,IBVres(ss+1,1));
        end;
        EIBVvec(i,1)=EIBV;
    end;

    % Gather data according to selected design 
    iimin2=find(EIBVvec<min(EIBVvec)+0.01);
    iimin2
    desSAMPLED(iimin2,1)=1;
    [nodesNEW,nodes1NEW,nodes2NEW]=findnodes(opstartVec(:,iimin2),opstartTime(iimin2),vel,domain,nv);
    Sel_line(ss+1)=iimin2;    

    yNEW=yREAL(nodesNEW,1);
    yn=[y;yNEW];
    y=yn;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nodes=[nodesOLD;nodesNEW];
    nodes1=[nodes1OLD;nodes1NEW];
    nodes2=[nodes2OLD;nodes2NEW];
   
    % Gaussian approximation (see separate file below)
    [P,z,u,Ri,margmean,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
    mOLD=margmean;
    srsOLD=srs;
    RiOLD=Ri;
    % Arrange mean value on the spatial grid
    marge=vec2mat(margmean,n1);

    % Arrange standard deviation on the spatial grid
    figure(11); 
    clf;
    plot(srs);
    margstd=sqrt(vec2mat(srs,n1));

    %%%%%%%%%%%%%%%%%%
    % Plotting of process
    %%%%%%%
    pl=1;
    if (pl==1)

        % Plot the predicted latent state
        figure(13)
        clf;
        imagesc(marge');
        title('Marginal predictions');
        colorbar;
        colormap(gray)
        axis xy;
        axis image;
        xlabel('East (m)'); ylabel('North (m)');

        figure(14)
        clf;
        %imagesc([min_east:delta:max_east],[min_north:delta:max_north],margstd');
        imagesc(margstd');
        title('Marginal predicted standard deviations');
        colorbar;
        colormap(gray);
        axis xy;
        axis image;
        xlabel('East (m)'); ylabel('North (m)');

        probspred=Probseval(margmean,srs,nodes,yREAL);
        figure(15)
        clf;
        imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(probspred,n1)');
        title('Predicted presence probabilities. (Last stage.)');
        colorbar;
        colormap(jet);
        axis xy;
        axis image;
        xlabel('Easting (m)'); ylabel('Northing (m)');

        figure(5)
        plot(min_east+delta*nodes1NEW,min_north+delta*nodes2NEW,'.k');
        title('Measurement locations');
        axis([min_east max_east min_north max_north]);

    end;

    
    IBVval=IBVeval(margmean,srs,nodes);
    disp(sprintf('Current IBV is %4.2f',IBVval));
    IBVres(ss+2,1)=IBVval;
    figure(8);
    plot([0:ss+1]',IBVres(1:ss+2,1));
    ylabel('IBV');
    xlabel('Sampling stage');

    figure(19); 
    subplot(1,3,1), plot(ss+1,IBVval,'ko'); ylabel('IBV');
    xlabel('Stage');
    MCrate=MisClassRateeval(margmean,srs,nodes,yREAL);
    figure(19); 
    subplot(1,3,2), plot(ss+1,MCrate,'ko'); ylabel('Integrated mis-classification probability');
    xlabel('Stage');
    logsc=logscore(margmean,srs,nodes,yREAL);
    figure(19); 
    subplot(1,3,3), plot(ss+1,logsc,'ko'); ylabel('Negative log score');
    xlabel('Stage');
    
end

disp('Scores:');
[IBVval MCrate logsc]

disp('Selection order:');
Sel_line

