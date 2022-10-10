function droploc_seq_designsSYNT;

%
% Calculation of optimal sequential design
% for drop locations. 
%



%
%
% Set domain and grid settings
[delta,min_east,max_east,min_north,max_north,domain,n1,n2,n,nv,vx,vy,vel]=setparSYNT;

% Set list of candidate designs
% Set design drop location and operation time
[opstartVec,opstartTime]=setdesignsSYNT;


%% Prior model
[sig,nu,ksi,mubeta,Sigmabeta]=sethyperparSYNT;

%% Covariate inputs
[H,ids]=findcovariatesSYNT(domain,nv);
figure(6); clf; imagesc(vec2mat(H(:,2),n1)'); colorbar;
axis xy; axis image;
% Form H_{x,D}
HSigma=H*Sigmabeta;

% Spatial covariance matrix
% Embed this on the torus for FFT calculations
% Construct toeplitz (see separate file below)
[cm,cmf]=toep_circ(sig,nu,delta,nv,ksi);


%% SIMULATE DATA
%%%%%%%%%%%%
zz=sample0(cmf);
bet=mubeta+chol(Sigmabeta)'*randn(size(mubeta,1),1);
xx=H*bet+zz;
yGEN=getdata(xx,[1:size(zz,1)]');

figure(2);
clf;
imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(yGEN,n1)');
title('Simulated truth');
colorbar;
colormap(gray)
axis xy;
axis image;
xlabel('East (m)'); ylabel('North (m)');
%print -deps xtruth

figure(31); 
clf; 
hold;
for ii=1:size(opstartVec,2),
    [nodes,nodes1,nodes2]=findnodes(opstartVec(:,ii),opstartTime(ii),vel,domain,nv);
    %plot(min_east+delta*nodes1,min_north+delta*nodes2,'.k');
    plot(nodes1,nodes2,'.k');
    title('Sampling designs');
    %text(min_east+delta*nodes1(1,1)-29,-60+min_north+delta*nodes2(1,1),sprintf('%d',ii));
    text(nodes1(1,1)-2,nodes2(1,1)-4,sprintf('%d',ii));
  %axis([min_east max_east min_north max_north]);
    %axis xy;
    axis image;
    %axis([0 60 0 90]);
    axis([0 150 0 225])
    %axis([min_east max_east min_north max_north]);
    xlabel('Easting '); ylabel('Northing');
end;
figure(32);
clf;
%imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(xx,n1)');
imagesc(vec2mat(xx,n1)');
title('Linear predictor');
colorbar;
%colormap(gray)
axis xy;
axis image;
xlabel('Easting'); ylabel('Northing');

figure(33);
clf;
%imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(yGEN,n1)');
imagesc(vec2mat(yGEN,n1)');
title('Simulated truth');
colorbar;
%colormap(gray)
axis xy;
axis image;
xlabel('Easting'); ylabel('Northing ');

figure(1); 
clf;
hold;

%% STAGE 1 OF DESIGN

disp('STAGE 1');
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
plot(0,IBVres(1,1),'ko');
hold;
ylabel('IBV');
xlabel('Sampling stage');

EIBVvec=zeros(length(opstartTime),1);
for i=1:length(opstartTime),
    [nodesi,nodes1i,nodes2i]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
    [EIBV]=corals_comp_desNEW(H,Sigmabeta,cmf,cm,mubeta,nodesi,nodes1i,nodes2i);
    EIBVvec(i,1)=EIBV;
end;
%figure(1);
%print -deps Statsurvey

figure(8);
plot(1*ones(length(opstartTime),1),EIBVvec,'rx');
ylabel('IBV ');
xlabel('Sampling stage');

% Gather data according to selected design 
iimin=min(find(EIBVvec<min(EIBVvec)+0.01));
[nodes,nodes1,nodes2]=findnodes(opstartVec(:,iimin),opstartTime(iimin),vel,domain,nv);
Sel_line(1)=iimin;

% Get class labels of selected deployment
y=yGEN(nodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmean=zeros(n,1); % Allocate mean
srs=zeros(n,1); % Allocate variance
    
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
plot([0:1]',IBVres(1:2,1),'k');
ylabel('IBV');
xlabel('Sampling stage');
figure(19); 
clf;
subplot(2,1,1), plot(1,IBVval,'o'); ylabel('IBV'); hold;
MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
figure(19); 
subplot(2,1,2), plot(1,MCrate,'o'); ylabel('MCrate'); hold;

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
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],marge');
    title('Marginal predictions');
    colorbar;
    %colormap(gray)
    axis xy;
    axis image;
    xlabel('Easting (m)'); ylabel('Northing (m)');
    %print -depsc m1
    
    figure(14)
    clf;
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],margstd');
    title('Marginal predicted standard deviations');
    colorbar;
    %colormap(gray);
    axis xy;
    axis image;
    xlabel('Easting (m)'); ylabel('Northing (m)');
    %print -depsc v1
    
    probspred=Probseval(margmean,srs,nodes,yGEN);
    figure(15)
    clf;
    %imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(probspred,n1)');
    imagesc(vec2mat(probspred,n1)')
    title('Marginal predicted probabilities');
    colorbar;
    %colormap(gray);
    axis xy;
    axis image;
    xlabel('Easting (m)'); ylabel('Northing (m)');

    figure(5)
    clf;
    plot(min_east+delta*nodes1,min_north+delta*nodes2,'.k');
    title('Measurement locations');
    %axis([min_east max_east min_north max_north]);
    hold;
end;

desSAMPLED=zeros(length(opstartTime),1);
desSAMPLED(iimin,1)=1;


GGG=2;
%% STAGE 2, 3,...
%for ss=1:length(opstartTime)-1,
for ss=1:GGG,
    
    disp(sprintf('STAGE %d',ss+1));
 
    indxx=find(desSAMPLED==1);
    opstartVecOLD=opstartVec(:,indxx);
    opstartTimeOLD=opstartTime(indxx);
    nodesOLD=nodes;
    nodes1OLD=nodes1;
    nodes2OLD=nodes2;

    figure(1); 
    clf;
    hold;

    EIBVvec=zeros(length(opstartTime),1);
    for i=1:length(opstartTime),
        if (desSAMPLED(i,1)==1)
            EIBV=10000000*100000;
        else
           [nodesNEWi,nodes1NEWi,nodes2NEWi]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
            [EIBVi]=corals_comp_SEQdesNEW(H,Sigmabeta,cmf,cm,mOLD,srsOLD,RiOLD,nodesNEWi,nodes1NEWi,nodes2NEWi,nodesOLD,nodes1OLD,nodes2OLD);
            EIBV=min(EIBVi,IBVres(ss+1,1));
            figure(8);
            plot(ss+1,EIBV,'rx');
            ylabel('IBV');
            xlabel('Sampling stage');
        end;
        EIBVvec(i,1)=EIBV;
    end;
 
    % Gather data according to selected design 
    iimin2=min(find(EIBVvec<min(EIBVvec)+0.01));
    iimin2
    Sel_line(ss+1)=iimin2;

    desSAMPLED(iimin2,1)=1;
    
    [nodesNEW,nodes1NEW,nodes2NEW]=findnodes(opstartVec(:,iimin2),opstartTime(iimin2),vel,domain,nv);

    yNEW=yGEN(nodesNEW,1);
    yn=[y;yNEW];
    y=yn;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nodes=[nodesOLD;nodesNEW];
    nodes1=[nodes1OLD;nodes1NEW];
    nodes2=[nodes2OLD;nodes2NEW];

    margmean=zeros(n,1); % Allocate mean
    srs=zeros(n,1); % Allocate variance

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
        imagesc([min_east:delta:max_east],[min_north:delta:max_north],marge');
        title('Marginal predictions');
        colorbar;
        %colormap(gray)
        axis xy;
        axis image;
        xlabel('Easting (m)'); ylabel('Northing (m)');
        %print -depsc m2

        figure(14)
        clf;
        imagesc([min_east:delta:max_east],[min_north:delta:max_north],margstd');
        title('Marginal predicted standard deviations');
        colorbar;
        %colormap(gray);
        axis xy;
        axis image;
        xlabel('Easting (m)'); ylabel('Northing (m)');
        %print -depsc v2

        probspred=Probseval(margmean,srs,nodes,yGEN);
        figure(15)
        clf;
        imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(probspred,n1)');
        title('Marginal predicted probabilities');
        colorbar;
        %colormap(gray);
        axis xy;
        axis image;
        xlabel('Easting (m)'); ylabel('Northing (m)');

        % Plot the predicted latent state
%         figure(23)
%         clf;
%         imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(margmeanOLD,n1)'-marge');
%         title('Change in predictions');
%         colorbar;
%         colormap(gray)
%         axis xy;
%         axis image;
%         xlabel('East (m)'); ylabel('North (m)');
%         %print -deps m2
% 
%         figure(24)
%         clf;
%         imagesc([min_east:delta:max_east],[min_north:delta:max_north],sqrt(vec2mat(srsOLD,n1)')-margstd');
%         title('Reduction in predicted standard deviations');
%         colorbar;
%         colormap(gray);
%         axis xy;
%         axis image;
%         xlabel('East (m)'); ylabel('North (m)');
%         print -deps v2
% 
        figure(5)
        plot(min_east+delta*nodes1NEW,min_north+delta*nodes2NEW,'.k');
        title('Measurement locations');
        axis([min_east max_east min_north max_north]);

    end;

    IBVval=IBVeval(margmean,srs,nodes);
    disp(sprintf('Current IBV is %4.2f',IBVval));
    IBVres(ss+2,1)=IBVval;
    figure(8);
    plot([0:ss+1]',IBVres(1:ss+2,1),'k');
    ylabel('IBV');
    xlabel('Sampling stage');

    figure(19); 
    subplot(3,1,1), plot(ss+1,IBVval,'o'); ylabel('IBV');
    MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
    subplot(3,1,2), plot(ss+1,MCrate,'o'); ylabel('MCrate');
    LS=logscore(margmean,srs,nodes,yGEN);
    subplot(3,1,3), plot(ss+1,LS,'o'); ylabel('Log score');

    gg=1;
    

end

disp('Scores:');
[IBVval MCrate LS]

disp('Selection order:');
Sel_line

figure(8);
ylabel('IBV realized');
xlabel('Sampling stage');
