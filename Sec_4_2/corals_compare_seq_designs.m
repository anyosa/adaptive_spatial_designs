function corals_compare_seq_designs;

%
% Calculation of optimal sequential design
%
% 


%
%
% 
[delta,min_east,max_east,min_north,max_north,domain,n1,n2,n,nv,vx,vy,vel]=setpar;

%% Prior model
[sig,nu,ksi,mubeta,Sigmabeta]=sethyperpar;

%% Covariate inputs
H=findcovariates(domain,nv);
figure(6); clf; imagesc(vec2mat(H(:,2),n1)); colorbar;
% Form H_{x,D}
HSigma=H*Sigmabeta;

% Spatial covariance matrix
% Embed this on the torus for FFT calculations
% Construct toeplitz (see separate file below)
[cm,cmf]=toep_circ(sig,nu,delta,nv,ksi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATE DATA
%%%%%%%%%%%%
zz=sample0(cmf);
bet=mubeta+chol(Sigmabeta)'*randn(size(mubeta,1),1);
xx=H*bet+zz;

figure(2);
clf;
imagesc([min_east:delta:max_east],[min_north:delta:max_north],vec2mat(xx,n1)');
title('Simulated truth');
colorbar;
colormap(gray)
axis xy;
axis image;
xlabel('East (m)'); ylabel('North (m)');
print -deps xtruth

figure(1); 
clf;
hold;
%opstartVec=[2200 7300;2300 6700]';
opstartVec=[2200 7300;2300 6700;2600 7300;...
    2900 7000;3100 6700;...
    3400 7000]';
opstartTime=[500 500 500 500 500 500];
%opstartTime=[500 500];

%% STAGE 1
EIBVvec=zeros(1,length(opstartTime));
for i=1:length(opstartTime),
    [EIBV]=corals_comp_des(opstartVec(:,i),opstartTime(i));
    EIBVvec(i)=EIBV;
end;
figure(1);
%print -deps Statsurvey

% Gather data according to selected design 
iimax=find(EIBVvec==max(EIBVvec));
[nodes,nodes1,nodes2]=findnodes(opstartVec(:,iimax),opstartTime(iimax),vel,domain,nv);

eta=xx(nodes,1); %H(nodes,:)*bet+xx(nodes,1);
y=0;
for i=1:size(nodes,1),
    y(i,1)=(rand<(1/(1+exp(-eta(i,1)))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmean=zeros(n,1); % Allocate mean
srs=zeros(n,1); % Allocate variance
    
% Gaussian approximation (see separate file below)
[P,z,u,R,mode,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
margmean=mat2vec(mode);
% Arrange mean value on the spatial grid
marge=mode;

% Arrange standard deviation on the spatial grid
figure(11); 
clf;
plot(srs);
srs=max(0,srs);
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
    print -deps m1
    
    figure(4)
    clf;
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],margstd');
    title('Marginal predicted standard deviations');
    colorbar;
    colormap(gray);
    axis xy;
    axis image;
    xlabel('East (m)'); ylabel('North (m)');
    print -deps v1
    figure(5)
    clf;
    plot(min_east+delta*nodes1,min_north+delta*nodes2,'.k');
    title('Measurement locations');
    axis([min_east max_east min_north max_north]);
    hold;
end;

opstartVecOLD=opstartVec(:,iimax);
opstartTimeOLD=opstartTime(iimax);
nodesOLD=nodes;
nodes1OLD=nodes1;
nodes2OLD=nodes2;


%% STAGE 2

EIBVvec=zeros(1,length(opstartTime));
for i=1:length(opstartTime),
    if (i==iimax)
        EIBV=0;
    else
        [EIBV]=corals_comp_SEQdes(margmean,opstartVec(:,i),opstartTime(i),opstartVecOLD,opstartTimeOLD);
    end;
    EIBVvec(i)=EIBV;
end;

% Gather data according to selected design 
iimax2=find(EIBVvec==max(EIBVvec));
iimax2
[nodesNEW,nodes1NEW,nodes2NEW]=findnodes(opstartVec(:,iimax2),opstartTime(iimax2),vel,domain,nv);

etaOLD=eta;
etaNEW=xx(nodesNEW,1); %H(nodes,:)*bet+xx(nodes,1);
eta=[etaOLD;etaNEW];
for i=size(etaOLD,1)+1:size(eta,1),
    y(i,1)=(rand<(1/(1+exp(-eta(i,1)))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes=[nodesNEW;nodesOLD];
nodes1=[nodes1NEW;nodes1OLD];
nodes2=[nodes2NEW;nodes2OLD];

mmean=zeros(n,1); % Allocate mean
srs=zeros(n,1); % Allocate variance
    
% Gaussian approximation (see separate file below)
[P,z,u,R,mode,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
margmean=mat2vec(mode);
% Arrange mean value on the spatial grid
marge=mode;

% Arrange standard deviation on the spatial grid
figure(11); 
clf;
plot(srs);
srs=max(0,srs);
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
    colormap(gray)
    axis xy;
    axis image;
    xlabel('East (m)'); ylabel('North (m)');
    print -deps m2
    
    figure(14)
    clf;
    imagesc([min_east:delta:max_east],[min_north:delta:max_north],margstd');
    title('Marginal predicted standard deviations');
    colorbar;
    colormap(gray);
    axis xy;
    axis image;
    xlabel('East (m)'); ylabel('North (m)');
    print -deps v2

    figure(5)
    plot(min_east+delta*nodes1NEW,min_north+delta*nodes2NEW,'.k');
    title('Measurement locations');
    axis([min_east max_east min_north max_north]);
   
end;

gg=1;
