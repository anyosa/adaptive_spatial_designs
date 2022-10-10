function [scorevecADAPT,scorevecP,scorevecSPAT]=ReplicateSynt;

% Number of replicates
B=100;

%% Case specifics

% Set domain and grid settings
[delta,min_east,max_east,min_north,max_north,domain,n1,n2,n,nv,vx,vy,vel]=setparSYNT;
% Set list of candidate designs
% Set design drop location and operation time
[opstartVec,opstartTime]=setdesignsSYNT;

%% REPLICATE RUNS
%%%%%%%%%%%%
for bb=1:B,

    % Prior model
    [sig,nu,ksi,mubeta,Sigmabeta,gamma]=sethyperparSYNT;

    %% Covariate inputs
    [H,ids]=findcovariatesSYNT(domain,nv);
    % Form H_{x,D}
    HSigma=H*Sigmabeta;

    %% Generate data
    [cm,cmf]=toep_circ(sig,nu,delta,nv,ksi);

    zz=sample0(cmf);
    bet=mubeta+chol(Sigmabeta)'*randn(size(mubeta,1),1);
    xx=H*bet+zz;
    rr=0; for ii=1:size(H,1), rr=rr+H(ii,:)*Sigmabeta*H(ii,:)'; end; 
    sig2reg=rr/size(H,1)
    sig2spat=sig^2

    %% IF DATA GENERATION
    yGEN=getdata(xx,[1:size(zz,1)]');
    fid=fopen(sprintf('ysimL%d.txt',bb),'w');
    fprintf(fid,'%d\n',yGEN);
    fclose(fid);
    %% IF DATA READING FROM FILES
    %fid=fopen(sprintf('ysimL%d.txt',bb),'r');
    %yGEN=fscanf(fid,'%d\n');
    %fclose(fid);

    
    %% Adaptive IEBV strategy
    
    % STAGE 1 OF DESIGN
    disp('STAGE 1');
    LSigma=chol(Sigmabeta)';
    HVV=H*LSigma;

    % Variances of x for these parameters sigma and nu
    %muvec=H*mubeta;
    %srs=cm(1,1)*ones(n,1)+sum(HVV.^2,2);

    EIBVvec=zeros(length(opstartTime),1);
    for i=1:length(opstartTime),
        [nodesi,nodes1i,nodes2i]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
        [EIBV]=corals_comp_desNEW(H,Sigmabeta,cmf,cm,mubeta,nodesi,nodes1i,nodes2i);
        EIBVvec(i,1)=EIBV;
    end;

    % Gather data according to selected design 
    iimin=min(find(EIBVvec<min(EIBVvec)+0.01));
    [nodes,nodes1,nodes2]=findnodes(opstartVec(:,iimin),opstartTime(iimin),vel,domain,nv);
    Sel_line(1)=iimin;

    % Get class labels of selected deployment
    y=yGEN(nodes,1);

    [P,z,u,Ri,margmean,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
    mOLD=margmean;
    srsOLD=srs;
    RiOLD=Ri;

    % Arrange mean value on the spatial grid
    marge=vec2mat(margmean,n1);

    IBVval=IBVeval(margmean,srs,nodes);
    IBVres(1)=IBVval;
          LS=logscore(margmean,srs,nodes,yGEN);
          MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
  
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

        EIBVvec=zeros(length(opstartTime),1);
        for i=1:length(opstartTime),
            if (desSAMPLED(i,1)==1)
                EIBV=10000000*100000;
            else
               [nodesNEWi,nodes1NEWi,nodes2NEWi]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
                [EIBVi]=corals_comp_SEQdesNEW(H,Sigmabeta,cmf,cm,mOLD,srsOLD,RiOLD,nodesNEWi,nodes1NEWi,nodes2NEWi,nodesOLD,nodes1OLD,nodes2OLD);
                EIBV=min(EIBVi,IBVres(ss,1));
                figure(8);
                plot(ss+1,EIBV,'rx');
                ylabel('IBV');
                xlabel('Sampling stage');
            end;
            EIBVvec(i,1)=EIBV;
        end;

        % Gather data according to selected design 
        iimin2=min(find(EIBVvec<min(EIBVvec)+0.01));
        Sel_line(ss+1)=iimin2;
        desSAMPLED(iimin2,1)=1;
        [nodesNEW,nodes1NEW,nodes2NEW]=findnodes(opstartVec(:,iimin2),opstartTime(iimin2),vel,domain,nv);
        yNEW=yGEN(nodesNEW,1);
        yn=[y;yNEW];
        y=yn;
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
        margstd=sqrt(vec2mat(srs,n1));

        IBVval=IBVeval(margmean,srs,nodes);
        IBVres(ss+1,1)=IBVval;

        LS=logscore(margmean,srs,nodes,yGEN);
         MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
 
    end

    disp('Scores:');
    scorevecADAPT(bb,:)=[IBVval MCrate LS]

    disp('Selection order:');
    Sel_lineADAPT(bb,:)=Sel_line;

 
    %% SPATIALLY BALANCED strategy
    %DesInd=[2 9 3 10];
    DesInd=[1 6 13];
    %DesInd=[1 2 3 8 9];
    
    % STAGE 1 OF DESIGN
    disp('STAGE 1');
    LSigma=chol(Sigmabeta)';
    HVV=H*LSigma;

    % Variances of x for these parameters sigma and nu
    muvec=H*mubeta;
    srs=cm(1,1)*ones(n,1)+sum(HVV.^2,2);

    % Gather data according to selected design 
    iimin=DesInd(1);
    [nodes,nodes1,nodes2]=findnodes(opstartVec(:,iimin),opstartTime(iimin),vel,domain,nv);
    Sel_lineSB(1)=iimin;

    % Get class labels of selected deployment
    y=yGEN(nodes,1);

    [P,z,u,Ri,margmean,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
    mOLD=margmean;
    srsOLD=srs;
    RiOLD=Ri;

    % Arrange mean value on the spatial grid
    marge=vec2mat(margmean,n1);

    IBVval=IBVeval(margmean,srs,nodes);
   
          LS=logscore(margmean,srs,nodes,yGEN);
          MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
      
    desSAMPLED=zeros(length(opstartTime),1);
    desSAMPLED(iimin,1)=1;

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

        % Gather data according to selected design 
        iimin2=DesInd(ss+1);
        Sel_lineSB(ss+1)=iimin2;

        desSAMPLED(iimin2,1)=1;

        [nodesNEW,nodes1NEW,nodes2NEW]=findnodes(opstartVec(:,iimin2),opstartTime(iimin2),vel,domain,nv);

        yNEW=yGEN(nodesNEW,1);
        yn=[y;yNEW];
        y=yn;
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
        margstd=sqrt(vec2mat(srs,n1));

        IBVval=IBVeval(margmean,srs,nodes);

          LS=logscore(margmean,srs,nodes,yGEN);
         MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
 

    end

    disp('Scores:');
    scorevecSPAT(bb,:)=[IBVval MCrate LS]

    disp('Selection order:');
    Sel_lineSPAT(bb,:)=Sel_lineSB;


    %% P=0.5 Adaptive strategy
    
    % STAGE 1 OF DESIGN
    disp('STAGE 1');
    LSigma=chol(Sigmabeta)';
    HVV=H*LSigma;

    % Variances of x for these parameters sigma and nu
    muvec=H*mubeta;
    srs=cm(1,1)*ones(n,1)+sum(HVV.^2,2);
    
    IPvec=zeros(1,length(opstartTime));
    for i=1:length(opstartTime),
              [nodesNEWi,nodes1NEWi,nodes2NEWi]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
              avtransect=AveVarTransect(muvec,srs,nodesNEWi);
              IPvec(1,i)=avtransect;
    end;

    % Gather data according to selected design 
    iimin=min(find(IPvec>max(IPvec)-0.001));
    [nodes,nodes1,nodes2]=findnodes(opstartVec(:,iimin),opstartTime(iimin),vel,domain,nv);
    Sel_lineP(1)=iimin;

    % Get class labels of selected deployment
    y=yGEN(nodes,1);

    [P,z,u,Ri,margmean,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
    mOLD=margmean;
    srsOLD=srs;
    RiOLD=Ri;

    % Arrange mean value on the spatial grid
    marge=vec2mat(margmean,n1);

    IBVval=IBVeval(margmean,srs,nodes);
    
          LS=logscore(margmean,srs,nodes,yGEN);
         MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);
  
    desSAMPLED=zeros(length(opstartTime),1);
    desSAMPLED(iimin,1)=1;

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

        IPvec=zeros(1,length(opstartTime));
        for i=1:length(opstartTime),
            if (desSAMPLED(i,1)==1)
                avtransect=-10;
            else
                [nodesNEWi,nodes1NEWi,nodes2NEWi]=findnodes(opstartVec(:,i),opstartTime(i),vel,domain,nv);
                avtransect=AveVarTransect(mOLD,srsOLD,nodesNEWi);
            end;
            IPvec(i,1)=avtransect;
        end;

        % Gather data according to selected design 
        iimin2=min(find(IPvec>max(IPvec)-0.001));
        iimin2
        Sel_lineP(ss+1)=iimin2;

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

        [P,z,u,Ri,margmean,srs]=Beval(H,HSigma,mubeta,Sigmabeta,cmf,cm,nodes,nodes1,nodes2,y);
        mOLD=margmean;
        srsOLD=srs;
        RiOLD=Ri;

        % Arrange standard deviation on the spatial grid
         margstd=sqrt(vec2mat(srs,n1));

        IBVval=IBVeval(margmean,srs,nodes);

         LS=logscore(margmean,srs,nodes,yGEN);
        MCrate=MisClassRateeval(margmean,srs,nodes,yGEN);

    end;

            disp('Scores:');
        scorevecP(bb,:)=[IBVval MCrate LS]

           disp('Selection order:');
    Sel_linePP(bb,:)=Sel_lineP;

    disp(sprintf('Replicate %d of %d',bb,B));
    
    gggggg=1;
    
    % Print to file
    % Not needed if stored before and re-read
    fileID=fopen(sprintf('Replic%d.csv',bb),'w');
    fprintf(fileID,'%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n',...
        [sig Sigmabeta(1,1) nu Sigmabeta(2,1) scorevecADAPT(bb,1) scorevecADAPT(bb,2) scorevecADAPT(bb,3) scorevecSPAT(bb,1) scorevecSPAT(bb,2) scorevecSPAT(bb,3) scorevecP(bb,1) scorevecP(bb,2) scorevecP(bb,3)]);
    fclose(fileID);

    
end;

Sel_lineADAPT


Sel_lineSPAT

Sel_linePP

disp('Adaptive mean and std');
mean(scorevecADAPT)
std(scorevecADAPT)

% disp('Pre-scripted mean and std');
% mean(scorevecPRESCRIPTED)
% std(scorevecPRESCRIPTED)
% 
disp('Spatially balanced mean and std');
mean(scorevecSPAT)
std(scorevecSPAT)

disp('P-choice mean and std');
mean(scorevecP)
std(scorevecP)

figure(21); 
clf; 
for ii=1:3,
    subplot(1,3,ii), hist([scorevecADAPT(:,ii) scorevecP(:,ii) scorevecSPAT(:,ii)]); colormap(gray);
    legend('EIBV adapt','P adapt','Spat bal');
end;
subplot(1,3,1), sx=xlabel('EIBV'); set(sx,'FontSize',16);
sy=ylabel('Fraction'); set(sy,'FontSize',16); 
subplot(1,3,2), sx=xlabel('Int mis-classif probability'); set(sx,'FontSize',16);
subplot(1,3,3), sx=xlabel('Negative log score'); set(sx,'FontSize',16);
fggdsfs=1;

% Differences
tstatdiff(scorevecADAPT,scorevecP,scorevecSPAT)
ranksscores(scorevecADAPT,scorevecP,scorevecSPAT);

figure(29); 
clf; 
hist([mat2vec(Sel_lineADAPT) mat2vec(Sel_lineSPAT) mat2vec(Sel_linePP)]);
legend('EIBV','SPAT','PredV');
xlabel('Sampling line'); ylabel('Times chosen'); 

ggg=1;



