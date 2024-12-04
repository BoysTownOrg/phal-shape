function [ifig,hfShape]=ProcessCanal5_JASA(ifig,imult,fnM,swNewPlot,...
  dirESS,swSave,fnShape,swSimpleTitle)
% ProcessCanal5_JASA, for code sharing simplified
% Revised from ProcessCanal5, 10/1/21. 
% Load MAT scan file with grid size of 0.25 mm of
% canal stalk rotated to approximate lie vertically along z axis.
% Iteratively calculate center line with 1 mm steps along z axis.
% Calculate modeled ellipse area and eccentricity along the center line.
% Calculate curvature and torsion of center line embedded in 3D space.

swVerbose=0; % 1 to debug after setting iterEx to desired value, usually to
if swVerbose
  % iterEC from a previous run for same MAT file. Use different values for
  % different ear tests!
  iterEx=5; % set iterEC<=Niter for sub-plots for swVerbose==1 to debug
  jjEx=10;
end
Niter=10;
disp(fnM);
load(fnM,'mVxec','delZm','NStep','Nimx','Nimy','Xm','Ym','Zm','az','el',...
  'STLfile','TriangulationV','ss','iis');
if delZm==0.25
  sdelZm='Dp25'; % default grid size in 2D slices; step size between slices =imult*delZm
elseif delZm==0.5
  sdelZm='Dp50'; % obsolete but still works if MAT file is consistent
elseif delZm==1
  sdelZm='D1'; % obsolete but still works if MAT file is consistent
else
  disp('Surprising value of delZm');
  return
end

if swNewPlot  % Repeated from end of do12
  ifig=PlotTriangulation(ifig,TriangulationV,az,el,ss,iis);
  ifig=Plot12(Xm,Ym,Zm,mVxec,ss,az,el,delZm,ifig);
end

izec=NStep:-1:1; % z coordinate of canal in voxels
if contains(fnM,'Dp25','IgnoreCase',true)
  %  imult=2; % select, 1,2,3,4,... default 2
  if imult==1
    msD=4;
  else
    msD=5;
  end
  izec=izec(1:imult:end); % z coordinate of canal in pixels
  NStep=length(izec);
  % delZmAx applies to Axial length along s or z
  delZmAx=imult*delZm; % e.g., 4 steps at 0.25 mm is axial step of 1 mm
end
if NStep-4<1
  error('NStep too small to fit canal model');
end

mVxec=uint8(mVxec);
jfig=101;
kfig=201;
co=colororder;
if strcmp(sdelZm,'Dp25')
  edgeAngle=(-pi:pi/6:pi); % plane of ellipse with 30 deg spacing (wedges)
  sNwedge='12';
else
  edgeAngle=(-pi:pi/4:pi); % plane of ellipse with 45 deg spacing (wedges)
  sNwedge='8';
end
rngEndNStep=[1,2,NStep-1,NStep];
zc=2:(NStep-1);
zcp=3:NStep;
zcm=1:(NStep-2);
zrngMain=3:(NStep-2);
NsDiff=length(zrngMain);
s=zeros(Niter,NStep,3); % init variables, discrete center-line curve
Ts=zeros(Niter,NStep,3); % Tangent vector to center line
ds=zeros(Niter,NStep,3); % discrete 2-point differences along s
sDiffRMS=zeros(1,Niter); % RMS distance between iterations along s curve
errCtr=zeros(Niter,NStep);
ResidualError=zeros(Niter,NStep);
cntOpen=zeros(Niter,NStep);
area=zeros(Niter,NStep); % store x,y center of ellipse
eccentricity=zeros(Niter,NStep); % store eccentricity of ellipse
iLastInterior=zeros(1,Niter);
for iter=1:Niter %#ok<FXUP> % for loop to calculate curved center axis
  if swVerbose
    disp(['iter = ',num2str(iter)]);
  end
  for jj=1:NStep
    izjj=izec(jj);
    if iter==1
      ctr=[1,1,izjj]; % not center in x y projection plane for initial iter
      [ecSlice,ecX,ecY,ecZ]=...
        obliqueslice(mVxec,ctr,[0,0,-1]);
    else
      ctr=siterSmooth(jj,:); % ctr voxel # on previous iteration
      % Orientation of obliqueslice: coordinates of output slice corner
      % [x(1,1),y(1,1),z(1,1)] are at the upper-left pixel location in the
      % image plane and [x(1,end),y(1,end),z(1,end)] are at the
      % upper-right pixel location in the imagei plane.
      % Here, x,y and z are the coordinates of the output slice.
      [ecSlice,ecX,ecY,ecZ]=myobliqueslice(mVxec,ctr,...
        squeeze(Ts(iter-1,jj,:))');
    end
    if swVerbose && iter==iterEx && ~isempty(ecSlice) && ...
        ~isempty(intersect(jj,jjEx)) % for AAS
      hfjj=Plotjj;
      disp(['jj = ',num2str(jj)]);
    end
    iSec=find(ecSlice==1); % points in ecSlice forming wall of ear canal
    if ~isempty(iSec)
      [iSecRow,iSecCol]=ind2sub(size(ecSlice),iSec);
      if swVerbose && iter==iterEx && ~isempty(intersect(jj,jjEx)) % for AAS
        PlotDebug1;
      else
        el1=do_fit_ellipse_Halir(iSecCol,iSecRow,0); % default
      end
      % model center [X0r,Y0r] bounded by array min (1) and max indices
      X0r=min(max(round(el1.X0),1),size(ecSlice,2));
      Y0r=min(max(round(el1.Y0),1),size(ecSlice,1));
      ecX0=ecX(Y0r,X0r); % s coordinates in original coor system
      ecY0=ecY(Y0r,X0r);
      ecZ0=ecZ(Y0r,X0r);
      s(iter,jj,:)=[ecX0,ecY0,ecZ0]; % initial coordinates of center axis s
      errCtr(iter,jj)=el1.errCtr; % pixel distance between model and mean ctrs
      ResidualError(iter,jj)=el1.ResidualError; % least-squares fit error
      area(iter,jj)=el1.area;
      eccentricity(iter,jj)=el1.eccentricity;
      rctr2=[el1.Y0,el1.X0]; % ctr of 2D data in slice, flipped x-y
      zSec=(iSecRow-rctr2(1))+1i*(iSecCol-rctr2(2));
      zangle=angle(zSec); % interval [-pi,pi]
      Npts=histcounts(zangle,edgeAngle); % counts in each wedge
      empties=find(Npts==0); % empty if pts present in all wedges
      if ~isempty(empties)
        if swVerbose && iter==iterEx
          figure(jfig);
          jfig=jfig+1;
          hp=plot(iSecRow,iSecCol,'x',el1.Y0,el1.X0,'o');
          title(['iter = ',num2str(iter),', jj = ',num2str(jj)]);
        end
        cntOpen(iter,jj)=length(empties);
      end
%     else % error condition
%       error(['Empty iter=',int2str(iter),', jj=',int2str(jj)]);
%       keyboard
    end
  end % jj
  if iter==1 % relaxation smoothing in this if statement
    siter=squeeze(s(iter,:,:));
  elseif iter==2
    siter=1/1.5*squeeze(s(iter,:,:)+0.5*s(iter-1,:,:));
  elseif iter==3
    siter=1/1.75*squeeze(s(iter,:,:)+0.5*s(iter-1,:,:)...
      +0.25*s(iter-2,:,:));
  else % iter>=4
    siter=1/1.875*squeeze(s(iter,:,:)+0.5*s(iter-1,:,:)...
      +0.25*s(iter-2,:,:)+0.125*s(iter-3,:,:));
  end
  siter(rngEndNStep,:)=squeeze(s(1,rngEndNStep,:)); % constrain end values
  if iter>1
    % linear preferable to spline in interp1 here
    siterGrid=interp1(siter(:,3),siter,izec,'linear','extrap');
    df=diff(siterGrid,1,1)';
    t=cumsum([0,sqrt([1 1 1]*(df.*df))]); % cumulative length along s
    h=(t(end)-t(1))/(NStep-1);
    p5=1/(1+5*h^3/6);  %    p1=1/(1+h^3/6);
    pp5=csaps(izec,siterGrid',p5);
    siterSmooth=fnval(pp5,izec,'l')'; % evaluate at discrete z values
  else
    siterSmooth=siter;
  end
  ds(iter,1,:)=(siterSmooth(2,:)-siterSmooth(1,:)); % first ds as forward difference
  ds(iter,zc,:)=0.5*(siterSmooth(zcp,:)-siterSmooth(zcm,:)); % ds as center difference
  ds(iter,NStep,:)=(siterSmooth(NStep,:)-siterSmooth(NStep-1,:)); % last ds as backwards difference
  dsM=vecnorm(squeeze(ds(iter,:,:)),2,2); % magnitude of ds using Euclidean norm
  % Calc Ts before rounding s
  Ts(iter,:,:)=squeeze(ds(iter,:,:))./dsM; % unit tangent vector to s curve
  %s(iter,:,:)=round(siterSmooth); % set to closest voxel value
  siterSmooth=round(siterSmooth); % set to closest voxel value
  s(iter,:,:)=siterSmooth; % set to closest voxel value
  if iter>1
    sDiff=squeeze(s(iter,zrngMain,1:2)-s(iter-1,zrngMain,1:2));
    sDiffM=vecnorm(sDiff,2,2);
%    sDiffRMS(iter)=delZmAx*sqrt(sum(sDiffM.^2))/NsDiff; % in mm
    sDiffRMS(iter)=delZmAx*sqrt(sum(sDiffM.^2))/length(zrngMain); % in mm
    if swVerbose
      disp(['Iteration ',int2str(iter),...
        ': RMS change in sDiff (mm): ',num2str(sDiffRMS(iter),'%5.2f')]);
    end
  end
  temp=find(cntOpen(iter,zrngMain),1,'first');
  if ~isempty(temp)
    iLastInterior(iter)=find(cntOpen(iter,zrngMain),1,'first');
  else
    iLastInterior(iter)=NStep-3; % 1 more than length(zrngMain), all non-empty
  end
end % iter
rngNI=4:Niter;
maxiLI=max(iLastInterior(rngNI));
imaxiLI=find(maxiLI==iLastInterior(rngNI));
if length(imaxiLI)==1
  iterEC=imaxiLI+3;
else % find least iter with smallest sDiffRMS for iter's in imaxiLI
  [~,iSmallest]=min(sDiffRMS(rngNI(imaxiLI)));
  iterEC=imaxiLI(iSmallest)+3;
end
disp(['Define s at best iteration number: ',int2str(iterEC)]);
areaEC=delZm^2*area(iterEC,:); % areaEC in sq mm
eccentricityEC=eccentricity(iterEC,:);

%PlotFitStatistics
iSliceEnd=PlotListCenterLine;
% Calculate discrete Frenet frame in voxel units
sEC=squeeze(s(iterEC,:,:));
dsEC=squeeze(ds(iterEC,1:(NStep-1),:));
dsECM=vecnorm(dsEC,2,2);
[T,N,B,Kappa,Torsion]=DiscreteFrenetFrame(dsEC,dsECM,NStep);
sCenterAll=[0,delZm*cumsum(dsECM')];
sCenter=sCenterAll(zrngMain); % range of Kappa and Torsion
dsECMmm=delZm*dsECM(zrngMain); % center length in mm
Kappa1=Kappa./dsECMmm; % curvature per unit length
Torsion1=Torsion./dsECMmm; % torsion per unit length
jjEntranceEC=find(cntOpen(iterEC,zrngMain),1,'first')+zrngMain(1)-1;
if isempty(jjEntranceEC)
  jjEntranceEC=zrngMain(end)+1;
end
% defines canal entrance on multiple measures of goodness
if ~isempty(iSliceEnd)
  jjEntranceEC=min(jjEntranceEC,iSliceEnd);
end
% Finally, check preceding value of areaEC for any entrance discontinuity
fracChangeArea=...
  abs((areaEC(jjEntranceEC)-areaEC(jjEntranceEC-1))/areaEC(jjEntranceEC-1));
if fracChangeArea>0.10
  % discontinuity-- >10% change in area over adjacent slices near entrance
  disp('Area discontinuity at determined entrance');
  jjEntranceEC=jjEntranceEC-1;
end
dzM=delZmAx*(jjEntranceEC-1); % magnitude of length in mm along z axis
sEntrancemm=sCenterAll(jjEntranceEC);

sLattice=0:delZmAx:(delZmAx*floor(sCenterAll(end)/delZmAx)); % no extrapolation on lattice
iElbowLattice=2; % skip initial lattice value 0 so start at 0.5 mm
% Find index of sCenterAll closed to sLattice(iElbowLattice)
[~,iElbowCenterAll]=min(abs(sCenterAll-sLattice(iElbowLattice)));
disp(['Total z-length (mm) to entrance = ',num2str(dzM,'%4.1f')]);
disp(['Total s-length (mm) to entrance = ',...
  num2str(sEntrancemm,'%4.1f')]);
sLength=sEntrancemm-sCenterAll(iElbowCenterAll);
ssLength=num2str(sLength,'%4.1f');
disp(['s-length (mm) from near-TM to entrance = ',ssLength]);

if swVerbose % enable for debugging shape code
  DebugShape(NStep,sEC,iterEC,T,B,N,Kappa1,Torsion1,areaEC,eccentricityEC);
end
hfShape=figure(ifig);
%set(hfShape,'Position',[100,100,520,624],'DefaultLineMarkerSize',msD);
set(hfShape,'Position',[100,100,560,640],'DefaultLineMarkerSize',msD);
ListPlotResults_JASA(delZmAx,sCenterAll,areaEC,jjEntranceEC,sEntrancemm,...
  iElbowCenterAll,ssLength,ss,eccentricityEC,sCenter,Kappa1,Torsion1,...
  swSimpleTitle);
% Use rounding approach here to preserve more data points
% Cubic spline interpolate individual-ear values of area, eccentricity,
% curvature, torsion to their values along the curved centerline with
% uniform spacing, i.e., on sLattice.
% Determine entrance location of each ear on sLattice for use in
% subsequent group averaging.
[~,iEntrance0]=min(abs(sLattice-sEntrancemm));
if iEntrance0<iElbowLattice
  disp('Error Lattice 1--only 0 or 1 space lattice points');
  areaShape0=[];
  eccentricity0=[];
  iLatticeRng=[];
  iLatticeRngKT=[];
else % default
  iLatticeRng=iElbowLattice:iEntrance0;
  ppShape=spline(sCenterAll,areaEC);
  areaShape0=ppval(ppShape,sLattice(iLatticeRng)); % sample area on lattice
  ppE=spline(sCenterAll,eccentricityEC);
  eccentricity0=ppval(ppE,sLattice(iLatticeRng));
  % find range of Kappa and Torsion on sLattice
  % if (iEntrance0-2)<(iElbowLattice+2)
  %   disp('Error Lattice 2--too few lattice points for curv. & torsion');
  %   return
  % else
  %   iLatticeRngKT=(iElbowLattice+2):(iEntrance0);
  % end
  if (iEntrance0-2)<(iElbowLattice+2)
    disp('Error Lattice 2');
    return
  end
  iLatticeRngKT=(iElbowLattice+2):(iEntrance0-2);
  ppK=spline(sCenter,Kappa1);
  ppT=spline(sCenter,Torsion1);
  % Calculate Kappa0 and Torsion0 on lattice for swSave==1 after plotting
  Kappa0=ppval(ppK,sLattice(iLatticeRngKT)); % in inverse length
  Torsion0=ppval(ppT,sLattice(iLatticeRngKT)); % in inverse length
end

if swSave
  % Set variables used in sound&shape comparisons and other uses not
  % relevant here to empty before saving. This code makes no use of sound
  % data.
  sLatticeSound=[];
  areaSound=[];
  bigecSlice=[];
  bigecX=[];
  bigecY=[];
  bigecZ=[];
  bigel1=[];
  bigCtr=[];
  save(fullfile([dirESS,filesep 'Shape' filesep 'EarResults_JASA'],fnShape),'fnM',...
    'fnShape','imult','delZm','sCenterAll','areaEC','eccentricityEC',...
    'sLattice','areaShape0','eccentricity0','iEntrance0','iLatticeRng',...
    'iLatticeRngKT','sLength','sEntrancemm','sLattice','dsECMmm',...
    'sCenter','Kappa0','Torsion0','Kappa1','Torsion1','sLatticeSound',...
    'iElbowLattice','areaSound','iElbowCenterAll','jjEntranceEC','Ts',...
    'mVxec','TriangulationV','az','el','bigecSlice','bigecX','bigecY',...
    'bigecZ','bigel1','bigCtr');
end

  function hfjj=Plotjj
    hfjj=figure(kfig);
    kfig=kfig+1;
    set(hfjj,'Position',[200,200,1050,350]);
    subplot(1,3,1),
    hfs=surf(ecX,ecY,ecZ,ecSlice,'EdgeColor','None',...
      'HandleVisibility','off');
    hold on;
    plot3(ecX(1,1),ecY(1,1),ecZ(1,1),'om','MarkerFaceColor','m');
    plot3(ecX(1,end),ecY(1,end),ecZ(1,end),'oc','MarkerFaceColor','c');
    title(['iter=',int2str(iter),' jj=',int2str(jj),...
      ': slice embedded in 3D']);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    daspect([1,1,1]);
%    view(gca,-37.5,50);
    view(gca,-60,-15);
    drawnow;

    subplot(1,3,2),
    imshow(ecSlice,[],'InitialMagnification',2000);
    size(ecSlice)
    drawnow; % need this line present after imshow
    hold on;
    plot(1,1,'ms','MarkerFaceColor','m');
    plot(size(ecSlice,2),1,'cs','MarkerFaceColor','c');
    title('2D slice');
    xlabel('Col-Y reverse');
    ylabel('Row-X');
    drawnow;
  end

  function PlotDebug1 % debugging code
    subplot(1,3,3),
    % rows and columns flipped in going from image to physical geom
    title('Data and ellipse model fit');
    el1=do_fit_ellipse_Halir(iSecCol,iSecRow,...
      ifig,[0,0,0.9],'o','x');
    X0r=min(max(round(el1.X0),1),size(ecSlice,2));
    Y0r=min(max(round(el1.Y0),1),size(ecSlice,1));
    ylim=get(gca,'YLim');
    ylim(2)=ylim(2)+7;
    ylim(1)=ylim(1)-2;
    set(gca,'YLim',ylim); % make room for text
    xlim=get(gca,'XLim');
    xlim(1)=xlim(1)-2;
    xlim(2)=xlim(2)+2;
    set(gca,'XLim',xlim); % make room for text
    text(el1.atxt(1,2),el1.atxt(2,2),['a=',num2str(el1.a*delZm,3),' mm'],...
      'HorizontalAlignment','center');
    text(el1.btxt(1,1),el1.btxt(2,1),['b=',num2str(el1.b*delZm,3),' mm'],...
      'HorizontalAlignment','center');
    sarea=num2str(el1.area*delZm^2,3);
    seccentricity=num2str(el1.eccentricity,2);
    text(xlim(1)+1,ylim(2),{['$$\textrm{Area=}',sarea,' \textrm{mm}^2$$'],...
      ['$$\textrm{Eccentricity=}',seccentricity,...
      '=\sqrt{1-(b/a)^2}$$']},'interpreter','latex',...
      'VerticalAlignment','bottom');

    subplot(1,3,1),
    plot3(ecX(Y0r,X0r),ecY(Y0r,X0r),ecZ(Y0r,X0r),'s',...
      'MarkerFaceColor','w');
    drawnow;
  end

  function iSliceEnd=PlotListCenterLine
    hfs=figure(ifig);
    ifig=ifig+1;
    set(hfs,'Position',[200,200,640,770],'DefaultLineMarkerSize',msD);
    sleg=cell(1,Niter);
    snum=[];
    if iterEC==1
      ix1=iterEC;
    else
      ix1=1+mod(iterEC-1,size(co,1));
    end
    for iter=1:Niter
      sleg{iter}=int2str(iter);
      s1=delZm*squeeze(s(iter,:,:)); % convert to mm
      if iter<=7
        smk='o';
      elseif iter<=14
        smk='s';
      else
        smk='d';
      end
      hp=plot3(s1(:,1),s1(:,2),s1(:,3),'Marker',smk);
      if iter==iterEC
        hp.MarkerFaceColor=co(ix1,:);
      end
      snum=[snum,'    s',int2str(iter),'      ','      '];
      hold on;
    end
    xlabel('x  (mm)');
    ylabel('y  (mm)');
    zlabel('z  (mm)');
    daspect([1,1,1]);
    legend(sleg,'AutoUpdate','off');
    s1=delZm*squeeze(s(iterEC,:,:)); %  plot on top in mm
    if iterEC<=7
      smk='o';
    elseif iterEC<=14
      smk='s';
    else
      smk='d';
    end
    plot3(s1(:,1),s1(:,2),s1(:,3),'Marker',smk,'Color',co(ix1,:),...
      'MarkerFaceColor',co(ix1,:));
    title({ss,['Best iter=',int2str(iterEC)]});
    for iter=1:Niter
      if iter==1
        sval=round(squeeze(s(1,:,:)));
      else
        sval=[sval,round(squeeze(s(iter,:,:)))];
      end
    end
    if swVerbose
      disp(sval);
    end

    %eCtr=sqrt(2); % mm,<2/27/2022
    eCtr=3; % mm,>=2/27/2022
    eCnt=0;
    iCtr=find(errCtr(iterEC,:)>eCtr/delZm); % model/data ctrs too far apart
    iCnt=find(cntOpen(iterEC,:)>eCnt);
    xlim=[0,NStep];

    hmef=figure(ifig);
    ifig=ifig+1;
    swDefault=1;
    if swDefault
      set(hmef,'Position',[200,200,560,760],'DefaultLineMarkerSize',msD);
      subplot(3,1,1),
    else
      set(hmef,'Position',[300,300,560,560],'DefaultLineMarkerSize',msD);
      subplot(2,1,1),
    end
    plot(delZm*errCtr(iterEC,:),'-o');
    hold on;
    iSliceEnd=[];
    if ~isempty(iCtr)
      plot(iCtr,delZm*errCtr(iterEC,iCtr),'o','MarkerFaceColor','k');
      iSliceEnd=iCtr(1);
    end
    ylim=[0,max(3.25,delZm*ceil(max(errCtr(iterEC,:))))];
    set(gca,'XLim',xlim,'YLim',ylim);
    plot(xlim,[eCtr,eCtr],'k:');
    xlabel('Slice #');
    ylabel('Center Distance Error (mm)');
    title({[ss,' Model Error Functions'],...
      ['Best iteration: ',num2str(iterEC)]});
    text(0.5,ylim(2),...
      'Distance between model ellipse center & mean of data in slice',...
      'VerticalAlignment','top');

    if swDefault
      subplot(3,1,2),
      plot(delZm*ResidualError(iterEC,:),'-o');
      set(gca,'XLim',xlim,'YLim',[0,...
        delZm*ceil(max(ResidualError(iterEC,:)))]);
      xlabel('Slice #');
      ylabel('Residual Error (mm)');
      title('Residual error of model to fit ellipse');
      subplot(3,1,3),
    else
      subplot(2,1,2),
    end
    plot(cntOpen(iterEC,:),'-o');
    if ~isempty(iCnt)
      hold on;
      plot(iCnt,cntOpen(iterEC,iCnt),'o','MarkerFaceColor','k');
      iCntCanal=iCnt(iCnt>2); % omit first 2 points here near TM
      if ~isempty(iCntCanal)
        iSliceEnd=min(iCntCanal(1),iSliceEnd);
      end
    end
    plot(xlim,[eCnt,eCnt],'k:');
    set(gca,'XLim',xlim,'YLim',[-0.25,0.25+max(cntOpen(iterEC,:))]);
    xlabel('Slice #');
    ylabel('# empty sectors');
    title(['# of empty sectors (total ',sNwedge,' sectors) around ellipse center: ']);
    if ~isempty(iSliceEnd)
      iSliceEnd=iSliceEnd-1; % last slice # with good properties
    end
  end

  function ListPlotResults_JASA(delZmAx,sCenterAll,areaEC,jjEntranceEC,...
      sEntrancemm,iElbowCenterAll,ssLength,ss,eccentricityEC,...
      sCenter,Kappa1,Torsion1,swSimpleTitle)
%    ListPlotResults_JASA, 10/2/24, differs from the older ListPlotResults
%    function by replacing old Kappa and Torsion by Kappa0 and Torsion0,
    %-----plot final results----------------------------------------------
    xlim=[-0.5*delZmAx,ceil(sCenterAll(end)+0.5*delZmAx)];
    subplot(4,1,1),
    plot(sCenterAll,areaEC,'-bo'); % area in sq pixels in original grid
    sxlabel='Center distance (mm) from medial end';
    ylim=[0,max(100,max(areaEC))];
    hold on;
    if ~isempty(jjEntranceEC)
      smm2=repmat(sEntrancemm,1,2);
      plot(smm2,ylim,'k:');
      text(sEntrancemm,ylim(1),'Entrance',...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    end
    if ~isempty(iElbowCenterAll)
      sElbow=sCenterAll(iElbowCenterAll);
      plot(repmat(sElbow,1,2),ylim,'k:');
      text(sElbow,ylim(2),'Near TM',...
        'HorizontalAlignment','left','VerticalAlignment','top');
    end
    set(gca,'XLim',xlim,'YLim',ylim);
    xlabel(sxlabel);
    ylabel('Area  (mm^2)');
    if swSimpleTitle==0
      title({ss,['Canal length: ',ssLength,' mm']});
    else
      title(['Canal length: ',ssLength,' mm']);
    end

    subplot(4,1,2),
    plot(sCenterAll,eccentricityEC,'-bo');
    xlabel(sxlabel);
    ylabel('Eccentricity');
    hold on;
    ylim=[0,1];
    if ~isempty(jjEntranceEC)
      plot(smm2,ylim,'k:');
      text(sEntrancemm,0,'Entrance',...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    end
    if ~isempty(iElbowCenterAll)
      plot(repmat(sElbow,1,2),ylim,'k:');
      text(sElbow,ylim(2),'Near TM',...
        'HorizontalAlignment','left','VerticalAlignment','top');
    end
    set(gca,'XLim',xlim,'YLim',ylim);

    subplot(4,1,3),
    plot(sCenter,Kappa1,'-bo');
    xlabel(sxlabel);
    ylabel('Curvature  (1/mm)');
    hold on;
    ylim=get(gca,'YLim');
    ylim=[0,ylim(2)];
    if ~isempty(jjEntranceEC)
      plot(smm2,ylim,'k:');
      text(sEntrancemm,0,'Entrance',...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    end
    if ~isempty(iElbowCenterAll)
      plot(repmat(sElbow,1,2),ylim,'k:');
    end
    set(gca,'XLim',xlim);

    subplot(4,1,4),
    plot(sCenter,Torsion1,'-bo');
    xlabel(sxlabel);
    ylabel('Torsion  (1/mm)');
    hold on;
    ylim=get(gca,'YLim');
    ylim=[0,ylim(2)];
    if ~isempty(jjEntranceEC)
      plot(smm2,ylim,'k:');
      text(sEntrancemm,ylim(2),'Entrance',...
        'HorizontalAlignment','right','VerticalAlignment','top');
    end
    if ~isempty(iElbowCenterAll)
      plot(repmat(sElbow,1,2),ylim,'k:');
    end
    set(gca,'XLim',xlim);
  end

  function wpad=Pad1(w)
    wpad=[0;0;w';0;0];
  end

  function wpad=Pad3(w,z3)
    wpad=[z3;z3;squeeze(w);z3;z3];
  end

  function DebugShape(NStep,sEC,iterEC,T,B,N,Kappa1,Torsion1,areaEC,eccentricityEC)
    disp(['Results at best iteration number: ',int2str(iterEC)]);
    z3=[0,0,0];
    Tpad=Pad3(T,z3);
    Bpad=Pad3(B,z3);
    Npad=Pad3(N,z3);
    Kappapad=Pad1(Kappa1');
    Torsionpad=Pad1(Torsion1');
    tnames={'Step','sEC','T','B','N','Curvature','Torsion','Area',...
      'Eccentricity'};
    tbl=table((1:NStep)',sEC,Tpad,Bpad,Npad,Kappapad,Torsionpad,...
      areaEC',eccentricityEC','VariableNames',tnames)
  end

end
