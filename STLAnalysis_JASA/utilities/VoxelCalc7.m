function [ifig,Vx,izMin,izMax,CmaxVoxel,CminVoxel,TriangulationV]=...
  VoxelCalc7(iis,ifig,delZm,Triangulation0,F,ss,az,el,zmmCut,iOver,iter)
% VoxelCalc7, 8/19/21, remove PeakDistance logic entirely. Deal differently
% with single peak case on both iterations by setting min z 20 mm below max
% z.
% VoxelCalc6, 6/29/21. Handle case of one extremum only at medial canal
% tip. Stalk is then rotated in do12 to vertically align the stalk. 
%   Reduce triangulation to canal stalk only, removing secondary peak of
%   superior conchal eminence.
%   Inputs:
%     iter, iteration #; ifig, current figure #
%     iis,
%     ifig, current figure #
%     delZm, z step size (in mm) for voxelization
%     Triangulation0, triangulation with vertices V0
%     F, connectivity list of triangles
%     ss, identifier for test ear
%     az, default azimuth; el, default elevation
%     zmmCut, if only 1 peak then voxels retained to zmmCut mm down
%   Outputs:
%     ifig, new current figure #
%     Vxec, binary voxel array of canal stalk for detailed calculations
%     Nstep, # of z axis steps; 
%     CmaxVoxel, max center
%     TriangulationV, new triangulation of canal stalk only after vol seg
swVerbose=0;
figbottom = 150;
V0=Triangulation0.Points;
SP0=bounds(V0,1);
V=V0+delZm-min(SP0); % Translate V0 re: its minimum value and delZm
clear V0 SP0
[SP,LP]=bounds(V,1); % LP is largest of all vertices V for all dimensions
%  easily interpret voxelized scan data.
TriangulationV=triangulation(F,V);
C=incenter(TriangulationV); % centers of all triangles
FN=faceNormal(TriangulationV); % normal vector to each triangle
[ifig,hf1]=PlotTriangulation(ifig,TriangulationV,az,el,ss,iis);

LPmax=delZm*ceil(max(LP)/delZm);
vxlmm=delZm:delZm:LPmax; % discrete each spatial dimension
Nvxl=length(vxlmm); % used below for spatial extent of voxelization
[X,Y] = meshgrid(vxlmm,vxlmm);
Z = gridtrimesh(F,V,X,Y); % non-MATLAB X,Y from meshgrid
[ZMAX,IMAX] = extrema2(Z); % local extrema from mesh, non-MATLAB
NMAX=length(ZMAX);
if NMAX>1
  ZMAX2=abs(ZMAX(2:NMAX)-ZMAX(1));
  a=find(ZMAX2>4.25,1,'first'); 
  if isempty(a)
    ZMAX=ZMAX(1);
    IMAX=IMAX(1);
    NMAX=1;
    disp('Only one peak');
  else
    ap1=a+1;
    if isempty(iOver) % default, ear-canal near TM has largest z value
      ZMAX=ZMAX([1,ap1]);
      IMAX=IMAX([1,ap1]);
    else % ear-canal near TM has 2nd largest z value, happens for child molds
      ZMAX=ZMAX([ap1,1]);
      IMAX=IMAX([ap1,1]);
    end
    NMAX=2;
  end
end
% use largest two peaks for ear canal and cavity unless only 1 peak
Pq=[X(IMAX),Y(IMAX),Z(IMAX)];
ID=nearestNeighbor(TriangulationV,Pq); %V(ID) % display nearest vertices
FA=vertexAttachments(TriangulationV,ID); % list of triangles
[CFA1zmax,i1max]=max(C(FA{1},3)); % identify what should be ear-canal peak (in mm)
Cmax=C(FA{1}(i1max),:); % Cmax is center of ear-canal triangle in mm with maximum z
CmaxVoxel=round(Cmax/delZm)-[0,0,1]; % transform distance (mm) to voxels
if NMAX>1
  [CFA2zmax,i2max]=max(C(FA{2},3)); % identify secondary peak (in mm)
  Cmin=C(FA{2}(i2max),:);
  CminVoxel=round(Cmin/delZm)+[0,0,1]; % transform distance (mm) to voxels (1 above min)
  figure(hf1);
  hold on; 
  lw=2;
  plot3([Cmin(1),Cmax(1)],[Cmin(2),Cmax(2)],[Cmin(3),Cmax(3)],...
    '-o','Color','k','MarkerSize',6,'LineWidth',lw);
else
  CminVoxel=CmaxVoxel;
  CminVoxel(3)=CminVoxel(3)-round(zmmCut/delZm); % zmmCut mm down if only 1 peak
end
izMin=CminVoxel(3);
% Cmin2Voxel corresponds to voxel coordinates of secondary peak if present
vxl=1:Nvxl; % in voxels
Nvxl3D=repmat(Nvxl,1,3);
FV.vertices=V/delZm; % FV is in voxel units
FV.faces=F;
% Convert mesh to voxelvolume
Vx=polygon2voxel(FV,Nvxl3D,'none'); % binary voxel array of TriangulationV

izMax=Nvxl;
CmaxVoxel3=CmaxVoxel(3); %z-axis range in voxels
swBreak=0;
while ~swBreak
  izStep=izMin:izMax; 
  Vxec=Vx(:,:,izStep); % reduce size of Vx for entire scan to Vxec of canal
  NStep=length(izStep);
  CmaxVoxel(3)=CmaxVoxel3+NStep-Nvxl; %z-axis range in voxels
  ecSlice=obliqueslice(Vxec,CmaxVoxel,[0,0,-1]);
  iSec=find(ecSlice==1,1); % points in ecSlice on wall of ear canal
  if 0 && isempty(iSec)
    izMax=izMax-1;
  else
    swBreak=1;
  end
end
CminVoxel(3)=CminVoxel(3)+NStep-Nvxl;

if iter==1
  [Xm,Ym,Zm]=meshgrid(vxl,vxl,1:NStep); % all 3 dimensions start at 1
  if ~swVerbose % default
    hf3=figure(ifig);
    ifig=ifig+1;
    hf3.Position=[200 figbottom 700 530];
    patch(isosurface(Xm,Ym,Zm,Vxec,0),'FaceColor','w');
    grid on;
    xlabel('x (voxels)');
    ylabel('y (voxels)');
    zlabel('z (voxels)');
    title({ss,['Iteration ',iis,', Surface in voxel units']});
    view(az,el);
    hold on;
  else
    ifig=Plot12(Xm,Ym,Zm,Vxec,ss,az,el,ifig);  % for debugging
  end
end
end


