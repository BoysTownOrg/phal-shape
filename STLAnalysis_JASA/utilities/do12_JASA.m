function fnM = do12_JASA(STLdir,dirMATdata,STLfile,SID,Ear)
% do12_JASA, 9/30/24, is based on do12, 8/20/21.
% This app is used to mask out ear scan parts except the canal stalk.

mypath=fileparts(mfilename('fullpath'));
currentpath=cd(mypath); %saves old path in currentpath, changes to mypath
fn = [STLdir STLfile];
STLfileOverride={'C004_L.STL'}; 
% Files listed have ear canal as 2nd largest z value not largest as in adults 
iOver=find(contains(STLfileOverride,STLfile)); % usually [] else used in VoxelCalc7

delZm=0.25; % sets mesh size in units (mm) of STL file, 0<delX<=1.0 mm.
if ~strcmp('C',SID(1)) % adult or anything except a child
  zmmCut=14;
else
  zmmCut=10;
end
% If STL has only 1 peak near TM, then the z axis distance in
% the voxelization is zmmCut in mm, i.e., round(zmmCut/delZm) in voxels.
% This only applies to iter==1 to rotate the canal to approximate verticaL.
if delZm>1 || delZm<=0
  disp('Error, Define delZm (mm) so that 0<delZm<=1');
  return
end
ifig=1;
if strcmpi('L',Ear)  
    az=110;
else
    az = 110-180;
end
STLsplit = strsplit(STLfile,'_');
if strcmp('M',SID(1))
  runnum = 0;
  STLdate = datestr(now,29);
  ss=[SID,', ',Ear,', Date ',STLdate]; % use in figures
else
    if length(STLsplit) < 3
        runnum = 0;
        STLdate = datestr(now,29);
    else
        runnum = STLsplit{3}(1:end-1);
        STLdate = STLsplit{4};
    end
    ss=[SID,', ',Ear, ', Run ',runnum,', Date ',STLdate]; % use in figures
end
el=10;
Niter=2;
Triangulation0=stlread(fn); % load initial data
F=Triangulation0.ConnectivityList; % Initialize for while loop
Triangulation=Triangulation0;
zLen=zeros(1,Niter);
sLen=zeros(1,Niter);
sMinMax=zeros(3,Niter);
iter=1;
iis=int2str(iter);
disp(['Iteration ',iis]);
[ifig,Vx,izMin,izMax,CmaxVoxel,CminVoxel]=...
  VoxelCalc7(iis,ifig,delZm,Triangulation,F,ss,az,el,zmmCut,iOver,iter);
co=colororder;
figure(ifig);
% el1Max,ellipse model struct for locations at max z;
% el1Min,ellipse model struct for locations at min z;
[el1Max,el1Min]=DoEndSlices(iter,Vx(:,:,izMin:izMax),CmaxVoxel,...
  CminVoxel,ifig,co,ss,iis);
ifig=ifig+1;
% sMinMax is current canal centerline for rotating canal to new z axis
% x and y are interchanged to restore the order associated with use of
% obliqueslice in DoEndSlices, and -x is used because this axis was
% inverted after obliqueslice (where it was the y axis!)
sMinMax(:,iter)=[el1Max.Ymean-el1Min.Ymean,...
  -el1Max.Xmean+el1Min.Xmean,el1Max.Z0-el1Min.Z0]; % in voxels
% Rotate canal stalk to vertical using phi and theta
costheta=sMinMax(3,iter)/norm(sMinMax(:,iter)); % angle between z axis and vMinMax
sintheta=sqrt(1-costheta^2);
thetaDeg=acos(costheta)*180/pi;
phi=atan2(sMinMax(2,iter),sMinMax(1,iter));
cosphi=cos(phi);
sinphi=sin(phi);
phiDeg=phi*180/pi;
disp(['Theta: ',num2str(thetaDeg,'%5.1f'),' deg']);
disp(['Phi: ',num2str(phiDeg,'%5.1f'),' deg']);

% Matrix rotation by azimuthal phi and polar theta angles
Bmat=[cosphi,sinphi,0;... % Product of 2 simple rotation matrices
  -costheta*sinphi,costheta*cosphi,sintheta;...
  sintheta*sinphi,-sintheta*cosphi,costheta];
Vrot=transpose(Bmat*transpose(Triangulation.Points));
%size(Triangulation.Points)
Triangulation=triangulation(F,Vrot); % Update for next iteration
% Rotate centers of min/max ellipses to find new z distance
% for next iteration in voxels.
iter=2;
iis=int2str(iter);
disp(['Iteration ',iis]);
[ifig,Vx,izMin,izMax,CmaxVoxel2,CminVoxel2,TriangulationV]=...
  VoxelCalc7(iis,ifig,delZm,Triangulation,F,ss,az,el,zmmCut,[],iter);

figure(ifig);
ifig=ifig+1;
DoEndSlices(iter,Vx(:,:,izMin:izMax),CmaxVoxel2,...
  CminVoxel2,ifig,co,ss,iis);

swStep1=0; % swStep1==1 if read file, plot, rotate canal to vertical, and 
% then return.  Else proceed to further steps to save a MAT file of
% processed ear-canal image from volumeSegmenter.   
if swStep1
  return
end
Vxec=Vx(:,:,1:izMax); % scan data for Vol Seg
% Proceed for use of volumeSegmenter app, see MATLAB documentation
% It's OK to try it out after the ReadEarSTL function works.
% Actions in volumeSegmenter:
% 0. Wait (a long time) until the app is fully open and an OK warning panel
% is displayed on top of the app GUI.Do not click the OK panel until 
% you've closed the volumeSegmenter App.  Left lower panel of app GUI 
% 3-D Display shows the complete ear scan.  You can click in this panel, 
% click it and rotate it.  When done, single-click main panel called Slice.
% This Slice panel is on slice 1, the most lateral or bottom slide normal
% to the z axis. The panel reports the number of slices, with the top slice
% the most medial that may not show any of the canal (i.e., all black).
% Steps below are for left ear.
% 1. In app GUI, click on DRAW tab.
% 2. Select Polygon tool 
% 3. At lower right of Slice panel, click the right arrow to advance until 
% you see a slice that is entirely distinct from any part of the concha.  
% 'Distinct' means that the black interior of the canal has a complete 
% wall, which is represented by white pixels (value 1), around which is one 
% or more black pixels.
% 4a. In that slice, click once in the lower right black to the right of any
% canal pixel, and one in the lower left black to the left of any canal
% pixel. These two clicks are just above the bounding rectance of the
% slice.
% 4b. Click once in the black to the upper left of the canal.  Then click 
% once in the black to the upper right of any white canal pixel.  This line 
% segment is entirely in black. If desired, click one more time in the 
% black in the upper part more or less above the initial click.  Or, hit
% Enter.  The Enter tells the app that the polygon is complete and it draws
% the final line segment back to the first click.
% 5. At lower right of Slice panel, click the arrow or in the horizontal
% bar to advance to the final slice.
% 6. At lower left of Slice panel, click the left arrow in the horizonal 
% bar to move to a preceding slice. Keep clicking until you encounter the
% first canal portion (1 or more white pixels) such that there is some
% inner area of black pixels. 
% 7. Repeat the actions of step 4 to draw a polygon around this canal slice. 
% 8. Click Auto-interpolate button in the upper row of the DRAW GUI. You
% will see a horizontal band of blue that highlights all of the slices of 
% the canal. Now you can check that the canal stalk is in the masked blue 
% area for all the intermediate slices.  Else, exit and redo.  See MATLAB
% instructions to use manual interpolate button, but probably no advantage.
% 9. In app GUI, click on SEGMENTER tab.
% 10. Click on Save Labels button.  If drop-down menu opens, select save as
% MAT file.  A Variable Type window pops up. Click OK to save segmentation
% as logical, i.e., as a binary image of 1's and 0's.  A Windows browser
% panel opens up.  Click Save to save it with the default names labels.mat in
% the default folder, which is that of the M file that is running.
% 11.  Close the VolumeSegmenter app by clicking X in the upper right of
% the window.
% 12.  Click OK in the small MATLAB window, which continues execution of
% the running M function.
% Explanation.  Fig. 4 shows the most medial and most lateral slices of the
% original canal scan (slice normal to z axis). 
% Fig. 7 shows the most medial and most lateral slices of the
% rotated canal scan (slice normal to z axis) that you saved in the MAT
% file.  This ends the processing in this M file.  The labels.mat is
% over-written each time you run this M file, but results for subsequent
% scan analyses are written into the fn file, which is in dirMATdata.
%    For right ear, the canal stalk in steps 3-4 tends to be in upper part 
%    of slice rather than lower.   
Vxec8=uint8(Vxec);
volumeSegmenter(Vxec8); % open new app from Image Processing Toolkit
% at max, choose largest slice number at which there is a non-zero space
% contained within the pixels.
CreateStruct.Interpreter='tex';
CreateStruct.WindowStyle='non-modal'; % non-modal so msgbox moves to background
waitfor(msgbox(...
  '\fontsize{12} Click OK after closing Volume Segmenter',...
  CreateStruct));

cd(currentpath); %change back to initial path
s=load(['.' filesep 'utilities' filesep 'labels']); % temp file created by volumeSegmenter
delete(['.' filesep 'utilities' filesep 'labels.mat']); % % delete temp file, no need to keep these
labels1=s.labels;
mVxec=Vxec8&labels1; % return masked canal stalk after dialog box closed
imV=find(mVxec); % points in ecSlice on wall of ear canal
if isempty(imV)
  disp('imV should not be empty');
  return
end
[imx,imy,imz]=ind2sub(size(mVxec),imV); % defines canal stalk
[imx1,imx2]=bounds(imx);
[imy1,imy2]=bounds(imy);
[imz1,imz2]=bounds(imz);
if imx2==imx1 || imy2==imy1 || imz2==imz1
  disp('lower and upper bounds should differ');
  return
end
mVxec=mVxec(imx1:imx2,imy1:imy2,imz1:imz2); % minimal cube for canal stalk
Nimx=imx2-imx1+1;
Nimy=imy2-imy1+1;
NStep=imz2-imz1+1;
[Xm,Ym,Zm]=ndgrid(1:Nimx,1:Nimy,1:NStep); % each dimension starts at 1
[imVx,imVy]=find(mVxec(:,:,NStep));
CmaxVoxel=[mean(imVx),mean(imVy),NStep]; % mean center of top z slice
ifig=Plot12(Xm,Ym,Zm,mVxec,ss,az,el,delZm,ifig);
[imVx,imVy]=find(mVxec(:,:,1));
CminVoxel=[mean(imVx),mean(imVy),1]; % mean center of top z slice
figure(ifig);
DoEndSlices(iter,mVxec,CmaxVoxel,CminVoxel,ifig,co,ss,iis);

if delZm==1
  sDel='D1';
else
  sDel=['Dp',int2str(round(delZm*100))];
end
t=datestr(now,30);
STLsplit = strsplit(STLfile,'_');
if strcmp('M',SID(1)) 
    runnum = '000';
    STLdate = datestr(now,29);
else
    if length(STLsplit) < 3
        runnum = '000';
        STLdate = datestr(now,29);
    else

        runnum = STLsplit{3}(1:end-1);
        STLdate = STLsplit{4};
    end
end
fnM=fullfile(dirMATdata,[SID,'_',Ear,'_',runnum,'_',sDel,'_',STLdate,'_',t]);
save(fnM,'mVxec','delZm','NStep','Nimx','Nimy','Xm','Ym','Zm','az','el',...
  'STLfile','TriangulationV','ss','iis','thetaDeg','phiDeg');
% The fnM MAT file holds the canal scan results used in the 2nd and final
% iteration to parameterize the canal image in terms of area distance
% function and related variables. fnM is time-stamped and includes
% SID, Ear, and voxel edge length delZm. 
% 8/20/21. Include thetaDeg and phiDeg in saved results as part of
% version=2, but need not rerun old data from version==2.
