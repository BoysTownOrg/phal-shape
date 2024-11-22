function [el1Max,el1Min]=DoEndSlices(iter,Vxec,CmaxVoxel,CminVoxel,...
  ifig,co,ss,iis)
% DoEndSlices. The goal of this function is to plot the most medial and
% most lateral slices of the ear canal image. This is used
% to confirm that these canal areas are approximately lined up.  However,
% obliqueslice produces a slice no longer lined up along the original z
% axis because of the axes rotations in the calling function.  A single
% call to DoEndSlices suffices to check that the canal is approximately
% vertically aligned along the rotated z axis.
ecSliceMax=obliqueslice(Vxec,CmaxVoxel,[0,0,-1]);
ecSliceMin=obliqueslice(Vxec,CminVoxel,[0,0,-1]);
iSecMax=find(ecSliceMax==1); % points in ecSlice on wall of ear canal
[iSecRowMax,iSecColMax]=ind2sub(size(ecSliceMax),iSecMax);

[el1Max,~,pMax]=do_fit_ellipse_Halir(iSecColMax,iSecRowMax,ifig,...
  co(1,:),'o','x');
iSecMin=find(ecSliceMin==1); % points in ecSlice on wall of ear canal
[iSecRowMin,iSecColMin]=ind2sub(size(ecSliceMin),iSecMin);
[el1Min,~,pMin]=do_fit_ellipse_Halir(iSecColMin,iSecRowMin,ifig,...
  co(2,:),'d','+');
if iter==1
  el1Max.Z0=CmaxVoxel(3);
  el1Min.Z0=CminVoxel(3);%el1Min.Z0=1;
end
grid on;
title({ss,['Iteration ',iis,...
  ', z-slices at ends of canal stalk']});
hl=legend([pMax,pMin],{'Medial','Lateral'});
hl.Location='NorthWest';