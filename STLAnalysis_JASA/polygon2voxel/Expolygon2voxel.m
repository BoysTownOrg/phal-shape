%function Volume=polygon2voxel(FV,VolumeSize,mode,Yxz)
% This function POLYGON2VOXEL will convert a Triangulated Mesh into a
% Voxel Volume which will contain the discretized mesh. Discretization of a 
% polygon is done by splitting/refining the face, until the longest edge
% is smaller than 0.5 voxels. Followed by setting the voxel beneath the vertice 
% coordinates of that small triangle to one.
%

%
%   % Load a triangulated mesh of a sphere
%   load sphere; 
%
%   % Show the mesh
%   figure, patch(FV,'FaceColor',[1 0 0]); axis square;
%
%   % Convert the mesh to a voxelvolume
%   Volume=polygon2voxel(FV,[50 50 50],'auto');
%
%   % Show x,y,z slices
%   figure,
%   subplot(1,3,1), imshow(squeeze(Volume(25,:,:)));
%   subplot(1,3,2), imshow(squeeze(Volume(:,25,:)));
%   subplot(1,3,3), imshow(squeeze(Volume(:,:,25)));
%
%   %  Show iso surface of result
%   figure, patch(isosurface(Volume,0.1), 'Facecolor', [1 0 0]);
%
close all,clear all
% Example2,
%
  % Make A Volume with a few blocks
  I = false(120,120,120);
  I(40:60,50:70,60:80)=1; I(60:90,45:75,60:90)=1;
  I(20:60,40:80,20:60)=1; I(60:110,35:85,10:60)=1;

  % Convert the volume to a triangulated mesh
  FV = isosurface(I);

  % Convert the triangulated mesh back to a surface in a volume
  J = polygon2voxel(FV,[120, 120, 120],'none'); 
  % Fill the volume
  J=imfill(J,'holes');

  % Differences between original and reconstructed
  VD = abs(J-I);

  % Show the original Mesh and Mesh of new volume
  figure, 
  subplot(1,3,1),  title('original')
    patch(FV,'facecolor',[1 0 0],'edgecolor','none'), camlight;view(3);
  subplot(1,3,2), title('converted');
    patch(isosurface(J,0.5),'facecolor',[0 0 1],'edgecolor','none'), camlight;view(3);
  subplot(1,3,3), title('difference');
    patch(isosurface(VD,0.8),'facecolor',[0 0 1],'edgecolor','none'), camlight;view(3); 
%
% Function is written by D.Kroon University of Twente (May 2009)
% last update (May 2019 at Demcon)
