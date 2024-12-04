% Shape_Part1, 10/7/24.  Main script to run Part 1 of the code to analyze
% ear-canal shape based on reading a STL file from a digital scan of the
% outer ear and ear canal.
% Adapted from batchdo12_JASA, 10/2/24.

close all, clear all
% External user needs a Shape folder that is appropriately defined. Present
% code is defined for use after Github download.
Shape = [fileparts(mfilename('fullpath'))] ; % parent folder name for ESS Shape
tmp = strfind(Shape, filesep);
Shape = Shape(1:tmp(end)) ;
% Shared code will have a browser to select STL file from folder STLs_JASA.
% Below are individual files of test ears used in paper.
%STLfile=...    % test ear A
%  'A021_Ax_004L_20210719_1509_training.training_L3DS-19-0043_5.2.0.801_Mo1.STL';
% STLfile=...     % test ear B
%   'A026_Ax_002R_20210818_1558_training.training_L3DS-19-0043_5.2.0.801_Mo1.STL';
% STLfile=...     % test ear C
%   'A024_Ax_003R_20210813_1520_training.training_L3DS-19-0043_5.2.0.801_Mo1.STL';
STLfile=...     % test ear D
   'A008_Ax_003R_20210521_1537_training.training_L3DS-19-0043_5.2.0.801_Mo1.STL';
answer = questdlg('Do you want to browse for STL file','FindSTL','Yes','No','Yes');
if strcmp(answer,'Yes')
    [STLfile,fileloc] = uigetfile(['..' filesep 'STLs_JASA' filesep '*.STL'],'Select an STL file');
end
SID=STLfile(1:4);
if ~strcmp(SID(1),'M') && ~strcmp(SID(1),'C')
  Ear=STLfile(12);
else
  Ear=STLfile(6);
end

% The input STL files analyzed by batchdo12_JASA are in folder 'STLs_JASA'.
% Shared code should allow external user to output MAT and PDF files output 
% from batchdo12, which is part 1 of the code that includes the manual step 
% of volume segmentation) to matFiles_External and PDF_External,
% respectively. 
% For testing purposes, this batch file now allows writing to
% matFiles_JASA and PDF_JASA, respectively, as these folders will be filled
% with results from the 58 valid test ears.
matFile=do12_JASA([Shape,'STLs_JASA' filesep],[Shape,'matFiles_JASA' filesep],...
  STLfile,SID,Ear);

% Add switch to shared code to control whether to output PDF&MAT file results.
%keyboard % For testing, uncomment to stop running before writing out MAT and PDF files

[~,fn]=fileparts(matFile);
pdfFile=fullfile([Shape,'PDF_JASA'],fn); % new PDF folder

hfAll=handle(sort(double(findobj(0,'type','figure'))));
N=length(hfAll);
for jj=1:N
  OutputPlot(4,hfAll(jj),[],pdfFile,jj); % 4 for multi-page PDF file
end
