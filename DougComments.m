% I’m inside the GitHub site in subfolder STLAnalysis_JASA of phal-shape.  I see three M files for Shape_Part1, Shape2_Part2 and InstallPathTmp in pal_shape and found other M functions 
% (we wrote) called by Part1 or Part2 in the utilities subfolder (and the 3rd party subfolders).
% 
% I’ll pretend I want to download everything and run Part1 first, as I cannot add anything to the Readme file if I don’t know how Part1 works. 
% So I’ve downloaded a ZIP file of phal-shape from GitHub and for some reason it is unzipped on my Mac and called phal-shap-main, 
% which is apparently the main branch of the repository.  So far so good.  I moved phase-shape-main to a permanent folder on my Mac, opened Matlab R2024b, 
% and navigated to the STLAnalysis_JASA folder as the home folder in the Matlab app.  The next step is to open do12_JASA function in the M editor as it’s called by Part1.
% 
% Method 1.  I edited the Matlab path to add this home folder and all its subfolders to its path. That was OK and I opened do12_JASA.m in the editor.
% Method 2.  I ran InstallPathTmp.m but that bombed.  The M file uses ‘\’ for the file separator for Windows, but this must be replaced by 
% filesep here and elsewhere in every M file so it will work on both Windows and Mac computers. 
% 
% Would you do that and make a new branch of the library for me to work with?   For example see lines 181-182 of do12_JASA.m. 
% Then I’d be happy to use InstallPathTmp if you think it is preferable to Method 1 (?).  
% I may find other platform differences later on but the file separator difference in the only one I can think of at the moment. 
% The Matlab functions ismac and ispc and other is* commands are available if that is helpful in the code.  
% Sorry about that but I think I need to be a pretend user on a Mac as that is the computer on which I have Matlab 
% (that’s currently approved for my use with Mathworks on the BTNRH license).  I’m also a little unclear on the scope of what’s needed in the Readme.mdb file, 
% but it seems like starting on its revision is the best way to proceed.   
% 
% I may refer users to the existing comments in one of our M files on how to use the Image Processing library function to work with the ear-canal STL file 
% rather than re-writing the explanation into the MDB file.  Do you think that would suffice? 
% The other question is my mind is the lack a license for subfolder gridtrimesh on Matlab File Exchange, although licenses are present there for polygon2voxel, extrema, 
% and altmany-export_fig-3.47.0.0.  I cannot log in to the Mathworks File Exchange without entering my BTNRH email address, 
% which I do not wish to do (as I cannot view these emails or at least assume I cannot).   
% Perhaps the MDB file can simply report that the code may be obtained from this File Exchange, what do you think?  
% I’m not sure it I asked you to ask Chris Stecker about that or not, as these external library issues do involve BTNRH.
% 
% You need not accomplish this before the holiday break as I have little time available the rest to this week, 
% but it would be helpful to have the revised M files early in December and I’ll continue testing on the M
% 
% Yes, you should make changes in M files (unless some new coding problem arises that needs input).  
% When ready, I’ll follow your git instructions on updating files on my Mac.   
% I was mistyping the name of the readme file coupled with Apple “helpfully" changing what I typed.  
% I’ll certainly refer to comments in M files by line numbers but perhaps you are saying that the readme file should contain all those instructions?  
% I will check out the files in utilities once I actually start to run coder just so we make it all as simple as possible to use for interested people.
% 
% Yeah, I guess we have to over-describe how properties of shape or geometry influence the sound field….  
% The Susan Voss ear-canal shape paper on which I was coauthor was rejected because it didn’t have any acoustics.  
% We had an acoustics/shape analysis in our manuscript and reviewer 1 still wants to reject it.