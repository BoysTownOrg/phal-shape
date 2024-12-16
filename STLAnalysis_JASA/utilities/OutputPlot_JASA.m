function OutputPlot_JASA(swPrint,h,fnroot,swFirst)
%OutputPlot_JASA, 12/10/24. Use exportgraphics to output multi-page PDF
%file.
if swPrint==0
  return
end
set(h,'Color','white');
fnroot=[fnroot,'.pdf'];
if exist(fnroot,'file')~=2 % create if doesn't exist
  exportgraphics(h,fnroot,'Resolution',600);
else % PDF file exists append
  if exist('swFirst','var')==1 && swFirst==1
    exportgraphics(h,fnroot,'Resolution',600);
  else
    exportgraphics(h,fnroot,'Append',true,'Resolution',600);
  end
end

