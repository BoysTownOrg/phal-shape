function OutputPlot(swPrint,h,fn,fnroot,swFirst)
%OutputPlot, 5/3/22. For swPrint==4, use exportgraphics for export_fig.
%OutputPlot, 2/28/19.  Generate output file based on swPrint value
% swPrint: 1 for EPS, 2 for EPS and PNG, 3 for EPS, PNG, multi-page PDF, 
%          4 for multi-page PDF, 5 for single-page PDF
if swPrint==0
  return
end
set(h,'Color','white');
if swPrint==5
  export_fig(h,fn,'-pdf','-r300');
elseif swPrint<4
  export_fig(h,fn,'-eps','-r600');
  if swPrint>1
    export_fig(h,fn,'-png','-r300');
    if swPrint==3
      if exist([fnroot,'.pdf'],'file')==1 % create
        export_fig(h,fnroot,'-pdf','-r300');
      else % append
        export_fig(h,fnroot,'-pdf','-r300','-append');
      end
    end
  end
elseif swPrint==4
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
elseif swPrint==6
  print(h,'-dpdf','-bestfit',fn);
elseif swPrint==7
  export_fig(h,fn,'-png','-r450'); % -r600 generates warnings
elseif swPrint==8
  export_fig(h,fn,'-jpg','-r450'); % -r600 generates warnings
end
if swPrint~=4
  savefig(h,fn);
end

