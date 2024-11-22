function ifig=Plot12(Xm,Ym,Zm,mVxec,sfn,az,el,delZm,ifig)
    hf12=figure(ifig);
    ifig=ifig+1;
    hf12.Position=[200 150 624 480];
    hpa=patch(isosurface(Xm,Ym,Zm,mVxec,0),'FaceColor','w');
    hpa.FaceAlpha=0.8;
    grid on;
    xlabel('x (voxels)');
    ylabel('y (voxels)');
    zlabel('z (voxels)');
    if ~isempty(sfn)
      title({sfn,'Canal surface in voxels'});
    else
      title('Canal surface in voxels');
    end
    view(az,el);
    set(gca,'DataAspectRatio',[1,1,1]);
    hold on;
    
    hf2=figure(ifig);
    ifig=ifig+1;
    set(hf2,'Visible','off');
    montage(~mVxec,'BackgroundColor',[1,1,1]);%    montage(mVxec);
    hf2.Visible='on';
    if ~isempty(sfn)
      title({sfn,['Iteration 1',...
        ', montage of 2D canal slices from bottom (medial end) to top (lateral end)'],...
        ['Slices vertically spaced along z-axis by ',num2str(delZm,2),...
        ' mm']});
    else
      title({...
        'Montage of 2D canal slices from bottom (medial end) to top (lateral end)',...
        ['Slices vertically spaced along z-axis by ',num2str(delZm,2),...
        ' mm']});
    end
    hf2.Position=[200,200,720,720];
    drawnow;
  end