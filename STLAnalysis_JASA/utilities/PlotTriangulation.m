function [ifig,hf1]=PlotTriangulation(ifig,Triangulation,az,el,ss,iis)
%PlotTriangulation, plot Triangulation data for canal

% Plot all data
hf1=figure(ifig);
ifig=ifig+1;
hf1.Position=[200 150 624 480];
trisurf(Triangulation,'EdgeColor','k');
view(az,el);
axis tight;
%title('From file');
colorbar;
xlabel('x  (mm)');
ylabel('y  (mm)');
zlabel('z  (mm)');
if ~isempty(ss)
  title({ss,['Iteration ',iis,', triangular mesh colored by z']});
else
  title('Triangular mesh colored by z');
end
end

