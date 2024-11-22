function [ellip,swxy,p] = do_fit_ellipse_Halir(x,y,ifig1,ecolor,mk,mk2)
% wrapper function for function fit_ellipse_Halir. swxy not used.
[A0,B0,C0,D0,E0,F0,ResidualError]=fit_ellipse_Halir(x,y);
if any(imag([A0,B0,C0,D0,E0,F0]))
  ellip=[];
  swxy=-1;
  p=[];
  return
end
temp=sqrt((A0-C0)^2+B0^2);
denom=4*A0*C0-B0^2;
a=sqrt(2*(A0*E0^2+C0*D0^2-B0*D0*E0-denom*F0)*(A0+C0+temp))/denom;
b=sqrt(2*(A0*E0^2+C0*D0^2-B0*D0*E0-denom*F0)*(A0+C0-temp))/denom;
X0=(-2*C0*D0+B0*E0)/denom;
Y0=(-2*A0*E0+B0*D0)/denom;
if B0~=0
  angle=atan2(C0-A0-temp,B0); % angle in [-pi,pi]
  if angle>=pi/2 % constrain -pi/2<angle<=pi/2
    angle=angle-pi;
  elseif angle<-pi/2
    angle=angle+pi;
  end
else
  if A0<C0
    angle=0;
  elseif A0>C0
    angle=pi;
  end % what if B0==0 & A0==C0? circle?
end
if a<b
  temp2=a;
  a=b;
  b=temp2;
  if angle>0  % angle contrained to [-pi/2,pi/2]
    angle=angle-pi/2;
  else
    angle=angle+pi/2;
  end
  swxy=1;
else
  swxy=0;
end
ca=cos(angle);
sa=sin(angle);
ellip.a=a;
ellip.b=b;
ellip.area=pi*a*b; % area of ellipse
ellip.EqDiameter=2*sqrt(a*b); % equiv diam of circle with same area as ellipse
ellip.eccentricity=sqrt(1-(b/a)^2); % eccentricity of ellipse
ellip.cos=ca;
ellip.sin=sa;
ellip.angleDeg=180*angle/pi;
ellip.X0=X0;
ellip.Y0=Y0;
centerFit=[X0;Y0];
ellip.Xmean=mean(x); % mean data
ellip.Ymean=mean(y);
ellip.errCtr=norm([ellip.X0-ellip.Xmean,ellip.Y0-ellip.Ymean]);
ellip.ResidualError=ResidualError;

if ifig1>0 && ~isempty(a) % else not a fitted ellipse
  R=[ca, -sa; sa, ca];
  theta_r         = linspace(0,2*pi);
  ellipse_x_r     = a*cos(theta_r );
  ellipse_y_r     = b*sin(theta_r );
  rotated_ellipse=centerFit+R*[ellipse_x_r;ellipse_y_r];
  plot(rotated_ellipse(1,:),rotated_ellipse(2,:),'Color',ecolor);
  set(gca,'DataAspectRatio',[1,1,1]);
  hold on;
  ver_line        = [ [0,0]; b*[-1 1] ];
  horz_line       = [ a*[-1 1]; [0,0] ];
  new_ver_line    = centerFit+R*ver_line;
  new_horz_line   = centerFit+R*horz_line;
  p=plot(new_ver_line(1,:),new_ver_line(2,:),'Color',ecolor);
  plot(new_horz_line(1,:),new_horz_line(2,:),'Color',ecolor);
  hp0=plot(X0,Y0,'Marker',mk,'MarkerFaceColor',ecolor);
%  hp1=plot(round(X0),round(Y0),'Marker',mk,'Color','k');
  hp2=plot(ellip.Xmean,ellip.Ymean,'Marker','s','Color','k');
  plot([ellip.Xmean,X0],[ellip.Ymean,Y0],'k:');
  ecolor2=ecolor+0.1; % lighten color
  msize=3;
  plot(x,y,'LineStyle','none','Color',ecolor2,'Marker',mk2,'MarkerSize',msize);
  set(gca,'YDir','reverse'); % to align with obliqueslice figure
%  xlabel('Col-Y');
%  ylabel('Row-X reverse');
  xlabel('Abscissa');
  ylabel('Ordinate');
  %hl=legend([hp0,hp2],{'Model mean','Data mean'},'Location','Northoutside');
 % set(hl,'Orientation','horizontal','AutoUpdate','off');
%   lpos=get(hl,'Position');
%   lpos=lpos+[0,0.10,0,0];
%   set(hl,'Position',lpos);
  ellip.btxt=centerFit+R*ver_line/2; % for AAS
  ellip.atxt=centerFit+R*horz_line/2;
%  ellip.hl=hl;
  ellip.hl=[];
else
  ellip.btxt=[];
  ellip.atxt=[];
  ellip.hl=[];
end
