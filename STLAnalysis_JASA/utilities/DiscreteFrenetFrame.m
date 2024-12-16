function [T,N,B,Kappa,Torsion]=DiscreteFrenetFrame(ds,dsM,Npts)
    % DiscreteFrenetFrame
    %   [T,N,B,Kappa,Torsion] = DiscreteFrenetFrame(ds,dsM,Npts);
    %   ds,  difference between adjacent points of discrete curve   
    %   dsM, norm of ds at each point
    %   Npts, # of points in discrete curve
    %   T, unit tangent vector along curve in discrete Frenet frame
    %   N, unit normal vector along curve in discrete Frenet frame
    %   B, unit binormal vector along curve in discrete Frenet frame
    %   Kappa = radian angle of curvature
    %   Torsion = radian angle of torsion
    if Npts<5  % need at least 5 data points
      T=[];
      N=[];
      B=[];
      Kappa=[];
      Torsion=[];
      return
    end
    Tlong=ds./dsM; % unit tangent vector to curve
    imI=1:(Npts-2);
    i0I=2:(Npts-1);
    T=Tlong(i0I,:);
    Tcross=cross(Tlong(imI,:),T,2);
    TcrossM=vecnorm(Tcross,2,2);
    inonzero=find(TcrossM>eps);
    izero=find(TcrossM<=eps);
    B=zeros(Npts-2,3);
    B(inonzero,:)=Tcross(inonzero,:)./TcrossM(inonzero);
    % Process rows with zero cross product.  Adjacent rows of ds with
    % collinear points are bad, producing zeros in TcrossM in denominator.
    if inonzero(1)>1
      B(1,:)=[-1,0,0];
      izero=izero(2:end);
    end % need to define initial B when initial Tcross is zero
    if ~isempty(izero)
      jj=1;
      while jj<=length(izero)
        kk=izero(jj);
        Tcross(kk,:)=Tcross(kk-1,:);
        TcrossM(kk,:)=TcrossM(kk-1);
        B(kk,:)=B(kk-1,:);
        jj=jj+1;
      end
    end
    N=cross(B,T,2);
    i0=2:(Npts-3);
    im=1:(Npts-4); % redefined
    ip=3:(Npts-2); % same as zrngMain in calling function
    % Constrain cosines to be in interval [-1,1], i.e., remove effects due
    % to finite precision arithmetic
    cosKappaAngle=max(-1,min(1,diag(T(im,:)*transpose(T(i0,:)))));
    cosTorsionAngle=max(-1,min(1,diag(B(i0,:)*transpose(B(ip,:)))));
    KappaAngle=acos(cosKappaAngle); % acos returns angle in [0,pi]
    TorsionAngle=acos(cosTorsionAngle);
    Kappa=KappaAngle./dsM(ip); % ip aligns Kappa with T etc.                                           
    Torsion=TorsionAngle./dsM(ip);
    T=T(i0,:);
    N=N(i0,:);
    B=B(i0,:); % all output arrays have length of Npts-4,
    % corresponding to interior values 3:(Npts-2) of input r
  end % DiscreteFrenetFrame
  