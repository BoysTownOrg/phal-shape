function [A,B,C,D,E,F,ResidualError]=fit_ellipse_Halir(x,y)
% Implementation from Radim Halir and Jan Flusser, "Nemerically stable
% direct least squares fitting of ellipses," Proc. 6th International Conf.
% in Central Europe on computer graphics and visualization (World Society 
% for Computer Graphics, Pizen, Czech Republic), p 125-132, 
%
% x,y are vectors of coordinates
Nx=length(x);
D1=[x.*x,x.*y,y.*y];                % quadratic part of the design matrix
D2=[x,y,ones(Nx,1)];                % linear part of the design matrix
S1=D1'*D1;                          % quadratic part of the scatter matrix
S2=D1'*D2;                          % combined part of the scatter matrix
S3=D2'*D2;                          % linear part of the scatter matrix
T=-inv(S3)*S2';                     % for getting a2 from a1
M=S1+S2*T;                          % reduced scatter matrix
M=[M(3,:)./2; -M(2,:); M(1,:)./2];  % premultiply by inv(C1)
[evec,eval]=eig(M);                 % solve eigensystem
cond=4*evec(1,:).*evec(3,:)-evec(2,:).^2; % evaluate a' C1 a
%a1=evec(:, find(cond>0));           % eigenvector for min. pos. eigenvalue
a1=evec(:, cond>0);                 % eigenvector for min. pos. eigenvalue
a2=T*a1;                            % eigenvector in lower part 
A=a1(1);                            % ellipse coefficients of:
B=a1(2);                            % A*x.^2+B*x.*y+C*y.^2+D*x+E*y+F=0
C=a1(3);
D=a2(1);
E=a2(2);
F=a2(3);
epsD=D1*a1+D2*a2; % algebraic error for each data point pair in conic
ResidualError=sqrt(sum(epsD.*epsD))/Nx;

