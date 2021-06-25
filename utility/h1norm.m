function f=h1norm(p,t,u)
% H1NORM Calculate the H1 norm of piecewise linear function u 
% on triangulation (p,t)

%       Lei Zhang 11-18-2011

if nargin == 0
    test_h1norm();
    return
end

np = size(p,1);
nt = size(t,1);

[ux,uy]=grad(p,t,u);
ar = artrg(p,t);

f = sum((ux.*ux+uy.*uy).*ar);
f = f + l2norm(p,t,u)^2;
f = sqrt(f);

function test_h1norm()

h = 1/200;
[X,Y] = meshgrid(0:h:1,0:h:1);
x = X(:); y = Y(:);
dt = DelaunayTri(x,y);
p = dt.X; t = dt.Triangulation;

np = size(p,1);
it1 = t(:,1); it2 = t(:,2); it3 = t(:,3);

ar = artrg(p,t);
aod = ar/12; % Off diagonal element
ad = 2*aod; % Diagonal element
M=sparse(it1,it2,aod,np,np);
M=M+sparse(it2,it3,aod,np,np);
M=M+sparse(it3,it1,aod,np,np);
M=M+M.';
M=M+sparse(it1,it1,ad,np,np);
M=M+sparse(it2,it2,ad,np,np);
M=M+sparse(it3,it3,ad,np,np);

ar = artrg(p,t);

u1 = ones(np,1);
[ux1,uy1] = grad(p,t,u1);

u2 = x + 2*y;
[ux2,uy2] = grad(p,t,u2);

u3 = x.^5 + y.^3;
[ux3,uy3] = grad(p,t,u3);

u4 = sin(x).*sin(y);
[ux4,uy4] = grad(p,t,u4);

f1 = sum((ux1.*ux1+uy1.*uy1).*ar);
f1 = u1'*M*u1 + f1;
f1 = sqrt(f1);
% quad2d(@(x,y) 1,0,1,0,1)

f2 = sum((ux2.*ux2+uy2.*uy2).*ar);
f2 = u2'*M*u2 + f2;
f2 = sqrt(f2)
a2 = quad2d(@(x,y)(x+2*y).^2 + 5,0,1,0,1);
a2 = sqrt(a2)
abs(f2-a2)/abs(a2)

f3 = sum((ux3.*ux3+uy3.*uy3).*ar);
f3 = u3'*M*u3 + f3;
f3 = sqrt(f3)
a3 = quad2d(@(x,y)(x.^5+y.^3).^2 + 25*x.^8 + 9*y.^4,0,1,0,1);
a3 = sqrt(a3)
abs(f3-a3)/abs(a3)

f4 = sum((ux4.*ux4+uy4.*uy4).*ar);
f4 = u4'*M*u4 + f4;
f4 = sqrt(f4)
a4 = quad2d(@(x,y) 1 - (cos(x).*cos(y)).^2,0,1,0,1);
a4 = sqrt(a4)
abs(f4-a4)/abs(a4)