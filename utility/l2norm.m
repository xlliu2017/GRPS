function f = l2norm(p,t,u)
% L2NORM Calculate the l2 norm of u over triangulation (p,t)

%       Lei Zhang 11-18-2011

if nargin == 0
    test_l2norm();
    return
end

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

f = u'*M*u;
f = sqrt(f);

function test_l2norm()

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

u1 = ones(np,1);
u2 = x + 2*y;
u3 = x.^5 + y.^3;
u4 = sin(x).*sin(y);

f1 = u1'*M*u1; f1 = sqrt(f1)
% quad2d(@(x,y) 1,0,1,0,1)
f2 = u2'*M*u2; f2 = sqrt(f2)
a2 = quad2d(@(x,y)(x+2*y).^2,0,1,0,1); a2 = sqrt(a2)
abs(f2-a2)/abs(a2)
f3 = u3'*M*u3; f3 = sqrt(f3)
a3 = quad2d(@(x,y)(x.^5+y.^3).^2,0,1,0,1); a3 = sqrt(a3)
abs(f3-a3)/abs(a3)
f4 = u4'*M*u4; f4 = sqrt(f4)
a4 = quad2d(@(x,y)(sin(x).*sin(y)).^2,0,1,0,1); a4 = sqrt(a4)
abs(f4-a4)/abs(a4)