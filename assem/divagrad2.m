function D = divagrad2(pf,tf,A)
%% DIVAGRAD2 compute the discrete operator for div a grad u 
% on the points pf and the triangular tf, were A is the discrete
% coefficient of on the mass centers of the triangular
% Usage: [pf, tf, ~] = squaregeom(Nf)
%         A= cala(opt, pf,tf);
%         D =divagrad2(pf,tf,A)
%
% Shengxin Zhu, May 2, 2015
% Reference : 
if nargin==0
    D=test_divagrad2();
    return;
end
if size(A,2)~=1;
    error('RPS:divagrad2 currently only works for scalar conductivity ')
end
narginchk(3,3);               % check the number of argin
if size(pf,2)~=2
    if size(pf,1)==2
        pf=pf';
    else
         error('RPS:gdvagrad2 :1st argument is invalid')
    end
end
if size(tf,2)~=3
    if size(tf,1)==3
        tf=tf';
    else
         error('RPS:gdvagrad2 :2nd argument is invalid')
    end
end
np = size(pf, 1);             % number of points 
%-----------------A1A2A3 is a triangular----------------------------------- 
 if1 = tf(:,1);               % index of A1, index 1 of the fine triangular
 if2 = tf(:,2);               % index of A2
 if3 = tf(:,3);               % index of A3
p1=pf(if1,:);                 % coordinate of A1
p2=pf(if2,:);                 % coordinate of A2
p3=pf(if3,:);                 % coordinate of A3

cc = circumCenter(p1, p2, p3);  % what's this 
pm1 = (p2+p3)/2; pm2 = (p1+p3)/2; pm3 = (p1+p2)/2;
r1 = sqrt(sum((pm1-cc).^2,2)); % this is Gamma1
r2 = sqrt(sum((pm2-cc).^2,2)); % this is Gamma2
r3 = sqrt(sum((pm3-cc).^2,2)); % this is Gamma3

n12=p2-p1; n12=bsxfun(@rdivide,n12, sqrt(sum(n12.^2,2)));
n23=p3-p2; n23=bsxfun(@rdivide,n23, sqrt(sum(n23.^2,2)));
n31=p1-p3; n31=bsxfun(@rdivide,n31, sqrt(sum(n31.^2,2)));
[arf, g1x, g1y, g2x, g2y, g3x, g3y] = artrg(pf,tf);
D = sparse(np, np);
D=D+sparse(if1,if2, A.*dot([g2x,g2y],  n12,2).*r3,np,np);
D=D+sparse(if1,if3, A.*dot([g3x,g3y], -n31,2).*r2,np,np);
D=D+sparse(if2,if3, A.*dot([g3x,g3y],  n23,2).*r1,np,np);
D=D+D';
v1=n12.*[r3 r3]-n31.*[r2 r2];
D=D+sparse(if1,if1, A.*dot([g1x,g1y],v1,2),np,np);
v2=n23.*[r1 r1]-n12.*[r3 r3]; 
D=D+sparse(if2,if2, A.*dot([g2x,g2y],v2,2),np,np);
v3=n31.*[r2,r2]-n23.*[r1,r1]; 
D=D+sparse(if3,if3, A.*dot([g3x,g3y],v3,2),np,np);

function D=test_divagrad2()
Nf = 16 + 1; hf = 1/Nf;
[Xf,Yf]=meshgrid(linspace(0,1,Nf),linspace(0,1,Nf));
xf = Xf(:); yf = Yf(:); npf = length(xf);
x0 = xf + sin(pi/6)*yf; y0 = cos(pi/6)*yf;
DTf = DelaunayTri(x0,y0);
triplot(DTf)
tf = DTf.Triangulation;
pf = [xf,yf];
% remove flat triangles
arf = polyarea(xf(tf),yf(tf),2);
% arf = polyarea(xf(tf), yf(tf),2);
tf = tf(find(arf>hf^2/2*1e-2),:);
ntf = size(tf,1);

if1 = tf(:,1); if2 = tf(:,2); if3 = tf(:,3);

x1 = xf(if1); y1 = yf(if1);
x2 = xf(if2); y2 = yf(if2);
x3 = xf(if3); y3 = yf(if3);

xm = (x1+x2+x3)/3; ym = (y1+y2+y3)/3;

a = ones(size(xm));
%a = 1 + 0.5*sin(19*xm).*sin(17*ym);

D = divagrad2(pf, tf, a);

% u = (xf.^2+yf.^2)/4;
u = xf;  
y = D*u/hf^2;
trisurf(tf, xf, yf, y);
