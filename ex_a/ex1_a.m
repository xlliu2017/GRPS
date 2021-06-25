function A = ex1_a(x, y,seed)
% EX1_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% random media
% Usage : a = ex1_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015
narginchk(2,3);
if nargin==2
    seed='default';
end
rng(seed);
m = 256;
g = 0.3 + 2*rand(m,m);
xmin = min(x);
xmax = max(x)+1e-8;
ymin = min(y);
ymax = max(y)+1e-8;
dx = (xmax - xmin)/m;
dy = (ymax - ymin)/m;

A = zeros(size(x));

for i = 1:length(x)
	j = floor((x(i)-xmin)/dx) + 1;
    k = floor((y(i)-ymin)/dy) + 1;
    A(i) = g(j,k);
end

end
