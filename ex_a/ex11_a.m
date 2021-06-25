function A = ex11_a(x, y)
% EX11_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% random checker board
% Usage : a = ex9_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

m = 64;
a = 10;
xmin = min(x);
xmax = max(x)+1e-8;
ymin = min(y);
ymax = max(y)+1e-8;
dx = (xmax - xmin)/m;
dy = (ymax - ymin)/m;
A=zeros(size(x)); % added by Zhu, preallocate for
rng('default')
r = binornd(1, 0.5, m^2, 1);
for i = 1:length(x)
	j = floor((x(i)-xmin)/dx) + 1;
    k = floor((y(i)-ymin)/dy) + 1;
    A(i) = r((j-1)*m+k)*(a-1/a) + 1/a;
end

end
