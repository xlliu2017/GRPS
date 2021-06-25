function A = ex9_a(x, y)
% EX9_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% checker board
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
for i = 1:length(x)
	j = floor((x(i)-xmin)/dx) + 1;
    k = floor((y(i)-ymin)/dy) + 1;
    A(i) = mod(j+k,2)*(a-1/a) + 1/a;
end

end
