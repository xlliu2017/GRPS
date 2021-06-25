function A = ex6_b(x, y)
% EX6_A returns the values of the coeficient a in the equation -\nabla a \nabla u=f
% in each triangulars. 
% percolation example 
% Usage : a = ex6_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

n = length(x);
m = floor(log2(sqrt(n)));

% establish a cartesian grid
xmin=min(x);
xmax=max(x)+1e-8;

ymin=min(y);
ymax=max(y)+1e-8;

dx=(xmax-xmin)/2^m;
dy=(ymax-ymin)/2^m;

kk = min(6, m);
R = rand(2^(kk));
c = 1 + 100*(R<0.5);
c = kron(c, ones(2^(m-kk)));

A = zeros(size(x));
for k=1:n
    i=floor((x(k)-xmin)/dx)+1;
    j=floor((y(k)-ymin)/dy)+1;
    A(k)=c(i,j);
end

end