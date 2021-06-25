function [A,se]= ex4_a(x, y,seed)
% EX4_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% random media 
% Usage : [A,se] = ex4_a(x, y,seed)
%         A=sex4_a(x,y);
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015
narginchk(2,3);         % modifed by Zhu
if nargin==2             % modifed by Zhu
    seed='default';   
end
n = length(x);
m=floor(log2(sqrt(n)))+1;
r=2;
se=rng(seed);           % add by Zhu, gurantee repeatable results
for k = 1:m
	a{k}=1/(1+r)+(1+r-1/(1+r))*rand(2^k);
end

% establish a cartesian grid
xmin=min(x);
xmax=max(x)+1e-8;
        
ymin=min(y);
ymax=max(y)+1e-8;
        
dx=(xmax-xmin)/2^m;
dy=(ymax-ymin)/2^m;
        
c=ones(2^m);
for i=1:m
	X=a{i};
	Y=ones(2^(m-i));
	c=c.*kron(X,Y);
end
A=zeros(size(x));     %modified by Zhu to    
for k=1:n
    i=floor((x(k)-xmin)/dx)+1;
    j=floor((y(k)-ymin)/dy)+1;
    A(k)=c(i,j);
end
A=A+1; % modified by Zhu to remove the case contrast is infinity
end
