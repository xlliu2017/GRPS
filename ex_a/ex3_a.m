function A = ex3_a(x, y,seed)
% EX1_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% example with non-separable scales
% Usage : a = ex3_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015
narginchk(2,3);   
if nargin==2      % added by zhu to gurantee repeateable results
    seed='default';
end
rng(seed);
k=2;
m=2^(k+1)+1;
[k1,k2]=meshgrid(-2^k:2^k,-2^k:2^k);
c=((k1.^2+k2.^2)<=2^(2*k));
a=.6*(rand(m)-.5);
b=.6*(rand(m)-.5);
A=zeros(length(x),1); % modefied by Zhu, to gurantee A is a column vector
for i=1:length(x)
	V=a.*c.*cos(2*pi*(x(i)*k1+y(i)*k2))+b.*c.*sin(2*pi*(x(i)*k1+y(i)*k2));
	V=sum(sum(V));
	A(i)=exp(V);
end

end
