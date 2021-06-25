function A = ex7_a(x, y,seed)
% EX7_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% channel example 
% Usage : a = ex7_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

if nargin==2
    seed='default';             % add by Zhu to gurantee repeatable result
end
n = length(x);
m = floor(log2(sqrt(n)))+1;
r = 2;
se=rng(seed);                  % add by Zhu to gurantee repeatable result
for k=1:m
    a{k}=1/(1+r)+(1+r-1/(1+r))*rand(2^k);
end

% establish a cartesian grid
xmin=min(x(:));
xmax=max(x(:))+1e-8;

ymin=min(y(:));
ymax=max(y(:))+1e-8;

dx=(xmax-xmin)/2^m;
dy=(ymax-ymin)/2^m;

c=ones(2^m);
for i=1:m
    X=a{i};
    Y=ones(2^(m-i));
    c=c.*kron(X,Y);
end
A1=zeros(n,1);
for k=1:n
    i=floor((x(k)-xmin)/dx)+1;   % here is a bug 
    j=floor((y(k)-ymin)/dy)+1;   
    A1(k)=c(i,j);
end

a1=zeros(size(x));
a2=zeros(size(x));

xs=0.2*2-1;
ys=0.4*2-1;
xe=0.7*2-1;
ye=0.8*2-1;
for i=1:length(x)
    if (x(i)>=xs) && (x(i)<=xe) && (abs(y(i)-ys)<=.05)
        a1(i)=1;
    elseif (y(i)>=ys-.05) && (y(i)<=ye) && (abs(x(i)-xe)<=.05)
        a1(i)=1;
    end
end

xs=2*0.15-1; ys=2*0.4-1;
xm=2*0.6-1; ym=2*0.2-1;
xe=2*0.8-1; ye=2*0.65-1;

width=0.05;
for i=1:length(x)
    dlsm=det([1 x(i) y(i);1 xs ys;1 xm ym])/sqrt((xs-xm)^2+(ys-ym)^2);
    dlme=det([1 x(i) y(i);1 xm ym;1 xe ye])/sqrt((xe-xm)^2+(ye-ym)^2);
    if (dlsm>=0)&&(dlsm<width)&&((xs-x(i))*(xs-xm)+(ys-y(i))*(ys-ym)>0)...
            &&((xm-x(i))*(xm-xs)+(ym-y(i))*(ym-ys)>0)
        a2(i)=1;
    end
    if (dlme>=0)&&(dlme<width)&&((xe-x(i))*(xe-xm)+(ye-y(i))*(ye-ym)>0)...
            &&((xm-x(i))*(xm-xe)+(ym-y(i))*(ym-ye)>0)
        a2(i)=1;
    end
end

A1= 1 + rand(size(x)) + max(A1)/4*a1.*(1+rand(size(x)));   %
A2 = 1 + rand(size(x)) + max(A1)/4*a2.*(1+rand(size(x)));  % is that right?

A = A1 + A2;

end
