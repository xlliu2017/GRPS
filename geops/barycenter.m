function [xb,yb]=barycenter(p,t)
%%
% p should be a nx2 array
% t should be a nt x 3 array 
%
xp=p(:,1);             % x-coordinate
yp=p(:,2);             % y-coordinate
xb=sum(xp(t),2)/3;  % barcenter coordinate
yb=sum(yp(t),2)/3;  % barcenter coordinate
%a=feval(funca,xb,yb);       % conductivity
%f=feval(funcf,xb,yb);       % function 
end