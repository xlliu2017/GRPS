function [A,se] = ex5_a(x, y,seed)
% EX5_A returns the values of the coeficient a in the equation -\nabla a \nabla u=f
% in each triangulars. 
% channel example 
% Usage : a = ex5_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

% choice: 1, align to the grid, 2, not
narginchk(2,3)
if nargin==2
    seed='default';
end
choice=2;
a=zeros(size(x));
if choice==1
    xs=0.2*2-1;
    ys=0.4*2-1;
    xe=0.7*2-1;
    ye=0.8*2-1;
    for i=1:length(x)
        if (x(i)>=xs) && (x(i)<=xe) && (abs(y(i)-ys)<=.05)
            a(i)=1;
        elseif (y(i)>=ys-.05) && (y(i)<=ye) && (abs(x(i)-xe)<=.05)
            a(i)=1;
        end
    end
elseif choice==2
    xs=2*0.15-1; ys=2*0.4-1;
    xm=2*0.6-1; ym=2*0.2-1;
    xe=2*0.8-1; ye=2*0.65-1;
    
    width=1/2^6;
    for i=1:length(x)
        dlsm=det([1 x(i) y(i);1 xs ys;1 xm ym])/sqrt((xs-xm)^2+(ys-ym)^2);
        dlme=det([1 x(i) y(i);1 xm ym;1 xe ye])/sqrt((xe-xm)^2+(ye-ym)^2);
        if (dlsm>=0)&&(dlsm<width)&&((xs-x(i))*(xs-xm)+(ys-y(i))*(ys-ym)>0)...
                &&((xm-x(i))*(xm-xs)+(ym-y(i))*(ym-ys)>0)
            a(i)=1;
        end
        if (dlme>=0)&&(dlme<width)&&((xe-x(i))*(xe-xm)+(ye-y(i))*(ye-ym)>0)...
                &&((xm-x(i))*(xm-xe)+(ym-y(i))*(ym-ye)>0)
            a(i)=1;
        end
    end
end
se=rng(seed);
A=1+rand(size(x))+a.*(100+rand(size(x)));

end
