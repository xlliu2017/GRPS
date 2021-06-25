function A = ex2_a(x, y)
% EX1_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% multiscale trignomeric example
% Usage : a = ex2_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

e1=1/5; e2=1/13; e3=1/17; e4=1/31; e5=1/65;
A=1/6*((1.1+sin(2*pi*x/e1))./(1.1+sin(2*pi*y/e1))...
    +(1.1+sin(2*pi*y/e2))./(1.1+cos(2*pi*x/e2))...
    +(1.1+cos(2*pi*x/e3))./(1.1+sin(2*pi*y/e3))...
    +(1.1+sin(2*pi*y/e4))./(1.1+cos(2*pi*x/e4))...
    +(1.1+cos(2*pi*x/e5))./(1.1+sin(2*pi*y/e5))+sin(4*x.^2.*y.^2)+1);

end

