function f = f_exponential(x, y)
% EX1_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% multiscale trignomeric example
% Usage : a = ex2_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

f=exp(-((x-0.5).^2+(y-0.5).^2)/0.1);
end