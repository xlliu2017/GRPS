function f = f_quasiDelta(x, y)
% EX1_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% multiscale trignomeric example
% Usage : a = ex2_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015
f= zeros(size(x));
f((abs(x-0.5)<0.02)&(abs(y-0.5)<0.02))=100;
end

