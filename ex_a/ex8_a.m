function A = ex8_a(x, y)
% EX5_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% interface example 
% Usage : a = ex5_a(x, y)
% 
% See also: TWOMESH_TEST
% Lei Zhang, May 02, 2015

A = 1+ 1000*(x<y);

end
