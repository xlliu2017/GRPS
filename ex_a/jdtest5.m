function A = jdtest5(x, y)
% EX5_A returns the values of the coefficient a in the equation -\nabla a \nabla u=f
% in each triangles. 
% interface example 
% Usage : a = ex5_a(x, y)
% 
% See also: TWOMESH_TEST
%  See Zhu Thesis F9

%A = 250*tanh(-3.*(.595576*(x+3.79761).^2 -y-10))+251;
A = 500*tanh(10*(y-x))+501;
end