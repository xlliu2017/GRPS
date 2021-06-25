function [l2error, h1error, linferror] = error_2d(p,t,u,u0)
% ERROR_2D Calculate error in L2, H1 and LINF norm
%
%        calculate the numerical error on the fine mesh
%        fine mesh solution is taken as the exact solution
%        coarse mesh solution has been projected over the 
%        fine mesh

l2error = l2norm(p,t,u-u0)/l2norm(p,t,u0);
h1error = h1norm(p,t,u-u0)/h1norm(p,t,u0);
nz = find(u0~=0);
linferror = max(abs(u(nz)-u0(nz)))/max(abs(u0(nz)));

