function [ux,uy]=grad(p,t,u)
%GRAD Compute the gradient of a PDE solution.
%
%       [UX,UY]=GRAD(P,T,U) returns grad(u) evaluated at
%       the center of each triangle.
%
%       The geometry of the PDE problem is given by the triangle data P,
%       and T.
%


%       A. Nordmark 12-22-94.
%       Copyright 1994-2003 The MathWorks, Inc.
%       $Revision: 1.8.4.2 $  $Date: 2010/10/08 17:14:47 $

np=size(p,1);
nt=size(t,1);
N=size(u,1)/np;

% Corner point indices
it1=t(:,1);
it2=t(:,2);
it3=t(:,3);

% Triangle geometries:
[ar,g1x,g1y,g2x,g2y,g3x,g3y]=artrg(p,t);

ux=u(it1,:).*g1x + u(it2,:).*g2x + u(it3,:).*g3x;
uy=u(it1,:).*g1y + u(it2,:).*g2y + u(it3,:).*g3y;

