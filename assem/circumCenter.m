function varargout = circumCenter(a, b, c)
%CIRCUMCENTER  Circumcenter of three points
%
%   CC = circumCenter(P1, P2, P3)
%
%   Example
%   circumcenter
%
%   Input: a, b, c are k x 2 matrix, the first dimension k is number 
%          of triangles, the second dimension is the coordinate of vertex.
%   Output: k x 2 matrix or [k vector, k vector]
%   See also
%     points2d, circumCircle, centroid
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-12-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.
% vectorized by Lei Zhang, 15 May 2012.
% the formula comes from Wikipedia

% pre-compute some terms
narginchk(3,3);
d=size(a,2);
if d~=2 && d~=3
    error('RPS;geops:circumCenter: input argument is invalid')
end
    
ah = sum(a .^ 2, 2);
bh = sum(b .^ 2, 2);
ch = sum(c .^ 2, 2);

dab = a - b;
dbc = b - c;
dca = c - a;

% common denominator
D = .5./(a(:,1).*dbc(:,2) + b(:,1).*dca(:,2) + c(:,1).*dab(:,2));
% D = D(:) % change D to columb vector

% center coordinates
xc =  (ah.*dbc(:,2) + bh.*dca(:,2) + ch.*dab(:,2)).*D;
yc = -(ah.*dbc(:,1) + bh.*dca(:,1) + ch.*dab(:,1)).*D;

if nargout <= 1
    varargout = {[xc yc]};
else
    varargout = {xc, yc};
end