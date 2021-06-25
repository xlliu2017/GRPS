function [p, t, dt] = squaregeom(n)
%% SQUAREGEOM create a regular triagulation for a square [0, 1] x [0, 1]
% USEAGE: [p, t, dt] = squaregeom(n)
% See also DELAUNAYTRI, POLYAREA,
% Lei Zhang Oct 04, 2012, modified and commented by Shengxin Zhu,18/11/2014

if nargin == 0
    test_squaregeom();
    return
end

h = 1/n; 

[X, Y] = meshgrid(linspace(0,1,n+1), linspace(0,1,n+1));
x = X(:); y = Y(:);
x0 = x + sin(pi/6)*y; y0 = cos(pi/6)*y;  %% why here

%dt = DelaunayTri(x0,y0);  %% why not x, y?  they produce the sameresult
dt=DelaunayTri(x,y);
p = [x ,y];
t = dt.Triangulation;

% remove flat triangles
ar = polyarea(x(t),y(t),2);   %% area

t = t(find(ar>h^2/2*1e-2),:);

end

function test_squaregeom()

[p, t] = squaregeom(10);
disp(['n=10'])
disp(['number of triangles are ', num2str(size(t,1))])
triplot(t, p(:,1), p(:,2));
end