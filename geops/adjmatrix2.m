function [adj, T2T]= adjmatrix2(dt,type,layer)
%%[adj,T2T]=]AJDMATRIX2(dt,type,layer) returns the adjacent 
% matrix for a vertices/edges to the simplices in the triangularization dt.
% INPUT:
%       dt: is the TRIREP class, 
%     type: can be 'E' or 'P'. 'E' results in the vertices to simplce 
%           adjacent matrix, and 'E' results in the edges to simplices 
%           adjacent matrix. 
%     layer: is used to indicate how many layers that containt each vertex and edge. 
%
% OUTPUT: adj: the adjacent matrix 
%         T2T: the adjacent matrix for simplex to simplex 
% 
% Usage: 
%         [p, t,dt] = squaregeom(40);
%         tic;[adj,t2t] = adjmatrix2(dt,'p',4);toc;
%         tic;[adj,t2t] = adjmatrix0(dt,'e',4);toc

%  See Also: ADMATRIX0
%
%
%  by Shengxin Zhu, LCP,IAPCM, MAY/2015/

% Reference:
if nargin == 0
    test_adjmatrix2();
    return
end
narginchk(1,3);
if ~isa(dt,'TriRep')
    error('RPS:geops:adjmatrix2:input argument is invalid')
end
inlayer=0;                  % by default increasing layers==0
if nargin>1
    if ~strcmpi(type,'P') &(~strcmpi(type,'E'))
        error('EPS:geops:adjmatrix2;the 2nd argument is invalid')
    end
end
if nargin==3
    if layer>0 & layer <10
        inlayer=layer-1;
    else
        error('RPS:geops:adjmatrix2: in increasing layers is invalid.')
    end
end
%%----------------Begin main operations-------------------------
switch type
    case 'P'
        adj=PointsInSimplex(dt);   % comput the point to simplex adjacent matrix
        if inlayer | nargout==2
            T2T=Simplex2Simplex(adj);
            if inlayer
            adj=IncreasingLayers2(adj,T2T,inlayer);
            end
        end
    case 'E'
        adj=EdgesInSimplex(dt);    % compute the point to simplex adjacent matrix
        if inlayer | nargout==2
            adj0=PointsInSimplex(dt);
            T2T=Simplex2Simplex(adj0);
            if inlayer
                adj=IncreasingLayers2(adj,T2T,inlayer);
            end
        end
end
%% ------------test functions
function adj=test_admatrix2()
tic;[p, t,dt] = squaregeom(60);t1=toc;
tic;[adj,T2T] = adjmatrix2(dt,'P',5);t2=toc;
fprintf('Test admatrix OK.\n')
fprintf('Test admatrix2 takes %9.7f\n',t2)
tic;[adj,T2T] = adjmatrix2(dt,'E',5);t2=toc;
