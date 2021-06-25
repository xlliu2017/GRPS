function P2P=Points2Points(T)
%% P2P=Points2Points(T) Return a sparse binary matrix of np*np,
% where np is the number of points, T is the trirep class

% 
% By Shengxin Zhu, LCP, IAPCM, May/06/2015
%
if nargin==0
    test_Points2Points();
    return
end
if ~isa(T,'TriRep')
    error('Points2Points:input argument should be a trirep class')
end
E=edges(T);
np=length(T.X);
P2P=sparse(E(:,1),E(:,2),1,np,np);
P2P=P2P+P2P'+speye(np,np);

function PAT=test_Points2Points()
[p,e,t]=initmesh('squareg','Hmax',0.7,'init','off');
%pdemesh(p,e,t) ;
T=t(1:3,:)';
PAT=PointsInSimplex(T);
dtc=TriRep(t(1:3,:)',p');
labelvertices(dtc);
P2P=Points2Points(dtc);
figure,spy(P2P)
