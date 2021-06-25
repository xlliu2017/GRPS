function adj=IncreasingLayers2(D2T,T2T,l)
%%%adj=IncreasingLayer(D2T,T2T,l) return a adjacency matrix of the degree
%%%of freedom to the tesselation, D2T(i,j)~=0 implies the degree of freedom
% i's (vertex, edge or triangle) to the triangle in the mesh.
% Input: D2T
%        T2T
% Output adj
% 
% See also:twomesh_jdtest5_v17
% 
% Modified by Shengxin Zhu, June 2017
if nargin==0
    test_IncreasingLayers2();
    return
end
narginchk(3,3);
[~,nt]=size(D2T);
nt2=size(T2T,1);
if nt~=nt2
    error('RPS:geops:INcreasingLayers:Inputargument is invalid');
end
D2T=double(D2T);
T2T=double(T2T);
for k=1:l-1
    D2T=D2T*T2T;
end
adj=bsxfun(@gt,D2T,0);


function PAT=test_IncreasingLayers2()
[p,e,t]=initmesh('squareg','Hmax',0.4,'init','off');
%pdemesh(p,e,t) ;
PAT=adjPinT(t);
T=t(1:3,:)';
PAT=PointsInSimplex(T);
T2T=Simplex2Simplex(PAT);
dtc=TriRep(t(1:3,:)',p');
labelvertices(dtc)
D2T=IncreasingLayers2(PAT,T2T,3);
for node=[4,6,21,31]
    vnode=D2T(node,:);
    plotneighbour(dtc,vnode);
end
Nc=6;
[pc, tc, dt] = squaregeom(Nc);   % initial geometry
labelvertices(dt)
PAT=PointsInSimplex(tc);
T2T=Simplex2Simplex(PAT);
D2T=IncreasingLayers2(PAT,T2T,3);
for node=[1 4 7 19 25]
    vnode=D2T(node,:);
    plotneighbour(dt,vnode);
end
E2T=EdgesInSimplex(dt);
labelEdges(dt);
D2T=IncreasingLayers2(E2T,T2T,3);
for EN=[1, 6, 19, 66, 67];
    Enode=D2T(EN,:);
    plotneighbour(dt,Enode);
end



