function T2T=Simplex2Simplex(P2T)
%%T2T=Simplex2Simplex(P2T) Return a sparse matrix of nt*nt T2T(i,j)=1 if
%%Simplex i and Simplex j have a common vertex
% See also: adjPinT
%
% By Shengxin Zhu, May/07/2017
if nargin==0
    test_Simplex2Simplex();
    return
end
if issparse(P2T)
  %  [np,nt]=size(P2T);
    T2T=P2T'*P2T;
    T2T=bsxfun(@gt,T2T,0);
end

function T2T=test_Simplex2Simplex()
[p,~,t]=initmesh('squareg','Hmax',0.7,'init','off');
%pdemesh(p,e,t) ;
P2T=adjPinT(t);
T2T=Simplex2Simplex(P2T);
dtc=TriRep(t(1:3,:)',p');
labelvertices(dtc)
