function [fidx,flag,bidx]=freedom3(T)
%%[fidx,flag,bidx]=freedom2(T) return the freedoms of the vertices which  is the
%(sub)-dmain consisits of the simplices T, which is a nt*3 arrary. 
%      flag(i):     is the number of simplices which includes the vertex i.
%                   flag can be used as the weights of area integral in the
%                   subdomain.  
%                   For this specific triangularization, if flag(i)<=3,
%                   then the ith vertex is on the boundary of the
%                   subdomain.
%      fidx:        is the index of thoses indices in the (sub)-doma
%      bidx:        is the boundary points
%USAGE: 
%         [p, t,dt] = squaregeom(12);
%         [adj,T2T] = adjmatrix0(dt,'P',4);
%         vnode=adj(10,:);
%         idx=find(vnode);
%         labelvertices(dt,idx);
%         [ff,flag]=freedom2(t(idx,:));
%
%See also:  ADJMATRIX2, 
%By Shengxin Zhu, IAPCM, MAY/2015
%
%Reference: 
if nargin==0
    test_freedom3();
    return;
end
[n1,n2]=size(T);
if n2~=3
    if n1==3
        T=T';
    else
    error('RPS:geops:subedges:input argument is invvalid')
    end
end
nm=max(max(T));
nt=length(T);
%flag=zeros(nm,1);
tt=1:nt;
pT=sparse(T(:,1),tt' ,1,nm,nt);
pT=sparse(T(:,2),tt',1,nm,nt)+pT;
pT=sparse(T(:,3), tt',1, nm,nt)+pT;
flag=sum(pT,2);
fidx=find(flag>3);     % only works for the regular triangularization
if nargout==3
    bidx=find( (flag>=1).*(flag <=3));
end

%%
function test_freedom3()
[p, t,dt] = squaregeom(32);
[adj,~] = adjmatrix2(dt,'P',4);
np=length(p);
% labelvertices(dt)
for i = 1:floor(np/2):np
     vnode=adj(i,:);
     idx=find(vnode);
     labelvertices(dt,idx);
     tic;[~,~]=freedom3(t(idx,:));t3=toc
     tic;[~,~]=freedom2(t(idx,:));t2=toc
end



