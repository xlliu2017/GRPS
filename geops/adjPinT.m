function [adj]= adjPinT(t)
% [adj,T2T]=adjPInT(tc) returns the adjacent matrix for vertices to simplex 
% Usage: the adjacent matrix can be the coarse vertices to the coarse
% tessellation or the fine vertices tot he fine mesh tesselation. 
%       
%
%  See Also: INITMESH
%
%
%  by Shengxin Zhu, LCP,IAPCM, MAY/2015/

% Reference:
if nargin == 0
    test_adjPInT();
    return
end
if size(t,1)~=4
    error('RPS:geops:adjPInT: input argument  should be a mesh data ')
end
nt=max(t(4,:));
npts=max(max(t(1:3,:)));
adj=spalloc(npts,nt,10*nt);
for i=1:3
    adj=adj+sparse(t(i,:),t(4,:),1,npts,nt);
end

%% ------------test functions
function adj=test_adjPInT()
Nc=8;J=2;aopt=@ex6_a;fopt=1;  
global N;
N=Nc;
% genterating coarse mesh
[pc, ec, tc] = initmesh('unitsquare_cg', 'Hmax', inf); % this will use N
ctnum=tc(4,:)';                                                              % coarse simplex number
[~,~,dtc]=pt2dt(pc,tc);                                                   % transfer to trirep data
labelvertices(dtc);
[MNT]=adjPinT(tc);                                  % find neigbours simplex of each points
