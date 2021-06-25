function[pc,tc,dt]=pt2dt(p,t)
%pt2dt  data structure transform function, 
%   [pc,tc,dt]=pt2dt(p,t) transfer mesh data from initialmesh or refinemesh 
%   to the data structure used in the TriRep class 
%
%   see also TriRep, initmesh
% 

%  Last modified by Shengxin.Zhu@xjtlu.edu.cn
%
narginchk(2,2)
if size(p,1)~=2
    error('RPS:geops:pt2dt: 1st input argument is invalid')
end
if size(t,1)~=4
      error('RPS:geops:pt2dt: 2nd input argument is invalid')
end
pc=p';
tc=t(1:3,:);
dt=TriRep(tc',pc(:,1),pc(:,2));
