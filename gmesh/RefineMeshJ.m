function[p,e,t]=RefineMeshJ(p,e,t,J)
%% [p,e,t]=RefinemeshJ(p,e,t,J) refine the meshdata p,e,t J times,
%  and the input and output data have the same name, 
%
% See also:CoarseMesh
%
%
% by Shengxin.Zhu@xjtlu.edu.cn
%
%
[p,e,t]=refinemesh('unitsquare_cg',p,e,t);
for k=2:J
    [p,e,t]=refinemesh('unitsquare_cg',p,e,t);
end