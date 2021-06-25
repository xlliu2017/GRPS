function[p,e,t]=RefineMeshJ2(p,e,t,J,Nc)
%% [p,e,t]=RefinemeshJ(p,e,t,J) refine the meshdata p,e,t J times,
%  and the input and output data have the same name, 
%
% See also:CoarseMesh
%
%
% by Shengxin.Zhu@xjtlu.edu.cn
%
%
% global N
% N = Nc;     % I doubt that if this is necessary, just make sure
[p,e,t]=refinemesh('unitsquare_cg',p,e,t);
for k=2:J
    [p,e,t]=refinemesh('unitsquare_cg',p,e,t);
end
% [p,e,t]=refinemesh('unitsquare_nonuniform_cg',p,e,t);
% for k=2:J
%     [p,e,t]=refinemesh('unitsquare_nonuniform_cg',p,e,t);
% end