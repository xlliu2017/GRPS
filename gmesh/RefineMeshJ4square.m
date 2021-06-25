function[p,e,t]=RefineMeshJ4square(p,e,t,J,Nc)
%% [p,e,t]=RefinemeshJ(p,e,t,J) refine the meshdata p,e,t J times,
%  and the input and output data have the same name, 
%
% See also:CoarseMesh
%
%
% by Shengxin.Zhu@xjtlu.edu.cn
%
%
global N
N = Nc;
[p,e,t]=refinemesh('square_cg',p,e,t);
for k=2:J
    [p,e,t]=refinemesh('square_cg',p,e,t);
end