function[pc,ec,tc]=CoarseMesh4nonuniform(Nc)
%%[]=CoarseMesh(N) generate the N x N mesh on the unitsquare which is
%%simplely the model a mesh in the followin type 
%      ---------------
%      |/|/|......|/|
%      |/|/|......|/|
% where pc,ec,tc are P, E, T in the meshdata 
%
% See also : initmesh
%
% Example:   [pc,ec,tc]=CoarseMesh(8);
%            figure,lable_cmesh(pc,ec,tc);
% 
%By Shengxin.Zhu@xjtlu.edu.cn
%
%

% Tricks: 
% This function wraps the global variable N in this function, so the scope
% of the global variable N range only this function.
%
% global N; 
% global Hmax;    
global j;     j=log2(Nc)-2;             % this variable wase used in the initmesh
% N=Nc; Hmax = 0.33; 
% [pc, ec, tc] = initmesh('unitsquare_cg', 'Hmax', inf); % this will use N

[pc, ec, tc] = initmesh('unitsquare_nonuniform_cg', 'Hmax', inf);