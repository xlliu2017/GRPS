function [E2p]=adjE2p(e,npts)
%% function [E2p]=adjp2E(ef) construct the adjacent matrices
% from the mesh data of a coarse edge to the fine vertices
%
%by Shengxin.Zhu@xjtlu.edu.cn
%npts=max(max(e(1:2,:)));                % number of fine vertices
nE=max(e(5,:));                         % number of coarse edges 
E2p=sparse(e(5,:)',e(1,:)',1,nE,npts);
E2p=sparse(e(5,:)',e(2,:)',1,nE,npts)+E2p;   
