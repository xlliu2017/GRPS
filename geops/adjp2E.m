function [p2E]=adjp2E(e,npts)
%% function [p2E]=adjp2E(e) construct the adjacent matrices between the fine
% veterices and the coarse edges. 
%  See also: subdomainfreedomE2
%
%
% 
%
%
nE=max(e(5,:));         % number of coarse edges 
p2E=sparse(e(1,:)',e(5,:)',1,npts,nE);
p2E=sparse(e(2,:)',e(5,:)',1,npts,nE)+p2E;
