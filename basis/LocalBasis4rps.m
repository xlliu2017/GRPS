function [ulocal]=LocalBasis4rps(Kd, A, B, centerDidx, sdpidx,sdcidx,sdbidx)
%% ulocal=computebasis(D,D2,BC,basisID, subfidx, subbidx,CpInFp)
% returns the basis function 
% with the freedom subfidx, where subidx is the index of the local basis in
% the whole domain
%Aims: to solve the local basis function of the form
%      
%       | A     B^T  |    u                |0 |
%       |            |                  =  |   |
%       | B      0   |   lambda      | e_i|
%
% INPUT:    D: the divagrad operator in the whole domain 
%                 D2: D2=D*D for the whole domain
%                 F: is the assemble 
%               subfidx: the index of the fine points/edges in the subdomain in the whole 
%                              fine mesh index
%               subcidx: the index of the fine points/edges in the subdomain in the 
%                              whole coarse points/mesh index
%                 BC: the constraint matrix
%
%OUTPUT: ulocal  is the lcoal basis function
%

Eb=zeros(length(A),1);
A=A(sdpidx,sdpidx);
M=Kd(sdpidx,sdbidx)*Kd(sdbidx,sdpidx);  %corector for A
A=A-M;
B=B(sdcidx,sdpidx);

Eb(centerDidx)=1;
Eb=Eb(sdcidx);
ab=A\B';
ulocal=ab*((B*ab)\Eb);


