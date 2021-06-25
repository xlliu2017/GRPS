function [E2T,T2T]=adjE2T(ec,tc)
%% function [E2T,T2T]=adjE2T(pc,ec,tc) construct the adjacent matrices
% from the mesh data
% 
%  
% by zhu_shengxin@iapcm.ac.cn
et=ec(5:7,:)';
et(:,2:3)=et(:,2:3)+1;
ne=length(ec);
nt=length(tc)+1;
adj=sparse(et(:,1),et(:,2),1,ne,nt);
adj=adj+sparse(et(:,1),et(:,3),1,ne,nt);
E2T=adj(:,2:end);
nt=length(tc);
npf=max(max(tc(1:3,:)));
P2T=spalloc(npf,nt,nt*3);
for k=1:3
P2T=P2T+sparse(tc(k,:),tc(4,:),1,npf,nt);
end
T2T=Simplex2Simplex(P2T);
