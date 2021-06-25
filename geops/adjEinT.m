
function EinT = adjEinT(e)

ne=max(e(5,:));
nt=max(max(e(6:7,:)))+1;    %plus subdomain 0
EinT=spalloc(ne,nt,2*ne);
for i=6:7
    EinT=EinT+sparse(e(5,:),e(i,:)+1,1,ne,nt); % because there's subdomain 0, which can't be the index of matrix
    
end
EinT = EinT(:,2:end);
end