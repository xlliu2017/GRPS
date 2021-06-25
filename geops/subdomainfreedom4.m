
function[sdfidx,sdbidx,sdcidx]=subdomainfreedom4(pinT,DinT,D2T,BasisID,...
    Dtype,nDidx1,nDidx2,nDidx3,ef,ec,E2T,p2E)
%%[sdpidx,sdpbidx,CpInFp]=subdomainfreedom(P,E,T,ADJ,BasisID) return the 
%freedom information on a subdomain which consist of the support  of current
%basis support.
%
% See also: PDESDP, ADJMATRIX2
% 
% if nargin==0
%     return
% end

if nargin==0
    test_subdomainfreedom();
    return;
end
% narginchk(5,6)
% if size(P,1)~=2 | size(E,1)~=7 | size(T,1)~=4
%     error('RPS;geops:subdomainfreedom: input argument P,E,T should be the mesh data in pdetoolbox')
% end

idx = find(D2T(BasisID,:));  % ADJ is the adjacent matrix on coarse mesh

% if nargin==5
%     idx = ctnum(idx);      % because the initial mesh ordering of the simplex is different with adjmatrix ordering
% end
% [sdfidx,sdbidx] = freedom5(pinT,idx,'P');
[sdcidx,~] = freedom5(DinT,idx,Dtype);

% p2E=adjp2E(ef,size(pinT,1));
% E2T = adjEinT(ec);
flagE = sum(E2T(:,idx),2);
bidx4E = flagE == 1;
flagp = sum(p2E(:,bidx4E),2);
sdbidx = find(flagp>0);

sdfidx = setdiff(find(sum(pinT(:,idx),2)),sdbidx);

%==============nested function==================================

    function [fidx,bidx] = freedom5(DT,idx,Dtype)
        
        switch Dtype
            case 'P'
                flag = sum(DT(:,idx),2);
                fidx = find(flag == 6);
                bidx = find((flag>=1).*(flag <=3));
            case 'E'
                flag = sum(DT(:,idx),2);
                fidx = find(flag == 2);
                bidx = find(flag == 1);
            case 'T'
                fidx = idx;
                bidx = [];
            case 'M'
                DinT2 = DT(size(DT,2)+1:end,:);
                flag = sum(DinT2(:,idx),2);
                fidx2 = find(flag == 2)+size(DT,2); fidx = [idx(:);fidx2];
                bidx = find(flag == 1);
            case 'D'
                DinT2 = DT(nDidx1+1:nDidx1+nDidx2,:); 
                DinT3 = DT(nDidx1+nDidx2+1:end,:);
                flag2 = sum(DinT2(:,idx),2); flag3 = sum(DinT3(:,idx),2);
                fidx2 = find(flag2 == 2)+nDidx1; 
                fidx3 = find(flag3 == 6)+nDidx1+nDidx2;
                fidx = [idx(:);fidx2;fidx3];
                bidx = find(flag2 == 1);    
        end        
    end
%===============================================================

end % subdomainfreedom2


% SDL=find(idx);
% if nargin==6
% SDL=ctnum(SDL);      % because the initial mesh ordering of the simplex is different with adjmatrix ordering
% end
% T=T';
% tsd=find(ismember(T(:,4),SDL));  % connection between coarse mesh and fine mesh through subdomain index
% [sdpidx,~,sdbidx]=freedom2(T(tsd,1:3));
% cidx=find(sdpidx<=size(ADJ,2));  % size(ADJ,2) means the maximum index of coarse points, the index of coarse points don't change after refinemesh
% CpInFp=sdpidx(cidx); % coarse points in freedom points of subdomain
% 
% % the above code can be simplified by getting rid of find

function test_subdomainfreedom()
 Nc=16;J=2; Btype = 'E';
 global N;
 N=Nc;
% genterating coarse mesh
[pc, ec, tc] = initmesh('unitsquare_cg', 'Hmax', inf); % this will use N




ntc=length(tc);
                                       % coarse simplex number
 
PinT=adjPinT(tc);

T2T=Simplex2Simplex(PinT);% transfer to trirep data


switch Btype
    case 'P'
        DinT = PinT;
        [Didx,~] = freedom4general(DinT,1:ntc,'P'); % here I use Didx to replace fdm_coarse
        nbasis = length(Didx);
    case 'E'
        DinT = adjEinT(ec);  % remain to complete
        [Didx,~] = freedom4general(DinT,1:ntc,'E'); % remain to be modified
        nbasis = length(Didx);
    case 'T'
        DinT = speye(ntc);
        Didx = 1:ntc;
        nbasis = ntc;
end

[pf,ef,tf]=refinemesh('unitsquare_cg',pc,ec,tc);
for k=2:J
    [pf,ef,tf]=refinemesh('unitsquare_cg',pf,ef,tf);
end
pT=adjPinT(tf);
p2E=adjp2E(ef,size(pT,1));
E2T = adjEinT(ec);
D2T = DinT;
D2T=IncreasingLayers2(D2T,T2T,3);

for bid=1:10:100  %length(Didx)/2
    [sdpidx,sdbidx,sdcidx]=subdomainfreedom4(pT,DinT,D2T,Didx(bid),Btype,...
        length(tc),length(ec),length(pc)...
        ,ef,ec,E2T,p2E);
    figure,pdemesh(pc,ec,tc)
    hold on
    switch Btype
        case 'P'
            plot(pf(1,sdpidx),pf(2,sdpidx),'ro','markersize',10)
            plot(pf(1,sdbidx),pf(2,sdbidx),'kx','markersize',10)
            plot(pf(1,sdcidx),pf(2,sdcidx),'gp','markersize',10)
        case 'E'
            plot(pf(1,sdpidx),pf(2,sdpidx),'ro','markersize',10)
            plot(pf(1,sdbidx),pf(2,sdbidx),'kx','markersize',10)
            plot(pf(1,[ec(1,sdcidx);ec(2,sdcidx)]),pf(2,[ec(1,sdcidx);ec(2,sdcidx)]),'LineWidth',5)
        case 'T'
            %plot(pf(1,sdpidx),pf(2,sdpidx),'ro','markersize',10)
            %plot(pf(1,sdbidx),pf(2,sdbidx),'kx','markersize',10)
            flag = sum(PinT(:,sdcidx),2);
            fidx = find(flag >= 1);
            plot(pf(1,fidx),pf(2,fidx),'gp','markersize',10)
            centerP = find(PinT(:,Didx(bid)));
             plot(pc(1,centerP),pc(2,centerP),'ro','markersize',15)
             centerP_f = find(pT(:,Didx(bid)));
             plot(pf(1,centerP_f),pf(2,centerP_f),'bo','markersize',5)
    end
    % xlim([0, 0.6])
    % ylim([0,0.6])
    title('Red(o) inner freedom, black(x) boundary points','fontsize',12)
end
end