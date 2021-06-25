function[psil]=ComptLBasis4general(Didx,D2T,DinT,pinT,Kd,Kd2,BH,Btype,normType,ntc,nec,npc,ef,ec,E2T,p2E)
%% [psil]=ComptLBasis(fdm_coarse,D2T,p2T,type) compute the local
% basis function of type, where type={'P','E','V'} are the type of the
% basis function
% Input: fdm_coarse: index of the coarse basis funcion in the
%                    coarse mesh data
%        D2T;    adjacent matrix of the basis function of
%        p2T:    adjacent matrix of the fine vertices to the tesselation
% Output: psil:  a sparse matrix to store the numerical value of the local
%                basis function
%
% See also: subdomainfreedom2, LocalBasis2
%
%

nbasis=length(Didx);
npf=size(pinT,1);

i4sparse = [];
j4sparse = [];
k4sparse = [];

switch normType
    case 'a'
        parfor bid=1:nbasis
            [sdpidx,~,sdcidx]=subdomainfreedom4(pinT,DinT,D2T,Didx(bid),Btype,ntc,nec,npc,ef,ec,E2T,p2E);
            ulocal = LocalBasis4generalNorm(Kd,BH,Didx(bid),sdpidx,sdcidx);
            i4sparse = [i4sparse;sdpidx];
            j4sparse = [j4sparse;bid*ones(size(sdpidx))];
            k4sparse = [k4sparse;ulocal];
        end
        
        
    case 'r'
        parfor bid=1:nbasis
            [sdpidx,sdbidx,sdcidx]=subdomainfreedom4(pinT,DinT,D2T,Didx(bid),Btype,ntc,nec,npc);
            ulocal = LocalBasis4rps(Kd,Kd2,BH,Didx(bid),sdpidx,sdcidx,sdbidx);
            i4sparse = [i4sparse;sdpidx];
            j4sparse = [j4sparse;bid*ones(size(sdpidx))];
            k4sparse = [k4sparse;ulocal];
        end
end

psil = sparse(i4sparse,j4sparse,k4sparse,npf,nbasis);

