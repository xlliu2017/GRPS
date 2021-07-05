%%Example
% This is a main file
%% Geometry parameters convention
%  D : degeree of freedom
%  T/t £º coarse/fine tesselation
%  E/e :  coarse/fine boundary edge
%  F/f:   coarse/fine internal (or + boundary ?) edges
%  P/p:   coarse/fine vertices
%  S/s:   coarse/fine time steps
%% Convenion for matrices
%  Kd: FVE stiffness matrix for div(a grad )
%  Kd2:  Kd2 = div(a grad ( div(a grad u))
%  K:  P1 FEM stiffness matrix for div(a grad )
%  M:  P1 FEM mass matrix
%  B/b: coarse/fine right handside f
%  C : Constraint (matrix)
%% Convenion for scalars
%  aopt: function handles for the conductivity coefficient
%  fopt: function handle or choice for the right handside
%  N :   global variables used in the initializationof coarse mesh
%  Nc;  Nc x Nc  coarse grid mesh
%  J:   the number of refinements
%  nl:   the number of maximal (ccoarse tesselation) layers will be used for local basis
%% Adjacent matrices
%  D2T: adjacent matrix for degrees of freedom (*coarse* basiss functions)
%       to tesselation/triangl, it can be P2T, T2T, F2T or E2T.
%  T2T: adjacent matrix for a triangle to another
%
%% Sep up paramters
% profile on

aopt=@(x,y)ex11_a(x,y); fopt=@(x,y) sin(x);  nl=8; Btypelist=['M', 'T']; normTypelist = ['a'];% set up parameters
Nc=[4 8 16 ];
J =[5 4 3];      % Caution!! when Btype is 'T', J can't have the number of 1

for normType = normTypelist
    for k=1:length(Nc)       
        tic;
        %% genterating coarse mesh and prepare coarse mesh information
        [pc, ec, tc] = CoarseMesh(Nc(k));                  % generate the coarse mesh
        npc = length(pc); ntc = length(tc); nec = length(ec);
        PinT=adjPinT(tc);                                  % construct D2T
        T2T=Simplex2Simplex(PinT);
        
        %% Refine the mesh J times and prepare fine mesh information
        [pf,ef,tf]=RefineMeshJ2(pc,ec,tc,J(k),Nc(k));             % refine the mesh J times
        npf=length(pf);
        pinT=adjPinT(tf);                                         % adjacent matrix p to T
        p2E=adjp2E(ef,size(pinT,1));
        E2T = adjEinT(ec);
        %% find the solution on the fine mesh
        [uf,Kd,K,Ff,fidx]=finemeshsolve3(pf,tf(1:3,:),aopt,fopt);
        Kd2=Kd*Kd;                                      % Compute the global Kd2
        % figure,H=pdeplot(pf,ef,tf,'xydata',uf,'zdata',uf,'mesh','off','colormap','Cool');
        % set(gca,'fontsize',20)
        toc; ct1 = toc; % time costed by generating mesh, geomety infomation and finesolve
        
        %% all the calculation above has nothing to do with the Btype
        
        for Btype = Btypelist
            tic;
            
            %% Coarse mesh freedome
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
                case 'M'
                    DinT1 = speye(ntc); DinT2 = adjEinT(ec); DinT = [DinT1;DinT2];
                    Didx1 = 1:ntc; Didx2 = freedom4general(DinT2,1:ntc,'E')+ntc; Didx = [Didx1(:);Didx2(:)];
                    nbasis = ntc+length(Didx2);
                case 'D'
                    DinT1 = speye(ntc); DinT2 = adjEinT(ec); DinT3 = PinT; DinT = [DinT1;DinT2;DinT3];
                    Didx1 = 1:ntc; Didx2 = freedom4general(DinT2,1:ntc,'E')+ntc; 
                    Didx3 = freedom4general(PinT,1:ntc,'P');
                    Didx = [Didx1(:);Didx2(:);Didx3(:)+nec+ntc];
                    nDidx1 = ntc; nDidx2 = length(Didx2); nDidx3 = length(Didx3);
                    nbasis = nDidx1+nDidx2+nDidx3;    
            end
            %[P2T,T2T,fdm_coarse,nbasis]=CoarseFreedomP(tc); % information for coarsemesh
            
            %% this is the constraint matrix
            switch Btype
                case 'P'
                    BH=sparse(Didx,Didx,1,npf,npf);
                case 'E'
                    BH=adjE2p(ef,npf)/2;    % edges constraint
                case 'T'
                    BH=pinT';BH=BH/3;
                case 'M'
                    BH2=adjE2p(ef,npf)/2;  BH1=pinT';BH1=BH1/3;
                    BH = [BH1;BH2];
                case 'D'
                    BH2=adjE2p(ef,npf)/2;  BH1=pinT';BH1=BH1/3; BH3=sparse(Didx3,Didx3,1,npc,npf);
                    BH = [BH1;BH2;BH3];
            end
            
            %% initial the storage for error at each level
            l2e=zeros(nl,1);h1e=l2e;lie=l2e;
            
            %% main loop for each layer
            D2T = DinT;
            for layer=[3:8,nl]
                if layer>=2
                    [D2T]=IncreasingLayers2(DinT,T2T,layer);
                end
                psil=ComptLBasis4general(Didx,D2T,DinT,pinT,Kd,Kd2,BH,Btype,normType,ntc,nec,npc,ef,ec,E2T,p2E);
                %
                Kc=psil'*K*psil;       % aseemble the coarse stiffness matrix
                Fc=psil'*Ff;           % assemble the load vectors
                alpha=(Kc\Fc);         %
                ucsol=psil*alpha;      % Solution on the coarse grid.
                %    figure,H=pdeplot(pf,ef,tf,'xydata',ucsol,'zdata',ucsol,'mesh','off','colormap','Cool');
                %    set(gca,'fontsize',20)
                [l2e(layer),h1e(layer),lie(layer)]=error_2d(pf',tf(1:3,:)',ucsol,uf);
            end
            disp([l2e,h1e,lie])
            toc;
            ct2 = toc;    % time costed by caculation basis and the coarse system
            name=sprintf('%s-%s-%sbasis-Nc%dJ%dnl%d-%s.mat',func2str(aopt),normType,Btype,Nc(k),J(k),nl,date)
            save(name,'l2e','h1e','lie','ct1','ct2','fopt')
            fprintf('-------------Done-----------------\n')
            % profile viewer
            % p=profile('info')
            % profsave(p,'Test5_profile_results4c')
        end % for j = 1:length(Btype)
    end % for k=1:length(Nc)
end

%%





