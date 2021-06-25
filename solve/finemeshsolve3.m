function [uf,D,K,Ff,fidx,M]=finemeshsolve3(pf,tf,aopt,fopt)
%%[uf,D,arf,a,K,Ff,fidx]=finemeshsolve2(pf,tf,aopt,fopt)s solve the probelm
% on the fine mesh .
%INPUT: pf is P in the meshdata 
%       tf is T in the meshdata 
%      aopt is a function handle for the conductivity
%      fopt is a function handle for the lefthandside
% Output
%     uf is the solution on the fine mesh 
%     D  is the stiffness matrix by the div a(grad u) by FVM
%     K  is the stiffness matrix by the FEM
%     Ff is the load vector on the right handside 
%     fidx is the index of degree of freedom of the vertives ( points in
%     the domain
%     arf 
%
%  
%
% Version modified from finemeshsolve2 by Shengxin.Zhu@xjtlu.edu.cn 
% commented by Shengxin.Zhu@xjtlu.edu.cn
narginchk(4,4)
pf=pf';
[nt, n3] = size(tf);    % nt is the number of triangulars
if n3~=3
    if nt~=3
        error('finesolve3: argument 3 is invalid')
    else
        tf=tf';
        nt=n3;                  % number of tiangulars
    end
end         
[xb,yb]=barycenter(pf,tf);
a=feval(aopt,xb,yb);
f=feval(fopt,xb,yb);
%[a,f]= cal_conductivityf(aopt,pf, tf,fopt);
[fidx]=freedom3(tf);           % This should be a finefreedom function by edge transform 
D = divagrad2(pf,tf,a);               % D is the matrix for div a grad
[K,Ff,M] = assemblematrixf(pf, tf, a,f);
Kf=K(fidx,fidx);
b=Ff(fidx);
npf=length(pf);
if npf <=2000
    uf=Kf\b;
end
if npf>2000
    L=ichol(Kf);
    [uf,flag]=pcg(Kf,b,.5/npf,npf,L,L');
    if flag
        switch flag
            case 1
                warning('RPS:solve:finemeshsolve2:pcg does not converge within maxmium iteration')
            case 2
                warning('PBS:solve:finemshsolve2: preconditioners for pcg was ill-conditioned')
            case 3
                warning('RPS:solve:finemeshsolve2:pcg stagnated')
            case 4
                warning('RPS:solve:finemeshsolve2: one of the scalar quantities calculated during pcg became too samll or too large')
        end
    end
end
u=zeros(npf,1);
u(fidx)=uf;
uf=u;
end
function [xb,yb]=barycenter(p,t)
%%
% p should be a nx2 array
% t should be a nt x 3 array 
%
xp=p(:,1);             % x-coordinate
yp=p(:,2);             % y-coordinate
xb=sum(xp(t),2)/3;  % barcenter coordinate
yb=sum(yp(t),2)/3;  % barcenter coordinate
%a=feval(funca,xb,yb);       % conductivity
%f=feval(funcf,xb,yb);       % function 
end