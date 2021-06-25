function [ib e]=findbdry(p, t)
% Subroutine to find the boundary nodes, namely, all
% the nodes with value 1 in the correpsonding row in the 
% adjacency matrix.
% It is also possible to find boundary edges, i.e.,
% find all the pairs with value 1 in A.
%
% Fri, 18 Nov, 2011, Lei Zhang


if nargin == 0
    test_findbdry();
    return
end
warning on   % modifiedy by zhu_shengxin@iapcm.ac.cn, bugs 
trep = TriRep(t, p);
%trep=t;
e = freeBoundary(trep)';

ib = indbdry(e);
end

function test_findbdry()
[p, t, dt] = squaregeom(60);
trep = TriRep(t, p);
[ib, e] = findbdry(p, t);
triplot(trep);
hold on ; plot(p(e, 1), p(e, 2), '-r', 'LineWidth',2) ; hold off ;
end

% function test_findbdry()
% x = rand(10,1);
% y = rand(10,1);
% dt = DelaunayTri(x,y);
% p = dt.X;
% t = dt.Triangulation;
% np = size(p,1);
% nt = size(t,1);
% % initialize adjacency matrix for nodes
% A = sparse(np,np);
% tmp = sparse(t(:,1),t(:,2),1,np,np); 
% A = A + tmp + tmp';
% tmp = sparse(t(:,2),t(:,3),1,np,np); 
% A = A + tmp + tmp';
% tmp = sparse(t(:,3),t(:,1),1,np,np); 
% A = A + tmp + tmp';
% A = prod(A-1,2);
% 
% triplot(dt);
% for i=1:np
%     text(x(i),y(i),sprintf('%d',i));
%     if A(i)==0
%         text(x(i),y(i)-0.02,'X')          
%     end
% end
% 
% x1=x(t(:,1));y1=y(t(:,1));  
% x2=x(t(:,2));y2=y(t(:,2));
% x3=x(t(:,3));y3=y(t(:,3));
% xm=(x1+x2+x3)/3;ym=(y1+y2+y3)/3;
% 
% for i=1:nt
%     text(xm(i),ym(i),sprintf('%d',i));
% end

% function ib=findbdry(p,t)
% 
% if nargin == 0
%     test_findbdry();
%     return
% end
% 
% np = size(p,1);
% nt = size(t,1);
% % initialize adjacency matrix for nodes
% A = sparse(np,np);
% tmp = sparse(t(:,1),t(:,2),1,np,np); 
% A = A + tmp + tmp';
% tmp = sparse(t(:,2),t(:,3),1,np,np); 
% A = A + tmp + tmp';
% tmp = sparse(t(:,3),t(:,1),1,np,np); 
% A = A + tmp + tmp';
% 
% A = prod(A-1,2);
% ib = find(A==0);