function A = jdtest3(x,y,m,contrast)
% EX6_A returns the values of the coeficient a in the equation -\nabla a \nabla u=f
% The original exa_6b depends on the length of x, whent he lenght of x
% changes, the coarse mesh and the fine mesh do not have the same
% conductivity. 
% Usage : a = ex6_a(x, y)
% 
% See also: exa_6b 
% Modifed by zhu_shengxin@iapcm.ac.cn
if nargin==0
    A=test_jdtest3();
    return
end
if nargin==2
m=8;
end
if nargin<=3
    contrast=500;
end

dx=1/m;
rng('default');
R=randi(contrast,m,m);
c = 1 + R;
A = zeros(size(x));
I=ceil(x/dx);I(I==0)=1;
J=ceil(y/dx);J(J==0)=1;
for k=1:numel(x)
    A(k)=c(I(k),J(k));
end

function z=test_jdtest3()
[x,y]=meshgrid(1/64:1/32:1-1/64);
z=jdtest3(x,y);
mesh(x,y,z)
set(gca,'fontsize',20)

end
end
