function A = ex6_a(x, y)
% EX6_A returns the values of the coeficient a in the equation -\nabla a \nabla u=f
% The original exa_6b depends on the length of x, whent he lenght of x
% changes, the coarse mesh and the fine mesh do not have the same
% conductivity. 
% Usage : a = ex6_a(x, y)
% 
% See also: exa_6b   this file have a bug related to the index 
% Modifed by zhu_shengxin@iapcm.ac.cn

m=8;
dx=1/8;
rng('default');         % to guarantee repeatable numerical results 
R=rand(m,m);
c = 1 + 100*(R<0.5);
A = zeros(size(x));
n=length(x);
for k=1:n
    i=floor(x(k)/dx)+1;
    j=floor(y(k)/dx)+1;
    A(k)=c(i,j);
end