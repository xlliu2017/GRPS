function[a]=spe10(x,y,layer)
%% a=spe10(x,y) return the a at the center of (x,y), where(x,y) in the 
% domain [0,220] x [0 60]. layer is an optional inputargument
% 
% Try:  spe10
%
% By Shengxin.Zhu@xjtlu.edu.cn, 2017, July 8th.
%
if nargin==0
    test_spe10()
    return
end
if nargin==2
    layer=1;
else
    if (layer>86) || (layer < 1)
        error('rps/ex_a/spe10: the 3ard argument should be an integer between 1 and 85')
    end
 x= 100*x; y = 100*y;    % rescale the geometry 
[ny,nx]=size(x);
if ny*nx<60*220
    error('The finemesh should be a refined mesh from 60*220')
end
load('K1.mat')
%K1=K1+4;
A=10.^K1(:,:,layer);
a=zeros(numel(x),1);
y(y==0)=eps;                        % to make sure ceil(y) ~=0;
x(x==0)=eps;                        % to make sure ceil(x) ~=0;
for k=1:length(x(:));
    a(k)=A(min(ceil(y(k)),60),min(ceil(x(k)),220));   % the  index may out of the range
end
end
%%
function test_spe10()
[x,y]=meshgrid(0.25:0.5:220-0.25, 0.25:0.5:60-0.25);
for layer=1:3:85-3
a=log(spe10(x,y,layer));
figure,axis equal, set(gca,'fontsize',20)
subplot(3,1,1),contour(x,y,reshape(a,size(x)))
str1=sprintf('layer %d',layer);
title(str1)
a=spe10(x,y,layer+1);a=log(a);
subplot(3,1,2),contour(x,y,reshape(a,size(x)))
str1=sprintf('layer %d',layer+1);
title(str1)
a=spe10(x,y,layer+2);a=log(a);
subplot(3,1,3),contour(x,y,reshape(a,size(x)))
str1=sprintf('layer %d',layer+2);
title(str1)

end
