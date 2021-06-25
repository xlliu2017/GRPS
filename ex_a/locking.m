function u=locking(x,y,delta)
%% function A=jdtest1(x,y,zl,zh) generate a domain on the unit 
%square with two kinds of material. Two square with center (0.25,0.75)
% (0.75,0.25) with length 0.25 with zh conductivity, other place with zl
% conductivity 
if nargin==0
    A=test_jdtest1();
    return
end
if nargin<=2
  delta=100000;
end
u=sin(2*pi*x).*exp(-2*pi.*y./sqrt(delta));
end

%%
function z=test_jdtest1()
delta=1e5;
[x,y]=meshgrid(1/64:1/32:1-1/64);
z=locking(x,y,delta);
meshc(x,y,z)
set(gca,'fontsize',20)

end

    
