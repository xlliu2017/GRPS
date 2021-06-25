function A=jdtest1(x,y,zl,zh)
%% function A=jdtest1(x,y,zl,zh) generate a domain on the unit 
%square with two kinds of material. Two square with center (0.25,0.75)
% (0.75,0.25) with length 0.25 with zh conductivity, other place with zl
% conductivity 
if nargin==0
    A=test_jdtest1();
    return
end
if nargin<=2
    zh=100;
    zl=1;
end
A=ones(size(x));
A=zl*A;
mask=find( (abs(x-0.25) <0.125) .* (abs(y-0.75) <0.125) );
A(mask)=zh;
mask=find( ((x-0.75).^2+(y-0.25).^2) <1/64);
A(mask)=zh;
end

%%
function z=test_jdtest1()
zh=100;zl=1;
[x,y]=meshgrid(1/64:1/32:1-1/64);
z=jdtest1(x,y,zl,zh);
meshc(x,y,z)
set(gca,'fontsize',20)

end

    
