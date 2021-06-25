function [C, F ]= cal_conductivityf(options,points, triangulars,optionf)
% CAL_CONDUCTIVITYF returns the coeficient of the conductivity and the right
% handside at given points and triangulars.
% 
% INPUTS:
%   options: can be a numerical value :0,1,2 or a function handle
%   optionf: can be a defualt value between 0 to 9 or a function handle
%   points and triangulars are the vetices and triangulars from a
%   triangularization 
% OUTPUT:
%   C is a matrix of length NT*1, or NT*2 or NT*3, corresponds with the 
%   coeficient to be a scalar function, diagonal matrix or 2 by 2 matrix
%   F is the load vector
% 
% Shengxin Zhu, May 2, 2015
%
if nargin==0
    C=test_cal_conductivityf();
    return;
end
narginchk(3,4);                  % check number of input arguments
[np, d2]= size(points);          % number of vertices 
if d2~=2
    error('cal_conductivity: argument 2 i invalid')
end
[nt, n3] = size(triangulars);    % nt is the number of triangulars
if n3~=3
    if nt~=3
        error('cal_conductivity: argument 3 is invalid')
    else
        triangulars=triangulars';
        nt=n3;                  % number of tiangulars
    end
end         
xp=points(:,1);                 % 
yp=points(:,2);                 %
x=sum(xp(triangulars),2)/3;     %  x=(x1+x2+x3)/3
y=sum(yp(triangulars),2)/3;     %  y=(y1+y2+y3)/3
if ~isempty(options)
C=ones(size(x));                %  preallocate C
if isnumeric(options)
    switch options
        case 0
            C=ones(size(x))*2;
        case  1
            C=1.5+.01*sin(x);
        case 2
            n=16;
            m=round(sqrt(n)/2);
            g=0.3+2*rand(m,1);
            xmin=min(x);
            xmax=max(x)+1e-8;
            dx=(xmax-xmin)/m;
            for i=1:n
                j=floor((x(i)-xmin)/dx)+1;
                C(i)=g(j);
            end
    end
end
if isa(options,'function_handle')
    C=feval(options,x,y);
end
else
    C=[];
end
if nargin==3
    optionf=0;
end
if isnumeric(optionf)
    switch optionf
         case 0
             F=zeros(size(x));
        case 1
            F=ones(size(x));
        case 2
            F=sin(pi*x).*sin(pi*y);
        case 3
            t=1;
            F=sin(2.4*x-1.8*y+2*pi*t); %% shoudl be 2*pi*t
        case 4   % is there someting incorrect on sigma 
            sigma=.05; 
            F=1/sqrt(2*pi)/sigma*exp(-1/2/sigma^2*(x.^2+(y-0.15).^2));
        case 5
            f0=30;sigma=.05;t=1;
            F=1/sqrt(2*pi)/sigma*exp(-1/2/sigma^2*(x.^2+(y-0.15).^2))...
            .*(1-2*pi^2*(f0*t-1).^2).*exp(-pi^2*(f0*t-1).^2);
        case 6 
            sigma=.05; F=1/sqrt(2*pi)/sigma*exp(-1/2/sigma^2*(x.^2+y.^2));
            u1=0;t=1;
            for k=1:10  
                u1 = u1 + 2*(1-(-1)^k)/(k*pi)*sin(k*pi*t/.5);  % ugly formula  
            end        
            u2=erfc((t-0.5)*8);
            F=F*u1*u2;
    end
end
if nargin==4 && isa(optionf,'function_handle')
    F=feval(optionf,x,y);
end
function C=test_cal_conductivityf()
    Nf = 16 + 1; hf = 1/Nf;
    [Xf,Yf]=meshgrid(linspace(0,1,Nf),linspace(0,1,Nf));
    xf = Xf(:); yf = Yf(:); npf = length(xf);
    x0 = xf + sin(pi/6)*yf; y0 = cos(pi/6)*yf;
    DTf = DelaunayTri(x0,y0);
    triplot(DTf)
    triangulars = DTf.Triangulation;
    points= [xf,yf];
    for prob=0:2
        C=cal_conductivityf(prob,points,triangulars);
        if isnumeric(C);
            fprintf('Test for aseembling conductivity problem %d passed.............\n',prob)
        end
    end
    prob1=@(x,y) 1.5+0.01.*sin(x);
    prob2=@(x,y)  [0.1*sin(x)+1, .2*sin(y)+2];
    C=cal_conductivityf(prob1,points,triangulars);
    if isnumeric(C);
        fprintf('Test for aseembling conductivity function handles passed......\n')
    end
    C=cal_conductivityf(prob2,points,triangulars);
    if isnumeric(C);
        fprintf('Test for assembling conductivity tensor coefcient passed......\n')
    end
    for probf=0:6
        [~, F]=cal_conductivityf(0,points,triangulars,probf);
        fprintf('Test for assembling RHS problem %d  passed....\n',probf)
    end
    probf1=@(x,y) sin(pi*x).*sin(pi*y);
    probf2=@(x,y) sin(2.4*x-1.8*y+2*pi*1);
    [~, F]=cal_conductivityf([],points,triangulars,probf1);
    fprintf('Test for assembling RHS problem %d function handle passed....\n',probf)
    [~, F]=cal_conductivityf(0,points,triangulars,probf2);
    fprintf('Test for assembling RHS problem %d function handle passed....\n',probf)
    
end
end