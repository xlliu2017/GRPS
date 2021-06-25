function [x,y]=rectangle_cg(bs,s)
%RECTANGLE_CG   Gives geometry data for the coarse geometry of rectangleg
%                 PDE model
%
%   NE=SQUAREG gives the number of boundary segment
%
%   D=SQUAREG(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left hand region.
%   Row 4 contains the number of the right hand region.
%
%   [X,Y]=SQUAREG(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.


% Copyright 1994-2003 The MathWorks, Inc.
% $Revision: 1.8.4.2 $

global N
N2 = N; N1 = 220/60*N;
% global pc ec tc fc;

% if nargin == 1
%     if strcmp(bs, 'test')
%         test_unitsquare_cg;
%         return
%     end
% end
rng('default') %% set default rand seed in poimesh
[pc, ec, tc] = poimesh('rectangleg',N1,N2);
dtc = TriRep(tc(1:3, :)', pc(1,:)', pc(2,:)');
fc = edges(dtc)';
ti = edgeAttachments(dtc,fc(1,:)',fc(2,:)');

nbs=N1*(N2+1)+N2*(N1+1)+N1*N2;%N*(3*N+2);

if nargin==0,
    x=nbs; % number of boundary segments
    return
end

d = zeros(4, nbs);
for ic = 1:nbs
    d(1, ic) = 0; % start parameter value
    d(2, ic) = 1; % end parameter value
    % d(3, ic) is left hand region
    % d(4, ic) is right hand region
    F1 = pc(:, fc(2, ic)) - pc(:, fc(1, ic));
    if length(ti{ic}) == 1
        thirdpoint = setdiff(tc(1:3, ti{ic}), fc(:, ic));
        F2 = pc(:, thirdpoint) - pc(:, fc(1, ic));
        if ccangle(F1, F2) < pi
            d(3, ic) = ti{ic};
            d(4, ic) = 0;
        else
            d(3, ic) = 0;
            d(4, ic) = ti{ic};
        end
    else
        thirdpoint = setdiff(tc(1:3, ti{ic}(1)), fc(:, ic));
        F2 = pc(:, thirdpoint) - pc(:, fc(1, ic));
        if ccangle(F1, F2) < pi
            d(3, ic) = ti{ic}(1);
            d(4, ic) = ti{ic}(2);
        else
            d(3, ic) = ti{ic}(2);
            d(4, ic) = ti{ic}(1);
        end
    end
end

bs1=bs(:)';

% if isempty(s)
%     x = d;
%     y = [];
%     return
% end

if find(bs1<1 | bs1>nbs),
    error(message('pde:squareg:InvalidBs'))
end

if nargin==1,
    x=d(:,bs1);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
    error(message('pde:squareg:SizeBs'));
end

if ~isempty(s)
    for i = 1:nbs
        % boundary segment i
        ii = find(bs==i);
        if length(ii)
            x(ii)=interp1([d(1,i),d(2,i)],[pc(1,fc(1,i)) pc(1, fc(2,i))],s(ii));
            y(ii)=interp1([d(1,i),d(2,i)],[pc(2,fc(1,i)) pc(2, fc(2,i))],s(ii));
        end
    end
end
    
end
    

function angle = ccangle(v1, v2)
% CCANGLE gives counterclockwise angle from vector v1 to vector v2
angle = mod(atan2(v1(1)*v2(2)-v2(1)*v1(2),v1(1)*v2(1)+v1(2)*v2(2)),2*pi); % Range: 0 to 2*pi radians
end

function test_unitsquare_cg

global N;
global pc ec tc fc;

N = 4;
[d, ~] = rectangle_cg('test', []);

pdemesh(pc, ec, tc);

for i = 1:size(pc, 2)
    text(pc(1, i), pc(2, i), num2str(i));
end

for i = 1:size(fc, 2)
    text((pc(1,fc(1,i))+pc(1,fc(2,i)))/2, (pc(2,fc(1,i))+pc(2,fc(2,i)))/2, ...
        num2str(i));
end

xm = (pc(1, tc(1, :)) + pc(1, tc(2, :)) + pc(1, tc(3, :)))/3;
ym = (pc(2, tc(1, :)) + pc(2, tc(2, :)) + pc(2, tc(3, :)))/3;

for i = 1:size(tc, 2)
    text(xm(i), ym(i), num2str(i));
end

end