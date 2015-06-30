function plot_mne(XY,h, options)

% plot_mne - plot of a vertex-based data for mne
%
%   plot_mne(XY,h, options);
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0;
% interpolation method
meth = getoptions(options, 'interp', 'natural');
% number of level-sel
nbrls = getoptions(options, 'nbrls', 6);
% color map
CM = getoptions(options, 'CM', parula(256));

n = 256;
t = linspace(-.1,1.1,n);
[Y,X] = meshgrid(t,t);

if 0
A = griddata(XY(1,:),XY(2,:), h, X(:),Y(:), meth);
A = reshape(A,[n,n]);
else
F = scatteredInterpolant(XY(1,:)',XY(2,:)', h(:), meth, 'linear');
A = reshape(F(X(:),Y(:)),[n,n]);
end
A = A';

% put a disk mask
N = size(XY,2);
m = mean(XY, 2);
U = XY - repmat(m, [1 N]);
r = 1.05*max(sqrt(sum(U.^2))); 
%
Mask = (X-m(2)).^2 + (Y-m(1)).^2 <= r^2;
Mask3 = repmat(Mask, [1 1 3]);

A(Mask==0)=0;
Acol = apply_colormap(A, CM);
Acol(Mask3==0)=1;

s = linspace(0,1,100);
ellipse = @(m,r)m(1) + r(1)*cos(2*pi*s) ... 
        + 1i * ( m(2) + r(2)*sin(2*pi*s) );

% remove outside disk
B = A;
B(Mask==0) = NaN;

lw = 2;

cval = linspace(0,1,nbrls);
hold on;
imagesc(t,t,Acol);
plot(XY(1,:), XY(2,:), '.k');
contour(t,t,B, cval, 'k', 'LineWidth', lw);
axis tight;
axis image; 
axis off;
% add a circle
plot( ellipse(m, [r r]), 'k', 'LineWidth', lw );
% add nose
q = [.47, .92]; 
h = .08;
plot( [q(1)-h, q(1)], [q(2) q(2)+h], 'k', 'LineWidth', lw );
plot( [q(1)+h, q(1)], [q(2) q(2)+h], 'k', 'LineWidth', lw );
% add hears
h1 = [.03 .09];
plot( ellipse( [m(1)-r-h1(1), m(2)], h1 ), 'k', 'LineWidth', lw );
plot( ellipse( [m(1)+r+h1(1), m(2)], h1 ), 'k', 'LineWidth', lw );
%
axis([m(1)-r-h, m(1)+r+h,m(2)-r,m(2)+r+h]);

end