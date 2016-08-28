%%
% Computation of gradient flow.

addpath('toolbox/');
addpath('data_images/');
addpath('lbfgsb/Matlab/');

%% 
% Global parameters.

name = 'twodisks';
name = 'randdisks';
name = 'randdots';
name = 'hibiscus';

% Width of the images. 
n = 255;


% TV Weight
lambda = .5; % very strong
lambda = .05;
lambda = .1;

% randdisks
lambda = 10;
lambda = 4;
lambda = 2; 

% hibiscus
lambda = .2; 

%%
% Helpers.

dotp = @(x,y)sum(x(:).*y(:));
mynorm = @(x)norm(x(:));
normalize = @(x)x/sum(x(:));
rifft2 = @(x)real(ifft2(x));
flat = @(x)reshape(x, [n*n size(x,3)]);
resh = @(x)reshape(x, [n n size(x,2)]);

%%
% Load data.

Xmin = 1e-3;
X0 = load_image(name, n);
X0 = crop(X0,n);
X0 = rescale(sum(X0,3))+Xmin; 
X0 = X0/sum(X0(:));

%% 
% Set up blur

gamma = 2; 
gamma = 1/4;  % image
gamma = 1; 
blur = load_filtering('imgaussian', n);
Kernel0 = @(x)blur(x,gamma);
Kernel = @(x)flat(Kernel0(resh(x)));

%%
% Set up the operators.

flatgrad = @(v)[flat(v(:,:,1)) flat(v(:,:,2))];
reshgrad = @(w)reshape( w, [n n 2] );
options.mygrad = @(f)flatgrad( grad(resh(f)) );
options.mydiv = @(w)flat(div( reshgrad(w) ));

% Projection on l^inf ball of radius lambda.
Amp = @(u)sqrt(sum(u.^2,2));
options.ProjInf = @(u,lambda) u .* repmat( min( 1, lambda./max(Amp(u),1e-10) ), [1 2]);

%%
% Generate weights. 


rep = ['results/images/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Run the iterations.

options.algorithm = 'fb';
options.algorithm = 'fista';
options.niter = 500;
options.fista_L = gamma/100;
options.fista_L = gamma/500;

niter = 20;

X = X0;
imwrite( rescale(X), sprintf('%s%s-flow-%d.png', rep, name, 0), 'png' );
for it=1:niter
    
    %%
    % Run algorithm.
    
    X1 = X;
    [X, v, Fval] = compute_tv_barycenters(1, flat(X1), Kernel, gamma, lambda, options);
    X = resh(X);
    v = resh(v);
    Constr = sqrt(sum(v.^2,3))/lambda;
            
    %%
    % Save
    
    imwrite( rescale(X), sprintf('%s%s-flow-%d.png', rep, name, it), 'png' );
    
    %%
    % Results
    
    figure(1);
    clf; 
    imageplot(X); drawnow;
    
    figure(2);
    imagesc(Constr);
    colormap parula(256);
    colorbar;  drawnow;
    
    figure(3);
    plot(Fval); axis tight;
    plot(log10(Fval-min(Fval))); axis tight;  drawnow;
    
end