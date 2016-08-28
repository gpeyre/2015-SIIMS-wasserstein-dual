%%
% Computation of barycenter of shapes.

addpath('toolbox/');
addpath('data_shapes/');
addpath('lbfgsb/Matlab/');

%% 
% Global parameters.

% Width of the images. 
n = 99;

% number of input densities
K = 2; 

% type of TV regularization
tvmode = 'iso';
tvmode = 'aniso';

% TV Weight
if not(exist('lambda'))    
    lambda = 20;
    lambda = 40;
end


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

X = load_densities(K,n);

%% 
% Set up blur

gamma = 2; 
gamma = 1; 
blur = load_filtering('imgaussian', n);
Kernel0 = @(x)blur(x,gamma);
Kernel = @(x)flat(Kernel0(resh(x)));
Kernel1 = @(X) cell2mat(cellfun(Kernel, mat2cell(X,n,n,ones(K,1)), 'UniformOutput', false));


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

q = 5; 
W = generate_weights(K,q);

rep = sprintf('results/shapes/bary%d-%s/lambda%d/', K, tvmode, round(lambda*10));

ext = '';
if exist('Wfixed') && not(isempty(Wfixed))
    W = Wfixed;
    rep = sprintf('results/shapes/bary%d-%s/lambda-varying/', K, tvmode);
    ext = ['-lambda' num2str(round(lambda*10))];
end

if not(exist(rep))
    mkdir(rep);
end

%%
% Run the iterations.

switch tvmode
    case 'iso'
        options.algorithm = 'fista';
        options.niter = 5000;
        options.niter = 5000*2;
    case 'aniso'
        options.algorithm = 'lbfgs';
        options.niter = 3000;
        options.niter = 500;
end

% W = ones(K,1)/K;

for iw=1:size(W, 2)
    w = normalize( W(:,iw) );
    
    %%
    % Run algorithm.
    
    [p, v, Fval] = compute_tv_barycenters(w, flat(X), Kernel, gamma, lambda, options);
    p = resh(p);
    v = resh(v);
    
    switch tvmode
        case 'iso'
            Const = sqrt(sum(v.^2,3))/lambda;
        case 'aniso'
            Constr = max(abs(v), [], 3)/lambda;
    end 
    
    %%
    % Save
    
    [~,i] = max(w);
    p0 = p(:,:,i);
    imwrite( rescale(-p0), sprintf('%sbary%d%s.png', rep, iw, ext), 'png' );
    
    %%
    % Results
    
    figure(1);
    clf; warning off;
    imageplot(mat2cell(p, n, n, ones(K,1))); warning on;
    
    figure(2);
    imagesc(Constr);
    colorbar;
    
    figure(3);
    plot(Fval); axis tight;
    plot(log10(Fval-min(Fval))); axis tight;
    
end