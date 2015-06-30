%%
% Computaiton of barycenter of MEG activity.

addpath('toolbox/');
addpath('lbfgsb/Matlab/');
addpath('toolbox_meg/');
addpath('data_meg/');

%%
% Helpers.

dotp = @(x,y)sum(x(:).*y(:));
mynorm = @(x)norm(x(:));
normalize = @(x)x/sum(x(:));
normalize1 = @(X)X ./ repmat(sum(X), [size(X,1) 1]);

%%
% Load MEG dataset 

load('mesh');
load('edges');
XY = vertices(:,1:2)';
E = double(edges'+1);
load('data_condensed');
Xall = { normalize1( data(labels==-1,:)' ), normalize1( data(labels==+1,:)' ) };

nE = size(E,2);
N = size(XY,2);

%%
% Helpers for visualization.

estim = @(x,i)x(i);
prctile = @(v,p)estim(sort(v(:)), round(length(v(:))*p) );
clamp = @(v,a,b)max(min(v,b),a);
clamp_prc = @(v,eta)clamp( v, prctile(v,eta),prctile(v,1-eta) );
eta = .04;
clampa = @(v)rescale( clamp_prc( v, eta ) );
opt.nbrls = 12;
visu = @(p)plot_mne(XY, clampa(p), opt );

x = Xall{1}(:,1);
clf;
visu(x);


%% 
% Global parameters.

% TV Weight
if not(exist('lambda'))   
    lambda = 40; 
    lambda = 0;
end



%% 
% Set up Gibb kernel

gamma = .2;
gamma = .09;
d = sqrt( compute_distance_matrix(XY) );
d = d / median(d(:));
Gibbs = exp( -(d/gamma).^2 );
Kernel = @(x)Gibbs*x;

% visu(Kernel(x));


%%
% Set up the grad/div operators.

% gradient sparse matrix
G = sparse( [1:nE 1:nE]', [E(1,:)'; E(2,:)'], [ones(nE,1); -ones(nE,1)] );
% convert gradient to NxQ matrices, nasty but imposed by the TV-regularization code
Q = ceil(nE/N);
fillGrad = @(g)reshape([g; zeros(Q*N-nE,1)], [N Q]);
options.mygrad = @(f)fillGrad( G*f );
options.mydiv = @(w)-G'*w(1:nE)';

% Projection on l^inf ball of radius lambda.
options.ProjInf = @(x,tau)max(min(x,tau),-tau);

%%
% Save observation.

repout = 'results/meg/input/';
if not(exist(repout))
    mkdir(repout);  
    figure(1);
    for k=1:2
        X = Xall{k};
        for i=1:size(X,2)
            clf;
            visu(X(:,i));
            saveas( gcf, sprintf('%s%d-%d.png', repout, k, i), 'png' );
        end
    end
end

%%
% Output directories. 

rep = 'results/meg/bary/';

if not(exist(rep))
    mkdir(rep);
end

%%
% Mean results


figure(1);
for icl=1:length(Xall)
    X = Xall{icl};
    clf;
    visu(mean(X,2)); drawnow;
    saveas( gcf, sprintf('%sbary%d-mean.png', rep, icl), 'png' );
end

%%
% Run the iterations.

options.niter = 500*8;
options.algorithm = 'fista';
options.algorithm = 'lbfgs';

% kernel .2, lamda=0
options.fista_L = gamma*50;
% kernel .1, lamda=0
options.fista_L = gamma;
% kernel .1, lambda=1
options.fista_L = gamma*1000;

% kernel .2, lambda=.1
options.fista_L = gamma*100; 

for icl=1:length(Xall)
    X = Xall{icl};
        
    K = size(X,2);
    w = ones(K,1)/K;
    
    %%
    % Run algorithm.
    
    [p, v, Fval] = compute_tv_barycenters(w, X, Kernel, gamma, lambda, options);
    Constr = abs(v(1:nE)')/max(lambda,1e-10);
    
    %%
    % Save
    
    [~,i] = max(w);
    p0 = p(:,i);
    figure(1);
    clf;
    visu(p0);
    saveas( gcf, sprintf('%sbary%d-lambda%d.png', rep, icl, round(lambda*1000)), 'png' );
   
    %%
    % Results
    
    figure(1);
    clf; 
    for i=1:min(K,6)
        subplot(2,3,i);
        visu(p(:,i));
    end
     drawnow;
        
    figure(2);
    clf;
    plot(Constr); axis tight; drawnow;
  
    figure(3); clf;
    plot(Fval); 
    axis tight;
    % plot(log10(Fval-min(Fval))); axis tight;
    
    
end