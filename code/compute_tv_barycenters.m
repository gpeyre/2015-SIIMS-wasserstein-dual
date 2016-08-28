function [p, v, Fval] = compute_tv_barycenters(w, X, Kernel, gamma, lambda, options)

% compute_tv_barycenters - OT barycenters with TV regularization
%
%   [p, Constr, Fval] = compute_tv_barycenters(w, X, Kernel, gamma, lambda, options)
%
%   Kernel is the kernel, with bandwidth gamma
%   lambda is the TV weight (=0 for no regularization)
%   X(:,i) is the ith input density
%   w is the vector of weights for barycenters
%   options.algorithm in {'fista' 'bfgs' 'gfb'} controls the used algorithm
%   options.niter controls number of iterations.
%
%   Copyright (c) 2015 Gabriel Peyre


%% 
% Sort weights.
[w,I] = sort(w, 'ascend');
X = X(:,I); 

options.null = 0;
algorithm = getoptions(options, 'algorithm', 'fista');
niter = getoptions(options, 'niter', 5000);
ProjInf = getoptions(options, 'ProjInf', @(x,lamba)max(min(x,lambda),-lambda));
mygrad = getoptions(options, 'mygrad', [], 1);
mydiv = getoptions(options, 'mydiv', [], 1);

N = size(X,1);
K = size(X,2);

% number of dimension of the gradient
Q = size(mygrad(zeros(N,1)),2);

global Fval;  
Fval = [];

u_e = @(U)U(:,1:K-1);
v_e = @(U)U(:,K:end);

switch algorithm
    case {'fista' 'fb'}
        ProxG = @(U,tau)cat(2, u_e(U), ProjInf(v_e(U),lambda) );
        options.rev_result = 0;
        options.flatify = 0;
        GradF = @(U)compute_barycenter_grad(U, w, X, Kernel, gamma,  options);
        L = getoptions(options, 'fista_L', gamma/50);
        opts.method = algorithm;
        opts.niter = niter;
        U0 = zeros(N,K+Q-1);
        [U,R] = perform_fb(U0, ProxG, GradF, L, opts); 
        
    case 'lbfgs'
        resh = @(U)reshape(U, [N K+Q-1]);
        flat = @(U)U(:);
        options.rev_result = 1;
        options.flatify = 1;
        GradF = @(U)compute_barycenter_grad( resh(U), w, X, Kernel, gamma, options);  
        
        if 0
        % CHECK GRADIENT %
        [e,deriv,deriv_fd] = check_gradient(GradF, N*(K+Q-1), 10, 1e-8, 1e-1);
        % plot([deriv;deriv_fd]'); axis tight;
        plot(e); 
        GradWS = @(x)compute_dual_wasserstein(Kernel,gamma, X(:,1), x, options );
        [e,deriv,deriv_fd] = check_gradient(GradWS, N, 50, 1e-12, 1e-5);
        % plot([deriv;deriv_fd]');
        plot(e); 
        end
        
        %%
        % bounds
        u  = flat( cat( 2, Inf+zeros(N,K-1), lambda*ones(N,Q) ) );
        l  = -u;
        % Request very high accuracy for this test:
        opts    = struct( 'factr', 1e4, 'pgtol', 1e-8); % precision high
        opts    = struct( 'factr', 1e2, 'pgtol', 1e-15); % precision high
        opts.m = 10; % memory
        opts.maxIts = niter;
        opts.printEvery  = Inf;
        % opts.errFcn = @(x)rand;
        % Run the algorithm:
        [U, f, info] = lbfgsb(GradF, l, u, opts );
        U = resh(U);
        Fval = info.err(:,1);       

    otherwise
        error('Unknown method.');
        
end

u = u_e(U); v = v_e(U);


% bsxfun(@times,v,areaWeights)
w1 = repmat( reshape(w(1:end-1),[1 K-1]) , [N 1]);
u(:,end+1) = mydiv(v)/w(end) - sum( u .* w1/w(end), 2);

p = [];
for i=1:K
    p(:,i) = compute_dual_wasserstein(Kernel,gamma, X(:,i), u(:,i) );
end

%%
% Re-order results
p(:,I) = p;

end