function [G,F] = compute_fista_grad(U,w, X, Kernel, gamma, options)

% compute_fista_grad - compute the dual energy of barycenter
%
%   [G,F] = compute_fista_grad(U,w, X, gamma,  options)
%
%   Notations:
%   U(:,:,1:K-1) are u's
%   U(:,:,K:K+1) is v
%
%   Compute the gradient G of 
%       F = sum_i w_i H_{X(:,:,i)}^*( U(:,:,i)/w_i );
%   where it is assumed that
%       U(:,:,K) = div(v)-sum_i U(:,:,i)
%
%   Set options.rev_result=1 if you want the output to be [F,G]
%   Set options.flatify=1 if you want the output to be F(:) instead of F. 
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0;
rev_result = getoptions(options, 'rev_result', 0);
flatify = getoptions(options, 'flatify', 0);

% l^inf projection
mygrad = getoptions(options, 'mygrad', [], 1);
mydiv = getoptions(options, 'mydiv', [], 1);


global Fval; 

K = length(w);
N = size(U,1); 
% number of dimension of the gradient
Q = size(mygrad(zeros(N,1)),2);

u_e = @(U)U(:,1:K-1);
v_e = @(U)U(:,K:end);

u = u_e(U);
v = v_e(U);

w1 = repmat( reshape(w(1:end-1),[1 K-1]) , [N 1]);
u(:,end+1) = mydiv(v)/w(end) - sum( u .* w1/w(end), 2);

F = 0; G = [];
for i=1:K
    [g,f] = compute_dual_wasserstein(Kernel,gamma, X(:,i), u(:,i) );
    G(:,i) = g;
    F = F + w(i) * f;
end
Fval(end+1) = F;

g = G(:,end);
G = G(:,1:end-1);
G = G - repmat(g, [1 size(G,2)]);
G = G .* w1; 
G(:,end+1:end+Q) = -mygrad(g);

if flatify
    G = G(:);
end
if rev_result % reverse the results
    [G,F] = deal(F,G);
end
    
end