function [G,F] = compute_dual_wasserstein(K,gamma,q,g, options)

% compute_dual_wasserstein - compute energy/gradient of dual Wasserstein distance
%
%   [G,F] = compute_dual_wasserstein(K,gamma,q,g);
%
%   Set options.rev_result=1 if you want the output to be [F,G]
%   Set options.flatify=1 if you want the output to be F(:) instead of F. 
%
%   Copyright (c) 2015 Gabriel Peyre


options.null = 0;
rev_result = getoptions(options, 'rev_result', 0);
flatify = getoptions(options, 'flatify', 0);

if size(q,3)>1
    % execute on each %
    G = []; F = [];
    for i=1:size(q,3)
        [G(:,end+1),F(end+1)] = compute_dual_wasserstein(K,gamma,q(:,i),g(:,i), options);
    end
    return;
end

mysum = @(x)sum(x(:));

alpha = exp(g/gamma);
q1 = q./K(alpha);
F = -gamma * mysum( q .* log( q1 ) );
G = alpha .* K( q1 );

if flatify
    G = G(:);
end
if rev_result % reverse the results
    [G,F] = deal(F,G);
end

end