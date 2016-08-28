function [x, R, z] = GeneralizedForwardBackward( gradF, proxGi, z, nIter, ga, la, verbose, report )
%
%	 [x, R, z] = GeneralizedForwardBackward( gradF, proxGi, z, nIter, ga, [la=1], [verbose=false], [report=None] )
%
% Solve
%		min_x F(x)+\sum_{i=1}^{n} Gi(x)
% where
%	F is differentiable and its gradient is 1/be-Lipschitz continuous,
%	prox_{ga*Gi}(x) = argmin_y 1/2||x-y||^2 + ga*Gi(x) is computable;
% with generalized forward-backward.
% INPUT:
%	'gradF': function handle: class(x) -> class(x)
%		computes the gradient of F
%		can be set to empty-array if F = 0, resulting in a relaxed version of Douglas-Rachford algorithm
%	'proxGi': n-long cell array of function handles : (class(x),double) -> class(x)
%		proxGi{i}(x,ga) computes prox_{ga*Gi}(x)
%	'z': [size(x)]-by-n matrix
%		the initial auxiliary variables
%	'nIter': non-negative integer
%		the number of iterations recquired
%	'ga': double
%		the size of forward step,
%		ga \in ]0,2*be[
%	'la' [default=1]: double
%		relaxation of the iterations, la \in ]0,min(2/3, 2/(1+(2be/ga))[
%	'verbose' [default=false]: logical
%		set the diaplay of iterations
%	'report' [default=None]: function handle : class(x) -> 1-by-m matrix
%		a user-defined report called at each iteration
% OUTPUT:
%	'x': vector or N-D matrix
%		the final minimizer
%	'R' [default='empty']: nIter-by-m matrix
%		the report sequence
%	'z': [size(x)]-by-n matrix
%		the final auxiliary variables
%
% Hugo Raguet 2011

if nargin < 6, la = 1; end
if nargin < 7, verbose = false; end
doReport = nargin >= 8;

n = length( proxGi );
catDim = ndims( z );
x = mean( z, catDim );
N = numel( x );
zi = zeros( size( x ) );
R = [];

for it=1:nIter
	if verbose, progressbar(it,nIter); end
	if doReport, R(it,:) = report( x ); end
	if isempty( gradF )
		forward = 0;
	else
		forward = ga*gradF(x);
	end
	for i=1:n
		idx = (i-1)*N+1:i*N;
		zi(:) = z(idx);	
		z(idx) = zi + la*( proxGi{i}( 2*x - zi - forward, n*ga ) - x );
	end
	x = mean( z, catDim );
end

end %GeneralizedForwardBackward
