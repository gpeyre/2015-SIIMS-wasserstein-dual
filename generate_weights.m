function W = generate_weights(p,q)

% generate_weights - generate a bunch of barycentric weights
%
%   W = generate_weights(p,q);
%
%   p in {1 2 3 4} is the simplex dimension
%   q controls number of weights
%
%   Copyright (c) 2015 Gabriel Peyre

switch p
    case 1
        W = 1;
    case 2
        % displacement interpolation
        t = linspace(0,1,q);
        W = [t;1-t];
    case 3
        % triangle interpolation
        W = [ ...
            [0, 0, 1]; ...
            [1, 0, 3]; [0, 1, 3]; ...
            [1,0,1]; [1,1,2]; [0,1,1]; ...
            [3,0,1]; [2,1,1]; [1,2,1]; [0,3,1]; ...
            [1,0,0]; [3,1,0]; [1,1,0]; [1,3,0]; [0,1,0] ...
            ]';
    case 4
        % bilinear interpolation
        t = linspace(0,1,q);
        [T,S] = meshgrid(t,t); S = S(:); T = T(:);
        W = [(1-S).*(1-T) S.*(1-T) (1-S).*T S.*T]';
end

end