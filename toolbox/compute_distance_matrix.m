function dist = compute_distance_matrix(X,Y)

% compute_distance_matrix - compute pairwise distance matrix.
%
%   D = compute_distance_matrix(X,Y, metric);
%   We have D(i,j)=|X(:,i)-Y(:,j)|^2.
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<2
    Y = X;
end

if size(X,1)>size(X,2)
    X = X';
end
if size(Y,1)>size(Y,2)
    Y = Y';
end

nX = size(X,2);
nY = size(Y,2);
X2 = sum(X.^2,1);
Y2 = sum(Y.^2,1);
dist = (repmat(Y2,nX,1)+repmat(X2.',1,nY)-2*X.'*Y);