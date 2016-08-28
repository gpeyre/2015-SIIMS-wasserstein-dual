function X = load_densities(p,n)

% load_densities - load some images`
%
%   X = load_densities(p,n);
%
%   p=#images
%   n=width
%   X(:,:,i) is an image.
%
%   Copyright (c) 2015 Gabriel Peyre

if not(exist('imresize'))
    imresize = @(x,s)image_resize(x,s);    
end

X = zeros(n,n,p);
for i=1:p
    f = imread(sprintf('shape%dfilled.png',i));
    f = rescale( double(sum(f,3)) );
    padding = abs(size(f,2)-size(f,1));
    pad1 = floor(padding/2);
    pad2 = padding-pad1;    
    if size(f,1)<size(f,2)
        f = [zeros(pad1,size(f,2)) ; f ; zeros(pad2,size(f,2))];
    elseif size(f,2) < size(f,1)
        f = [zeros(size(f,1),pad1) f zeros(size(f,1),pad2)];
    end    
    f = 1-imresize(f,[n n]);
    f(f<0) = 0; f = f>.01;    
    f = f + 1e-3; f = f/sum(f(:));
    X(:,:,i) = f;
end

end