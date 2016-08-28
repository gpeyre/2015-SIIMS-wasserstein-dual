a = dir('*.png');
for i=1:length(a)
    s = a(i).name;
    f = imread(s);
    imwrite(255-f, s);
end