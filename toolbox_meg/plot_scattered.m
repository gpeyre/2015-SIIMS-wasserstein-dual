function plot_scattered(XY, ps, col)

hold on;
if size(XY,1)==2
    for i=1:size(XY,2)
        plot(XY(1,i), XY(2,i), '.', 'MarkerSize', ps(i), 'color', col(:,i) );
    end
elseif size(XY,1)==3
    for i=1:size(XY,2)
        plot3(XY(1,i), XY(2,i), XY(3,i), '.', 'MarkerSize', ps(i), 'color', col(:,i) );
    end
else
    error('Problem');
end

end