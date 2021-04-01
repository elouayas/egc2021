function [xv, yv] = add_border_point(xv, yv, point)
for i = 1:length(point)/2
    xv = [xv, point(2*i-1)];
    yv = [yv, point(2*i)];
end