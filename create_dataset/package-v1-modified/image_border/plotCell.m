%% to plot cells of three types
function plotCell(X)
color = ['b', 'r', 'g'];
m = size(X,1);
for i = 1:m
    plot(X(i,1), X(i,2), 'color', color(X(i,3)), 'marker', '.');
    hold on;
end
hold off;
