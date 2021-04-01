%% to plot C from contour function
% only plot one height
function plot_C(C, color)
h = C(1,1);
interval = find(C(1,:) == h);
if ~exist('color', 'var') || isempty(color)
    color = 'b';
end
for i = 1:length(interval)-1
    x = C(1, (interval(i)+1):(interval(i+1)-1));
    y = C(2, (interval(i)+1):(interval(i+1)-1));
    plot(x,y,color);
    hold on;
end
i = i+1;
x = C(1, interval(i)+1:end);
y = C(2, interval(i)+1:end);
plot(x,y,color);
hold off;