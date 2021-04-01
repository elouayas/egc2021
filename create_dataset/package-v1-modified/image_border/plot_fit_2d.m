% to plot fit_2d
%load data
%filename = '60021_1.mat';
filename = '60021_2.mat';

%% construct data
load(filename); %X
M = 200;
%construct data and fit
m = size(X,1);
yb = zeros(m, 1);
yr = zeros(m, 1);
b = find(X(:,3) == 1); %% for blue cells
r = find(X(:,3) == 2); %% for red cells
yb(b) = 1;
yr(r) = 1;
X = X(:,1:2);
Yb = fit_2d(X, yb, M);
Yr = fit_2d(X, yr, M);
%% plot
[xq,yq] = meshgrid(linspace(min(X(:,1)), max(X(:,1)), 200)',...
linspace(min(X(:,2)), max(X(:,2)), 200)');
figure(1);
contour(xq,yq,Yb,[0.5 0.5],'k');
hold on;
contour(xq,yq,Yr,[0.5 0.5],'m');
load(filename); %X
color = ['b', 'r', 'g'];
m = size(X,1);
i = 1;
plot(X(i,1), X(i,2), 'color', color(X(i,3)), 'marker', '.');
hold on;
for i = 2:m
    plot(X(i,1), X(i,2), 'color', color(X(i,3)), 'marker', '.');
end
hold off;
%% plot for left-right inverted
[xq,yq] = meshgrid(linspace(min(X(:,1)), max(X(:,1)), 200)',...
linspace(min(X(:,2)), max(X(:,2)), 200)');
figure(2);
contour(-xq,yq,Yb,[0.5 0.5],'k');
hold on;
contour(-xq,yq,Yr,[0.5 0.5],'m');
load(filename); %X
color = ['b', 'r', 'g'];
m = size(X,1);
i = 1;
plot(-X(i,1), X(i,2), 'color', color(X(i,3)), 'marker', '.');
hold on;
for i = 2:m
    plot(-X(i,1), X(i,2), 'color', color(X(i,3)), 'marker', '.');
end
hold off;