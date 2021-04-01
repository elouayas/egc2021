%% feature extraction
%to extract features of cancer images. E.g, proportion of other cells in
%lymphocytes region; perimeter/size of stroma/lympho region

%% load data and calculate C
filename = '60021_1.mat';
load(filename); %X
M = 200;
%construct data and fit
m = size(X,1);
yb = zeros(m, 1);
yr = zeros(m, 1);
b = find(X(:,3) == 1); %% for blue cells (lympho)
r = find(X(:,3) == 2); %% for red cells  (stroma)
yb(b) = 1;
yr(r) = 1;
x = X(:,1:2);
Yb = fit_2d(x, yb, M);
Yr = fit_2d(x, yr, M);
[xq,yq] = meshgrid(linspace(min(X(:,1)), max(X(:,1)), 200)',...
linspace(min(X(:,2)), max(X(:,2)), 200)');
Cb = contour(xq,yq,Yb,[0.5 0.5]);
Cr = contour(xq,yq,Yr,[0.5 0.5]);
%% find cells within a certain region
%lympho region
In_b = cells_in_region(X, Cb);
In_b_T = find(In_b == 1);
In_b_Lymph = length(find(X(In_b_T,3) == 1));
In_b_Tumor = length(find(X(In_b_T,3) == 3));
In_b_Stroma = length(In_b_T) - In_b_Lymph - In_b_Tumor;
%5170, 66, 418 for lymph before considering circle in other circle

%stroma region
In_r = cells_in_region(X, Cr);
In_r_T = find(In_r == 1);
In_r_Lymph = length(find(X(In_r_T,3) == 1));
In_r_Tumor = length(find(X(In_r_T,3) == 3));
In_r_Stroma = length(In_r_T) - In_r_Lymph - In_r_Tumor;

%verify correct detection
%confirmed, saved as 'cells_in_region_confirmation.jpg'
%contour(xq,yq,Yr,[0.5 0.5],'m');
%hold on;
%color = ['b', 'r'];
%for i = 1:m
%    plot(X(i,1), X(i,2), 'color', color(In_r(i)+1), 'marker', '.');
%end
%hold off;
%% calculating perimeter
peri_b = polygon_perimeter(Cb);
peri_r = polygon_perimeter(Cr);
