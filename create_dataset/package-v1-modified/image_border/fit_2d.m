%% written by tutu
function output = fit_2d(X, y, M)
%M: number of meshes
xi = X(:,1)';
yi = X(:,2)';
zi = y';
xmin = min(xi);ymin = min(yi);
xmax = max(xi);ymax = max(yi);
sigma_square = 4*(xmax-xmin)*(ymax-ymin)/size(X,1);
%sigma_square = 888; 
%sigma_square = 2*(xmax-xmin)*(ymax-ymin)/size(X,1);
%887.8272 for image 60021_1
%sigma_square = (xmax-xmin)*(ymax-ymin)/size(X,1);
if ~exist('M', 'var') || isempty(M)
    M = round(sqrt(size(input,1)));
end
%M = round(sqrt(size(input,1)));
xo = linspace(xmin,xmax,M);
yo = linspace(ymin,ymax,M);

output = zeros(M,M);

for i = 1:M
for j = 1:M
dx = xi-xo(i);
dy = yi-yo(j);
dist = dx.*dx + dy.*dy;
weight = exp(-dist/sigma_square);
weight = weight/sum(weight);
output(i,j) = sum(weight.*zi);
end
end
output = output';
