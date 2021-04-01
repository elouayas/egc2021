%% to add coordinates of "blank cells"
function blank_region = find_blank_region(X, M)
%M: number of meshes
xi = X(:,1)';
yi = X(:,2)';
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

blank_region = [];

for i = 1:M
    for j = 1:M
        dx = xi-xo(i);
        dy = yi-yo(j);
        dist = dx.*dx + dy.*dy;
        weight = exp(-dist/sigma_square);
        if weight < exp(-1)
            blank_region = [blank_region; xo(i), yo(j), 4];
        end
    end
end