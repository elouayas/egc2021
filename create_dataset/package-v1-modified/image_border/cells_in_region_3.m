%% to decide if cells are within certain region
function in_out = cells_in_region_3(X,C)
%calculate all regions, odd is in and even is out
if size(C, 2)>0 %C exists
%X saves coordinates, C saves contour results
h = C(1,1);
Start = find(C(1,:) == h);
End = [Start(2:end)-1 size(C,2)];
n = length(Start);
m = size(X,1);
in_out = zeros(m,1);
xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));

%% region modification
% add borders
add_border = zeros(n,5);
for i = 1:n
    S = C(:, Start(i)+1:End(i));
    cross_bound = length(find(S(1,:) == xmin |...
        S(1,:) == xmax | S(2,:) == ymin | S(2,:) == ymax));
    % 0 means not crossing boundary
    if cross_bound == 0
        continue;
    else
         cross_y = find(S(1,:) == xmin| S(1,:) == xmax);
         cross_x = find(S(2,:) == ymin | S(2,:) == ymax);
         if ~isempty(cross_y) && ~isempty(cross_x)
             % add one point
             %[x2,y1] = [S(1, cross_y), S(2, cross_x)]
             %x1 = S(1, cross_x);
             y1 = S(2, cross_x);
             x2 = S(1, cross_y);
             %y2 = S(2, cross_y);
             add_border(i,1) = 1;
             add_border(i,2) = x2;
             add_border(i,3) = y1;
         else
             cross_y_1 = find(S(1,:) == xmin);
             cross_y_2 = find(S(1,:) == xmax);
             cross_x_1 = find(S(2,:) == ymin);
             cross_x_2 = find(S(2,:) == ymax);
             if (length(cross_y_1)==2 || length(cross_y_2)==2 ||...
                     length(cross_x_1)==2 || length(cross_x_2)==2)
                 continue;
             else
                 disp('cross two x/y axises detected');
                 add_border(i,1) = 2;
                 if ~isempty(cross_y)
                     %[xmin, y1], [xmax, y2]
                     %[S(1,1),S(2,1)],[S(1,end),S(2,end)]
                     %add [xmin, ymin], [xmax, ymin]
                     add_border(i,2) = S(1,end);
                     add_border(i,3) = ymin;
                     add_border(i,4) = S(1,1);
                     add_border(i,5) = ymin;
                 else
                     %add [xmin, ymin], [xmin, ymax]
                     add_border(i,2) = xmin;
                     add_border(i,3) = S(2,end);
                     add_border(i,4) = xmin;
                     add_border(i,5) = S(2,1);
                 end
                 
             end
         end %end if isempty(cross_x)/isempty(cross_y)
    end
end %end for loop for adding borders

%% determine if cell within the region
for i = 1:m
    o = 0;
    for j = 1:n
        xv = C(1, Start(j)+1:End(j));
        yv = C(2, Start(j)+1:End(j));
        if add_border(j,1)==1
            [xv, yv] = add_border_point(xv,yv,add_border(j,2:3));
        elseif add_border(j,1)==2
            [xv, yv] = add_border_point(xv,yv,add_border(j,2:5));
        end
        %in_out(i) = inpolygon(X(i,1),X(i,2),xv,yv);
        if inpolygon(X(i,1),X(i,2),xv,yv) == 1
           o = o+1;
        end
    end
    if mod(o,2)
        %odd, in
        in_out(i) = 1;
    end
end
else
    m = size(X,1);
    in_out = zeros(m,1);
end