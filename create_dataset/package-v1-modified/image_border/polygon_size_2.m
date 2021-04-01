%% to calculate polygon size
%for input C, result from contour function
%considered borders and corrections
function poly_size = polygon_size_2(C, r_corr, X)
%r_corr indicated if corrected
if size(C, 2)>0 %C exists
    h = C(1,1);
    Start = find(C(1,:) == h);
    End = [Start(2:end)-1 size(C,2)];
    n = length(Start);
    poly_size = 0;
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
                 %disp('cross two x/y axises detected');
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

%% detecting regions within another region
region_inside = zeros(n,1);
for i = 1:n
     x = C(1, Start(i)+1);
     y = C(2, Start(i)+1);
     num_in = 0;
     temp = [];
     for j = 1:n
        if j == i
            continue;
        end
        xv = C(1, Start(j)+1:End(j));
        yv = C(2, Start(j)+1:End(j));
        if add_border(j,1)==1
            [xv, yv] = add_border_point(xv,yv,add_border(j,2:3));
        elseif add_border(j,1)==2
            [xv, yv] = add_border_point(xv,yv,add_border(j,2:5));
        end
        if inpolygon(x, y, xv, yv) == 1
            num_in = num_in + 1;
            temp(num_in) = j;
            %region_inside(i) = j; %to save the region number which includes this region
            %break;
        end
     end
     if num_in == 1
         region_inside(i) = temp; 
     elseif num_in == 3
         %disp('num_in == 3 found');
         for k = 1:3
             num_in = 0;
             x = C(1, Start(temp(k))+1);
             y = C(2, Start(temp(k))+1);
             for l = 1:3
                 if l == k
                     continue;
                 end
                 xv = C(1, Start(temp(l))+1:End(temp(l)));
                 yv = C(2, Start(temp(l))+1:End(temp(l)));
                 if add_border(temp(l),1)==1
                     [xv, yv] = add_border_point(xv,yv,add_border(temp(l),2:3));
                 elseif add_border(temp(l),1)==2
                     [xv, yv] = add_border_point(xv,yv,add_border(temp(l),2:5));
                 end
                 if inpolygon(x, y, xv, yv) == 1
                     num_in = num_in + 1;
                     %region_inside(i) = j; %to save the region number which includes this region
                     %break;
                 end
             end %end l loop
             if num_in == 2
                 region_inside(i) = temp(k);
                 break;
             end
         end %end k loop
     end %end if num_in == 1/3
end
    %% calculate parameters & borders
    %if ~r_corr %no correction
        for i = 1:n
            S = C(:, Start(i)+1:End(i));
            if add_border(i,1)==1
                [xv,yv] = add_border_point(S(1,:),S(2,:),add_border(i,2:3));
            elseif add_border(i,1)==2
                [xv,yv] = add_border_point(S(1,:),S(2,:),add_border(i,2:5));
            else
                xv = S(1,:);
                yv = S(2,:);
            end
            if region_inside(i)
                poly_size = poly_size - polyarea(xv,yv);
            else %not inside another region
                poly_size = poly_size + polyarea(xv,yv);
            end %end if region_inside
        end %end for loop
        
if r_corr %if with correction
    poly_size = (xmax - xmin)*(ymax - ymin) - poly_size;
end
else
    poly_size = 0;
end