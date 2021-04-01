%% to calculate polygon perimeter
%for input C, result from contour function
%considered borders and corrections
function perimeter = polygon_perimeter(C, r_corr, X)
%r_corr indicated if corrected
if size(C, 2)>0 %C exists
    h = C(1,1);
    Start = find(C(1,:) == h);
    End = [Start(2:end)-1 size(C,2)];
    n = length(Start);
    perimeter = 0;
    border = 0;
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

    %% calculate perimeters & borders
    %if ~r_corr %no correction
        for i = 1:n
            S = C(:, Start(i)+1:End(i));
            E = [S(:,2:end) S(:,1)];
            E_S = E-S;
            x = E_S(1,:);
            y = E_S(2,:);
            cross_bound = length(find(S(1,:) == xmin |...
                S(1,:) == xmax | S(2,:) == ymin | S(2,:) == ymax));
            % 0 means not crossing boundary
            if ~cross_bound
                perimeter = sum(sqrt(x.^2+y.^2)) + perimeter;
            else
                if region_inside(i)
                    %if inside another region
                    cross_y = find(S(1,:) == xmin| S(1,:) == xmax);
                    cross_x = find(S(2,:) == ymin | S(2,:) == ymax);
                    if ~isempty(cross_y) && ~isempty(cross_x)
                        %cross x & cross y
                        %[x1,y1] = [S(1, cross_x), S(2, cross_x)]
                        %[x2,y2] = [S(1, cross_y), S(2, cross_y)]
                        %[x2,y1] = [S(1, cross_y), S(2, cross_x)]
                        x1 = S(1, cross_x);
                        y1 = S(2, cross_x);
                        x2 = S(1, cross_y);
                        y2 = S(2, cross_y);
                        border_temp = abs(x1-x2) + abs(y2 - y1);
                        %perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            %- border - sqrt((x1-x2)^2 + (y1-y2)^2);
                        perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            - sqrt((x1-x2)^2 + (y1-y2)^2);
                        border = border - border_temp;
                    elseif ~isempty(cross_y)
                        %cross y axis
                        perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            - abs(S(2, cross_y(1)) - S(2, cross_y(2)));
                        border = border - abs(S(2, cross_y(1)) - S(2, cross_y(2)));
                    else
                        %cross x axis
                        perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            - abs(S(1, cross_x(1)) - S(1, cross_x(2)));
                        border = border - abs(S(1, cross_x(1)) - S(1, cross_x(2)));
                    end
                else
                    %if not inside another region
                    cross_y = find(S(1,:) == xmin| S(1,:) == xmax);
                    cross_x = find(S(2,:) == ymin | S(2,:) == ymax);
                    if ~isempty(cross_y) && ~isempty(cross_x)
                        %cross x & cross y
                        x1 = S(1, cross_x);
                        y1 = S(2, cross_x);
                        x2 = S(1, cross_y);
                        y2 = S(2, cross_y);
                        border_temp = abs(x1-x2) + abs(y2 - y1);
                        perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            - sqrt((x1-x2)^2 + (y1-y2)^2);
                        border = border + border_temp;
                    elseif ~isempty(cross_y)
                        %cross y axis
                        perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            - abs(S(2, cross_y(1)) - S(2, cross_y(2)));
                        border = border + abs(S(2, cross_y(1)) - S(2, cross_y(2)));
                    else
                        %cross x axis
                        perimeter = sum(sqrt(x.^2+y.^2)) + perimeter...
                            - abs(S(1, cross_x(1)) - S(1, cross_x(2)));
                        border = border + abs(S(1, cross_x(1)) - S(1, cross_x(2)));
                    end
                end
            end %end if cross bound               
        end %end for loop
if ~r_corr %no correction
    perimeter = perimeter + border;
else %if r_corr %corrected
    perimeter = perimeter + 2*(ymax-ymin+xmax-xmin) - border;
end
else
    perimeter = 0;
end
