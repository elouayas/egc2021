%to plot Contour matrix with indicator
function plot_contour_matrix(C, indicator)
switch nargin
    case 2
        h = C(1,1);
        Start = find(C(1,:) == h);
        End = [Start(2:end)-1 size(C,2)];
        n = length(Start);
        figure(1);
        hold on;
        color = ['b','r'];
        for i = 1:n
            l = End(i) - Start(i); %length of this region
            for j = 1:l
                plot(C(1,Start(i)+j),C(2,Start(i)+j),...
                    'color', color(indicator(Start(i)+j)+1),...
                    'marker','.');
            end
        end
        hold off;
    case 1
        h = C(1,1);
        Start = find(C(1,:) == h);
        End = [Start(2:end)-1 size(C,2)];
        n = length(Start);
        figure(1);
        hold on;
        for i = 1:n
                plot(C(1,(Start(i)+1):End(i)),C(2,Start(i)+1:End(i)), 'b');
        end
        hold off;
end