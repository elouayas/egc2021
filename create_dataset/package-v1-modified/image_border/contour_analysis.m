%% contour analysis
%% area size/perimeter histogram
%after plot_fit_2d.m
%C = contour(xq,yq,Yr,[0.5 0.5],'k');
h = C(1,1);
interval = find(C(1,:) == h);
interval2 = [interval(2:end) size(C,2)+1];
perimeter = interval2 - interval-1;
%hist(perimeter, 0:2:max(perimeter));
%xlim([0,max(perimeter)]);
%% remove circles with length < 8
l_8 = find(perimeter < 8);
C_new = [];
j = 1;
for i = 1:size(perimeter,2)
    if i == l_8(j)
        if j<length(l_8)
            j = j+1;
        end
    else
        C_new = [C_new C(:,interval(i):(interval(i)+perimeter(i)))];
    end
end
figure(2);
% clabel(C_new) %no function for plotting C
plot_C(C,'k');
hold on;
plot_C(C_new, 'b');