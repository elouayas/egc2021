%to transfer contour matrix C to another matrix without the decollator
function M = contour_2_matrix(C)
if size(C, 2)>0 %C exists
    M = [];
    h = C(1,1);
    Start = find(C(1,:) == h);
    End = [Start(2:end)-1 size(C,2)];
    n = length(Start);
    for i = 1:n
        M = [M C(:, (Start(i)+1):End(i))];
    end
else
    M = [];
end