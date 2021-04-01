function border_length = border_length_calculation(S, indi, xmin, xmax, ymin, ymax)
%% find indi = 1 sections
%indi_1 = find(indi == 1);
%number_sections = 0;
i = 1;
S = [S S(:,1)];
indi = [indi indi(:,1)];
border_length = 0;
while i <= length(indi)
    if indi(i)
        Start = i;
        while i<=length(indi)&&indi(i)
            i = i+1;
        end
        if (i-Start) >= 2
            % is a section and can calculate length
            L = S(:, Start:(i-1));
            cross_bound = length(find(L(1,:) == xmin |...
                L(1,:) == xmax | L(2,:) == ymin | L(2,:) == ymax));
            if cross_bound == 3
                L = L(:, 1:(end-1));
            end
            E = L(:,2:end);
            L = L(:,1:(end-1));
            E_L = E-L;
            x = E_L(1,:);
            y = E_L(2,:);
            border_length = border_length + sum(sqrt(x.^2+y.^2));  
        end
    else
        i = i+1;
    end
    %disp(i)
end