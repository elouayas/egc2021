%% to analyze cancer images with cell labeling
%% 1) get filepath
function main(DataFolder)
%cd('/Users/shidan/Desktop/perl/package/image_border/');
%cd([WorkingDir, '/image_border/']);
%DataFolder = '/Users/shidan/Desktop/perl/package/examples/image01';
fileFolder = fullfile([DataFolder, '/cell_Info']);
featureFolder = [DataFolder, '/extracted_features'];
dirOutput = dir(fullfile(fileFolder,'*.xlsx'));
fileNames = {dirOutput.name}';
number_files = length(fileNames);
% setting parameters
M = 200;
if exist([DataFolder, '/extracted_feature.mat'])
    load([DataFolder, '/extracted_feature.mat']);
else
    extracted_feature = cell(number_files+1,7);
    extracted_feature{1,1} = 'filename';
    extracted_feature{1,2} = 'peri_blue';
    extracted_feature{1,3} = 'peri_red';
    extracted_feature{1,4} = 'peri_green';
    extracted_feature{1,5} = 'size_blue';
    extracted_feature{1,6} = 'size_red';
    extracted_feature{1,7} = 'size_green';
    save([DataFolder, '/extracted_feature.mat'],'extracted_feature');
end

for number_file = 1:number_files %process every file
    %% 2) load file and set parameters
    fprintf('started %d\n',number_file);
    file = fileNames{number_file};
    filename = [fileFolder, '/',file];
    save_name = file(1:end-5);
    disp(save_name);
    mat_name = [featureFolder, '/', save_name, '.mat'];
    if exist(mat_name)
        disp('processed');
    else
    %X = xlsread(filename); not worked on server for xlsx writed by R... 
    %So using readtable(filename) instead.
    X = readtable(filename);
    X = X(:, 2:4);
    X = table2array(X);
    m = size(X,1);
    %% 3) plot_fit_2d and save image
    %save images to '~/cell_border'
    %M = 200;
    %construct data and fit
    %m = size(X,1);
    %save_name = file(1:end-5);
    jpg_name = [DataFolder, '/cell_border/', save_name, '.jpg'];
    
    b = find(X(:,3) == 1); %% for blue cells
    r = find(X(:,3) == 2); %% for red cells
    g = find(X(:,3) == 3); %% for green cells
    w = find_blank_region(X, M); %% add blank region "cells"
    length_w = size(w,1);
    yb = zeros(m+length_w, 1);
    yr = zeros(m+length_w, 1);
    yg = zeros(m+length_w, 1);
    yw = zeros(m+length_w, 1); 
    yb(b) = 1;
    yr(r) = 1;
    yg(g) = 1;
    yw((m+1):end) = 1;
    x = X(:,1:2);
    X = [X; w];
    Yb = fit_2d(X, yb, M);
    Yr = fit_2d(X, yr, M);
    Yg = fit_2d(X, yg, M);
    Yw = fit_2d(X, yw, M);
    [xq,yq] = meshgrid(linspace(min(X(:,1)), max(X(:,1)), M)',...
    linspace(min(X(:,2)), max(X(:,2)), M)');
    if ~exist(jpg_name)
        %plot (right/left inverted)
        disp('plotting cells and borders...');
        figure(1);
        contour(-xq,yq,Yb,[0.5 0.5],'k');
        hold on;
        contour(-xq,yq,Yr,[0.5 0.5],'m');
        contour(-xq,yq,Yg,[0.5 0.5],'g');
        contour(-xq,yq,Yw,[0.5 0.5],'y');
        color = ['b', 'r', 'g'];
        %m = size(X,1);
        for i = 1:m
            plot(-X(i,1), X(i,2), 'color', color(X(i,3)), 'marker', '.');
        end
        hold off;
        saveas(1, jpg_name);
    end
    %% 4) region determination
    Cb = contour(xq,yq,Yb,[0.5 0.5]);
    Cr = contour(xq,yq,Yr,[0.5 0.5]);
    Cg = contour(xq,yq,Yg,[0.5 0.5]);
    Cw = contour(xq,yq,Yw,[0.5 0.5]);
    
    %find cells within a certain region    
    %lympho region
    In_b = cells_in_region_3(X, Cb);
    %interpolation
    tubaobao = interp2(xq,yq,Yb,X(1,1),X(1,2),'linear');
    if(tubaobao >= 0.5)
        tubaobao = 1;
    else
        tubaobao = 0;
    end
    danbaobao = In_b(1);
    if(tubaobao ~= danbaobao)
        L_corr = 1;%indicating corrected
    else
        L_corr = 0;
    end
    
    %stroma region
    In_r = cells_in_region_3(X, Cr);
    %interpolation
    tubaobao = interp2(xq,yq,Yr,X(1,1),X(1,2),'linear');
    if(tubaobao >= 0.5)
        tubaobao = 1;
    else
        tubaobao = 0;
    end
    danbaobao = In_r(1);
    if(tubaobao ~= danbaobao)
        S_corr = 1;
    else
        S_corr = 0;
    end
    
    %cancer region
    In_g = cells_in_region_3(X, Cg);
    %interpolation
    tubaobao = interp2(xq,yq,Yg,X(1,1),X(1,2),'linear');
    if(tubaobao >= 0.5)
        tubaobao = 1;
    else
        tubaobao = 0;
    end
    danbaobao = In_g(1);
    if(tubaobao ~= danbaobao)
        T_corr = 1;
    else
        T_corr = 0;
    end
    %% 5) calculating perimeter
    disp('calculating region perimeter...');
    peri_b = polygon_perimeter(Cb, L_corr, X);
    peri_r = polygon_perimeter(Cr, S_corr, X);
    peri_g = polygon_perimeter(Cg, T_corr, X);
    %% 6) calculating size
    disp('calculating region size...');
    size_b = polygon_size_2(Cb, L_corr, X);
    size_r = polygon_size_2(Cr, S_corr, X);
    size_g = polygon_size_2(Cg, T_corr, X);
    %% 7) save to mat
    disp('writing features...');
    extracted_feature{number_file+1,1} = save_name;
    extracted_feature{number_file+1,2} = peri_b;
    extracted_feature{number_file+1,3} = peri_r;
    extracted_feature{number_file+1,4} = peri_g;
    extracted_feature{number_file+1,5} = size_b;
    extracted_feature{number_file+1,6} = size_r;
    extracted_feature{number_file+1,7} = size_g;
    
    fprintf('finished %d\n',number_file);
    single_file = extracted_feature(number_file+1,:);
    save(mat_name, 'single_file');
    end
    save([DataFolder, '/extracted_feature.mat'],'extracted_feature');
end