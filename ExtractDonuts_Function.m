%% Extracts output (images, pixel values, mean spectra, mean sope and pixel slope vlaues for each bin of given plant)
% Meant for files where center of plant was measured from cropped images
% and not from original uncropped image

%%% Inputs
% NDVI, BW_plant, data_norm_t, h, Wavelengths = Output from from readHyperspecFile3.m function
% M = matrix with center coordinates for plants A, B and C of each image file processed
% name = orig file name to be used for naming new files (date, genotype, stress, etc.)
% date = date of imaging in mmddyyyy format as named in folder and images
% plant = A, B or C
% folder_data = folder where all the donut output will be stored. The folders Imaging/ Mean Slopes/ Mean SPectra and Pixel Values must be contained in the folder specified by this path 
% donuts = number of donuts desired
%%% Outputs
% DonutPixelValues = cell matrix with all pixel values for each donut for given plant. 
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}
    
% DonutMeanSlopes = Average slope of spectra for all pixels in given donut. 
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}
    
% DonutAllSlopes = Slope of spectra for all pixels in given donut.  
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}
    
% DonutMeanSpectra =  Average spectra for all pixels in given donut.  
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}

% xx = wavelength breaks for DonutMeanSlopes
 
function [DonutPixelValues, DonutMeanSlopes, DonutAllSlopes, DonutMeanSpectra, xx] = ExtractDonuts_Function(NDVI, BW_plant, data_norm_t, h, Wavelengths, M, name, date, plant, folder_data)
    
    donuts = 10; %10 leaf segments
    
    labBW = labelmatrix(h);

    s = regionprops(labBW, 'BoundingBox');
    num = length(s);
    BoundingBox = zeros(num, 1);

    for i = 1: num
        BoundingBox(i,1) = s(i).BoundingBox(:,1);
    end

    idx_A = find(BoundingBox <= 500 & BoundingBox >= 100);
    
    idx_C = find(BoundingBox >= 938);

    bw_A = ismember(labBW, idx_A);
    
    bw_C = ismember(labBW, idx_C);

    h_A = bwconncomp(bw_A);
   
    h_C = bwconncomp(bw_C);

    % create matrix with coordinates matching those of image
    dims = size(BW_plant);
    [x,y] = meshgrid(1:dims(2), 1:dims(1));

    %% plant 
   
    h_p = eval(strcat("h_", plant));
    xcent = strcat("plant", plant, "x");
    ycent = strcat("plant", plant, "y");

    cent = zeros(2,1);
    %cent = regionprops(h_p, 'centroid');
    %cent = cat(1, cent.Centroid);
    
    
    cent(1) = eval(xcent); 
    cent(2) = eval(ycent); 

    % extrema? (top-left and bottom-right seem fine, but could use others)
    ex = regionprops(h_p, 'extrema');
    ex = cat(1, ex.Extrema);
    ex_tl = ex(1, :);
    ex_br = ex(5, :);

    % find distance between center and two leaf extrema, continue with longer
    % leaf

    leaf_length1 = sqrt((ex_br(1)-cent(1))^2+(ex_br(2)-cent(2))^2);
    leaf_length2 = sqrt((ex_tl(1)-cent(1))^2+(ex_tl(2)-cent(2))^2);
    
    
    if leaf_length1 < leaf_length2
        x1 = linspace(ex_tl(1), cent(1), 11);
        y1 = linspace(ex_tl(2), cent(2), 11);
        extrema = ex_tl;

   
    else
        x1 = linspace(ex_br(1), cent(1), 11);
        y1 = linspace(ex_br(2), cent(2), 11);
        extrema = ex_br;
    end
    

    %% calculate radii based off these coordinates (third plant)
    rad_1 = sqrt((x1(1)-cent(1))^2+(y1(1)-cent(2))^2);
    rad_2 = sqrt((x1(2)-cent(1))^2+(y1(2)-cent(2))^2);
    rad_3 = sqrt((x1(3)-cent(1))^2+(y1(3)-cent(2))^2);
    rad_4 = sqrt((x1(4)-cent(1))^2+(y1(4)-cent(2))^2);
    rad_5 = sqrt((x1(5)-cent(1))^2+(y1(5)-cent(2))^2);
    rad_6 = sqrt((x1(6)-cent(1))^2+(y1(6)-cent(2))^2);
    rad_7 = sqrt((x1(7)-cent(1))^2+(y1(7)-cent(2))^2);
    rad_8 = sqrt((x1(8)-cent(1))^2+(y1(8)-cent(2))^2);
    rad_9 = sqrt((x1(9)-cent(1))^2+(y1(9)-cent(2))^2);
    rad_10 = sqrt((x1(10)-cent(1))^2+(y1(10)-cent(2))^2);

    x0=cent(1);
    y0=cent(2);

    % Rings between `r1` and `r2`
    f = @(r1,r2) (x-x0).^2+(y-y0).^2<=r2^2 & ... 
        (x-x0).^2+(y-y0).^2>=r1^2;

    % binary mask for each section
    R_1 = logical(f(0,rad_10));
    R_2 = logical(f(rad_10,rad_9));
    R_3 = logical(f(rad_9,rad_8));
    R_4 = logical(f(rad_8,rad_7));
    R_5 = logical(f(rad_7,rad_6));
    R_6 = logical(f(rad_6,rad_5));
    R_7 = logical(f(rad_5,rad_4));
    R_8 = logical(f(rad_4,rad_3));
    R_9 = logical(f(rad_3,rad_2));
    R_10 = logical(f(rad_2,rad_1));

    %clear rad_1 rad_2 rad_3 rad_4 rad_5 rad_6 rad_7 rad_8 rad_9 rad_10

    if leaf_length1 < leaf_length2
        leaf_coord(:,1) = sort(repmat(linspace(1, cent(2), cent(2)), 1, dims(2)));
        leaf_coord(:,2) = repmat(linspace(1, dims(2), dims(2)), 1, cent(2));
    else
        leaf_coord(:,1) = sort(repmat(linspace(cent(2), dims(1), dims(1)-cent(2)+1), 1, dims(2)));
        leaf_coord(:,2) = repmat(linspace(1, dims(2), dims(2)), 1, dims(1)-cent(2)+1); 
    end

    mask_1 = R_1 & BW_plant;
    mask_2 = R_2 & BW_plant;
    mask_3 = R_3 & BW_plant;
    mask_4 = R_4 & BW_plant;
    mask_5 = R_5 & BW_plant;
    mask_6 = R_6 & BW_plant;
    mask_7 = R_7 & BW_plant;
    mask_8 = R_8 & BW_plant;
    mask_9 = R_9 & BW_plant;
    mask_10 = R_10 & BW_plant;
    
    % select for one leaf before or after cc?
    mask_1 = bwselect(mask_1, leaf_coord(:,2), leaf_coord(:,1));
    mask_2 = bwselect(mask_2, leaf_coord(:,2), leaf_coord(:,1));
    mask_3 = bwselect(mask_3, leaf_coord(:,2), leaf_coord(:,1));
    mask_4 = bwselect(mask_4, leaf_coord(:,2), leaf_coord(:,1));
    mask_5 = bwselect(mask_5, leaf_coord(:,2), leaf_coord(:,1));
    mask_6 = bwselect(mask_6, leaf_coord(:,2), leaf_coord(:,1));
    mask_7 = bwselect(mask_7, leaf_coord(:,2), leaf_coord(:,1));
    mask_8 = bwselect(mask_8, leaf_coord(:,2), leaf_coord(:,1));
    mask_9 = bwselect(mask_9, leaf_coord(:,2), leaf_coord(:,1));
    mask_10= bwselect(mask_10,leaf_coord(:,2), leaf_coord(:,1));

    % choose one region, closest to center?
    
    % find center of each region, if >2, for the one further away from
    % plant center, make zeros in original mask
    mask_1_L = bwlabel(mask_1);
    mask_2_L = bwlabel(mask_2);
    mask_3_L = bwlabel(mask_3);
    mask_4_L = bwlabel(mask_4);
    mask_5_L = bwlabel(mask_5);
    mask_6_L = bwlabel(mask_6);
    mask_7_L = bwlabel(mask_7);
    mask_8_L = bwlabel(mask_8);
    mask_9_L = bwlabel(mask_9);
    mask_10_L = bwlabel(mask_10);
    

     %% one region donut 2
     mask = cell(donuts,1);
     
     for region =1:donuts
         maskdonut = eval(strcat("mask_", num2str(region)));
         mask_L = bwlabel(maskdonut);
         stats = regionprops(maskdonut, 'centroid');   
         Dist = [];
         
         for j= 1:length(stats)
            pt(1) = stats(j).Centroid(1);
            pt(2) = stats(j).Centroid(2);
            v1 = cent;
            v2 = extrema;
        
            stats(j).Dist = point_to_line_distance(pt, v1, v2); %calculate distance of region center point to closest point on line fit through extrema and plant center
         end
         
         %%
         MinDist = min([stats.Dist]);
         Dist = [stats.Dist];
         tf1 = Dist == MinDist;
         idx = find(tf1);
         
         mask{region} = ismember(mask_L, idx);
    
     end
  %%  
   mask_1 = mask{1};
   mask_2 = mask{2};
   mask_3 = mask{3};
   mask_4 = mask{4};
   mask_5 = mask{5};
   mask_6 = mask{6};
   mask_7 = mask{7};
   mask_8 = mask{8};
   mask_9 = mask{9};
   mask_10 = mask{10};
    
    %% make colored image
    t1 = cat(3, 166*uint8(mask_1), 206*uint8(mask_1), 227*uint8(mask_1));    % light blue
    t2 = cat(3, 31*uint8(mask_2), 120*uint8(mask_2), 180*uint8(mask_2));     % dark blue
    t3 = cat(3, 178*uint8(mask_3), 223*uint8(mask_3), 138*uint8(mask_3));    % light green
    t4 = cat(3, 51*uint8(mask_4), 160*uint8(mask_4), 44*uint8(mask_4));      % dark green
    t5 = cat(3, 251*uint8(mask_5), 154*uint8(mask_5), 153*uint8(mask_5));    % pink
    t6 = cat(3, 227*uint8(mask_6), 26*uint8(mask_6), 28*uint8(mask_6));      % red
    t7 = cat(3, 253*uint8(mask_7), 191*uint8(mask_7), 111*uint8(mask_7));    % light orange
    t8 = cat(3, 255*uint8(mask_8), 127*uint8(mask_8), 0*uint8(mask_8));      % dark orange
    t9 = cat(3, 202*uint8(mask_9), 178*uint8(mask_9), 214*uint8(mask_9));    % light purple
    t10 = cat(3, 106*uint8(mask_10), 61*uint8(mask_10), 154*uint8(mask_10)); % dark purple

    figure1 = figure('visible', 'off', 'Position', [10 10 2000 1500]);
    imshow(t1+t2+t3+t4+t5+t6+t7+t8+t9+t10,[])
    hold on
    plot(x0, y0, 'b*')
    hold off
    title('plant A applied donuts')
    %saveas(figure1, strcat(folder_data, 'Images/', date, '/' ,name,'_plant',plant), 'png')
%% Matrix of Pixel Values
    clear t1 t2 t3 t4 t5 t6 t7 t8 t9 t10
    clear mask_1_L mask_2_L mask_3_L mask_4_L mask_5_L mask_6_L mask_7_L mask_8_L mask_9_L mask_10_L
    clear stats_1 stats_2 stats_3 stats_4 stats_5 stats_6 stats_7 stats_8 stats_9 stats_10
    clear PixelValues_1 PixelValues_2 PixelValues_3 PixelValues_4 PixelValues_5 PixelValues_6 PixelValues_7 PixelValues_8 PixelValues_9 PixelValues_10
    
    mask_1 = bwconncomp(mask_1);
    mask_2 = bwconncomp(mask_2);
    mask_3 = bwconncomp(mask_3);
    mask_4 = bwconncomp(mask_4);
    mask_5 = bwconncomp(mask_5);
    mask_6 = bwconncomp(mask_6);
    mask_7 = bwconncomp(mask_7);
    mask_8 = bwconncomp(mask_8);
    mask_9 = bwconncomp(mask_9);
    mask_10 = bwconncomp(mask_10);

    for i2 = 1:size(data_norm_t, 3) % i2 = 290 bands
        PixelValues_1(:, :, i2) = regionprops(mask_1, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_2(:, :, i2) = regionprops(mask_2, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_3(:, :, i2) = regionprops(mask_3, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_4(:, :, i2) = regionprops(mask_4, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_5(:, :, i2) = regionprops(mask_5, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_6(:, :, i2) = regionprops(mask_6, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_7(:, :, i2) = regionprops(mask_7, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_8(:, :, i2) = regionprops(mask_8, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_9(:, :, i2) = regionprops(mask_9, data_norm_t(:, :, i2), 'PixelValues');
        PixelValues_10(:, :, i2) = regionprops(mask_10, data_norm_t(:, :, i2), 'PixelValues');
    end

    
    PixelValues_1 = permute(PixelValues_1, [1 3 2]);
    PixelValues_1 = struct2cell(PixelValues_1); %convert to cell: end up with 1x25x290
    PixelValues_1 = permute(PixelValues_1, [2 3 1]); %convert back to 25x290
    
    PixelValues_2 = permute(PixelValues_2, [1 3 2]);
    PixelValues_2 = struct2cell(PixelValues_2); %convert to cell: end up with 1x25x290
    PixelValues_2 = permute(PixelValues_2, [2 3 1]); %convert back to 25x290

    PixelValues_3 = permute(PixelValues_3, [1 3 2]);
    PixelValues_3 = struct2cell(PixelValues_3); %convert to cell: end up with 1x25x290
    PixelValues_3 = permute(PixelValues_3, [2 3 1]); %convert back to 25x290

    PixelValues_4 = permute(PixelValues_4, [1 3 2]);
    PixelValues_4 = struct2cell(PixelValues_4); %convert to cell: end up with 1x25x290
    PixelValues_4 = permute(PixelValues_4, [2 3 1]); %convert back to 25x290

    PixelValues_5 = permute(PixelValues_5, [1 3 2]);
    PixelValues_5 = struct2cell(PixelValues_5); %convert to cell: end up with 1x25x290
    PixelValues_5 = permute(PixelValues_5, [2 3 1]); %convert back to 25x290

    PixelValues_6 = permute(PixelValues_6, [1 3 2]);
    PixelValues_6 = struct2cell(PixelValues_6); %convert to cell: end up with 1x25x290
    PixelValues_6 = permute(PixelValues_6, [2 3 1]); %convert back to 25x290

    PixelValues_7 = permute(PixelValues_7, [1 3 2]);
    PixelValues_7 = struct2cell(PixelValues_7); %convert to cell: end up with 1x25x290
    PixelValues_7 = permute(PixelValues_7, [2 3 1]); %convert back to 25x290

    PixelValues_8 = permute(PixelValues_8, [1 3 2]);
    PixelValues_8 = struct2cell(PixelValues_8); %convert to cell: end up with 1x25x290
    PixelValues_8 = permute(PixelValues_8, [2 3 1]); %convert back to 25x290

    PixelValues_9 = permute(PixelValues_9, [1 3 2]);
    PixelValues_9 = struct2cell(PixelValues_9); %convert to cell: end up with 1x25x290
    PixelValues_9 = permute(PixelValues_9, [2 3 1]); %convert back to 25x290

    PixelValues_10 = permute(PixelValues_10, [1 3 2]);
    PixelValues_10 = struct2cell(PixelValues_10); %convert to cell: end up with 1x25x290
    PixelValues_10 = permute(PixelValues_10, [2 3 1]); %convert back to 25x290
  
    list1 = {PixelValues_1, PixelValues_2, PixelValues_3, PixelValues_4, PixelValues_5, ...
        PixelValues_6, PixelValues_7, PixelValues_8, PixelValues_9, PixelValues_10};
    
    for m=1:size(list1,2)
        for i=1:size(list1{m},1)
            for k=1:length(list1{m})
                test(:,k) = list1{m}{i,k};
            end
            %csvwrite(strcat(folder_data, 'Pixel Values/', 'data_', name, '_', string(m),'_region',string(i),'_plant',plant, '.csv'), test)
        
             DonutPixelValues{m,1} = test(:,:);
            clear test
        end
    end

%% Cell array of Mean Slopes
    
for donut = 1:10

k = eval(strcat("PixelValues_", num2str(donut)));
[num,nc] = size(k); %number of plants

plants = cell(1,num);

for i = 1:num; 
    plants{1,i} = [k{i,:}];
end

%%
SplineCurve = cell(num,1); %Getting slope line per pixel in each plant
for i = 1:num;
%A = zeros(length(plants{1,i}), 290);



y = (transpose(Wavelengths));
x = plants{1,i};

%% Breaks interpolated from data
SplineCurve{i,1} = splinefit(y,x,20);  % 20 breaks, 19 pieces
end

% Construct lines for all objects into matrix
xx = linspace(400,1000);
a = length(SplineCurve);
b = length(xx);
lines = cell(a,1);

for i = 1:a;
lines{i,:} = ppval(SplineCurve{i,1},xx);
end

lines2 = transpose(lines);
xx_cell = cell(1,num);

for p = 1:num;
xx_cell{1, p} = xx;
end

%cellfun(@plot, lines2, xx_cell);

% get slope for each plant smoothed spectra
a = length(SplineCurve);
b = length(xx);

slopes = cell(a,1);
mean_slope = zeros(a,b);

for i = 1: a;
[FY] = gradient(lines{i,1});
[FX] = gradient(xx);

    slope = zeros(length(FY(:,1)),length(FX(1,:)));

    for z = 1:length(FX(1,:));
        for l = 1:length(FY(:,1));
            slope(l,z) = FY(l,z)/FX(1,z);
        end
        mean_slope(i,z) = mean(slope(:,z));
        
    end
    DonutMeanSlopes{donut,1} = mean_slope(i,:);
    DonutAllSlopes{donut,1} = slope(:,:);
    
end

%csvwrite(strcat(folder_data, 'Mean Slopes/', 'slope_', name, '_', num2str(donut),'_region',string(i),'_plant',plant, '.csv'), mean_slope)


end

    
%% Getting the Mean
 clear PixelValues_1 PixelValues_2 PixelValues_3 PixelValues_4 PixelValues_5 PixelValues_6 PixelValues_7 PixelValues_8 PixelValues_9 PixelValues_10
 clear leaf_mask leaf_coord leaf_length1 leaf_length2 cent ex ex_br ex_tl x1 x0 y1 y0 m k i
 
 for i2 = 1:size(data_norm_t, 3) % i2 = 290 bands
        PixelValues_1(:, :, i2) = regionprops(mask_1, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_2(:, :, i2) = regionprops(mask_2, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_3(:, :, i2) = regionprops(mask_3, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_4(:, :, i2) = regionprops(mask_4, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_5(:, :, i2) = regionprops(mask_5, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_6(:, :, i2) = regionprops(mask_6, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_7(:, :, i2) = regionprops(mask_7, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_8(:, :, i2) = regionprops(mask_8, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_9(:, :, i2) = regionprops(mask_9, data_norm_t(:, :, i2), 'MeanIntensity');
        PixelValues_10(:, :, i2) = regionprops(mask_10, data_norm_t(:, :, i2), 'MeanIntensity');
    end

    %clear mask_1 mask_2 mask_3 mask_4 mask_5 mask_6 mask_7 mask_8 mask_9 mask_10
    
    PixelValues_1 = permute(PixelValues_1, [1 3 2]);
    PixelValues_1 = struct2cell(PixelValues_1); %convert to cell: end up with 1x25x290
    PixelValues_1 = permute(PixelValues_1, [2 3 1]); %convert back to 25x290
    
    PixelValues_2 = permute(PixelValues_2, [1 3 2]);
    PixelValues_2 = struct2cell(PixelValues_2); %convert to cell: end up with 1x25x290
    PixelValues_2 = permute(PixelValues_2, [2 3 1]); %convert back to 25x290

    PixelValues_3 = permute(PixelValues_3, [1 3 2]);
    PixelValues_3 = struct2cell(PixelValues_3); %convert to cell: end up with 1x25x290
    PixelValues_3 = permute(PixelValues_3, [2 3 1]); %convert back to 25x290

    PixelValues_4 = permute(PixelValues_4, [1 3 2]);
    PixelValues_4 = struct2cell(PixelValues_4); %convert to cell: end up with 1x25x290
    PixelValues_4 = permute(PixelValues_4, [2 3 1]); %convert back to 25x290

    PixelValues_5 = permute(PixelValues_5, [1 3 2]);
    PixelValues_5 = struct2cell(PixelValues_5); %convert to cell: end up with 1x25x290
    PixelValues_5 = permute(PixelValues_5, [2 3 1]); %convert back to 25x290

    PixelValues_6 = permute(PixelValues_6, [1 3 2]);
    PixelValues_6 = struct2cell(PixelValues_6); %convert to cell: end up with 1x25x290
    PixelValues_6 = permute(PixelValues_6, [2 3 1]); %convert back to 25x290

    PixelValues_7 = permute(PixelValues_7, [1 3 2]);
    PixelValues_7 = struct2cell(PixelValues_7); %convert to cell: end up with 1x25x290
    PixelValues_7 = permute(PixelValues_7, [2 3 1]); %convert back to 25x290

    PixelValues_8 = permute(PixelValues_8, [1 3 2]);
    PixelValues_8 = struct2cell(PixelValues_8); %convert to cell: end up with 1x25x290
    PixelValues_8 = permute(PixelValues_8, [2 3 1]); %convert back to 25x290

    PixelValues_9 = permute(PixelValues_9, [1 3 2]);
    PixelValues_9 = struct2cell(PixelValues_9); %convert to cell: end up with 1x25x290
    PixelValues_9 = permute(PixelValues_9, [2 3 1]); %convert back to 25x290

    PixelValues_10 = permute(PixelValues_10, [1 3 2]);
    PixelValues_10 = struct2cell(PixelValues_10); %convert to cell: end up with 1x25x290
    PixelValues_10 = permute(PixelValues_10, [2 3 1]); %convert back to 25x290
  
    list1 = {PixelValues_1, PixelValues_2, PixelValues_3, PixelValues_4, PixelValues_5, ...
        PixelValues_6, PixelValues_7, PixelValues_8, PixelValues_9, PixelValues_10};
    
    for m=1:size(list1,2)
        for i=1:size(list1{m},1)
            for k=1:length(list1{m})
                test(:,k) = list1{m}{i,k};
            end
            %csvwrite(strcat(folder_data, 'Mean Spectra/', 'mean_', name, '_', string(m),'_region',string(i),'_plant',plant, '.csv'), test)
            DonutMeanSpectra{donut,1} = test(1,:);
 
            clear test
        end
    end


%DonutMeanSpectra = list1;


  

