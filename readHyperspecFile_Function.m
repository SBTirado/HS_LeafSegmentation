% Can handle data from different binning procedures from our pipeline
% Impelments trimming function to mainting consistent trimming practices


function [NDVI, BW_plant, data_norm_t, h, Wavelengths] = readHyperspecFile_Function (file, hdrfile, file_white, hdrfile_white, file_dark, hdrfile_dark); 

%% Read file with multibandread function, and assign to variable 'data'.(filename , size [height, width, N - total number of bands], precision, offset, interleave, byteorder)
% [lines, roi_width, bands]

file = file;
hdrfile = envihdrread(hdrfile);
lines = hdrfile.lines;
roi_width = hdrfile.roi_width;
bands = hdrfile.bands;

%Wavelengths = strsplit(hdrfile.Wavelength, ',');
Wavelengths = hdrfile.Wavelength;
Wavelengths = Wavelengths(2:end-1);
Wavelengths = str2num(Wavelengths);

data = multibandread(file, [lines, roi_width, bands], 'uint16', 0, 'bil', 'ieee-le');


% % % References

%% Read white reference file with multibandread function, and assign to variable 'data_w'.
file_white = file_white;
hdrfile = envihdrread(hdrfile_white);
lines = hdrfile.lines;
roi_width = hdrfile.roi_width;
bands = hdrfile.bands;

data_w = multibandread(file_white, [lines, roi_width, bands], 'uint16', 0, 'bil', 'ieee-le');

% Read dark reference file with multibandread function, and assign to variable 'data_w'
file_dark = file_dark;
hdrfile = envihdrread(hdrfile_dark);
lines = hdrfile.lines;
roi_width = hdrfile.roi_width;
bands = hdrfile.bands;

data_d = multibandread(file_dark, [lines, roi_width, bands], 'uint16', 0, 'bil', 'ieee-le');


% Average dark reference into 2 dimensions (width and wavelengths).
data_davg = mean(data_d, 1);

% Average white reference into 2 dimensions (width and wavelengths).
data_wavg = mean(data_w, 1);

% Normalize data to white and dark references.
data_norm_num = bsxfun(@minus, data, data_davg);
data_norm_denom = data_wavg - data_davg;  
data_norm = bsxfun(@rdivide, data_norm_num, data_norm_denom);

% Transpose data
data_norm_t_whole = permute(data_norm, [2 1 3]);

%% NDVI Mask

data_norm_t = data_norm_t_whole(70:650,90:1635, :);

if length(Wavelengths)> 290
    c = data_norm_t(:,:,267);
    d = data_norm_t(:,:,354);
else
    c = data_norm_t(:, :, 134); % 134 (~678.8)
    d = data_norm_t(:, :, 177); % 177 (~779.8)
end

NDVI = ((d-c)./(d+c));
%figure, imshow(NDVI, [])

BW = roicolor(NDVI, 0.35, 1);
BW_plant = bwareafilt(BW, [1000 30000]);
h = bwconncomp(BW_plant); % separte connected pixels into plants

%imshow(BW)
%labBW = labelmatrix(h);
%figure, imshow(label2rgb(labBW), []);

 % trim noisy wavelengths off data (could probably trim more)
if length(Wavelengths)>290
        data_norm_t = data_norm_t(:,:,41:540); %41-540 or 21-270
        Wavelengths = Wavelengths(41:540);
else
        data_norm_t = data_norm_t(:,:,21:270); %41-540 or 21-270
        Wavelengths = Wavelengths(21:270);

end

% normalize all data using L2-norm method
for i=1:size(data_norm_t,1) %581
    for j=1:size(data_norm_t,2) %1546
        [s,n] = sumsqr(data_norm_t(i,j,:));
        s2 = sqrt(s);
        for k=1:size(data_norm_t,3) %190
            data_norm_t_sos(i,j,k) = data_norm_t(i,j,k)/s2;
        end
    end 
end


data_norm_t = data_norm_t_sos;

 if length(Wavelengths)> 290
    c = data_norm_t(:,:,267);
    d = data_norm_t(:,:,354);
else
    c = data_norm_t(:, :, 134); % 134 (~678.8)
    d = data_norm_t(:, :, 177); % 177 (~779.8)
end

NDVI_2 = ((d-c)./(d+c));


end


