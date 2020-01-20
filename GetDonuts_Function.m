function [DonutMeanSpectra_1] = GetDonuts_Function (date1, donuts)

folder = strcat(date1, '/');

ext_white = 'white_';
ext_dark = 'dark_';
ext_raw = '.raw';
ext_hdr = '.hdr';

file_white1 = strcat(folder, ext_white,date1, ext_raw);
hdrfile_white1 = strcat(folder, ext_white,date1, ext_hdr);
file_dark1 = strcat(folder, ext_dark,date1, ext_raw);
hdrfile_dark1 = strcat(folder, ext_dark,date1, ext_hdr);

% Read in all files from folder
imds = imageDatastore(folder,...
'IncludeSubfolders',true,'FileExtensions','.raw');

filenum = length(imds.Files)-2;

M1 = readtable(strcat(folder, 'Center_coordinates_', date1, '.csv'), 'Delimiter', ',', 'ReadRowNames', true);

%%

for i = 1:filenum;  %%%%Change back to 1
file_i = imds.Files(i);
[filepath,name1,ext] = fileparts(file_i{1,1});

file_1 = strcat(folder,name1,ext); 
hdrfile_1 = strcat(folder, name1,ext_hdr); 

[NDVI_1, BW_plant_1, data_norm_t_1, h_1, Wavelengths_1] = readHyperspecFile_Function (file_1, hdrfile_1, file_white1, hdrfile_white1, file_dark1, hdrfile_dark1);

%% Plant A
plant_1 = 'A';


% Getting donut hyperspectral data from each of the two plants
folder_data = "Data Analysis/Donuts/"; %Have folder named: Images /  Mean Slopes / Mean Spectra / Pixel Values


[DonutPixelValues_1, DonutMeanSlopes_1, DonutAllSlopes_1, DonutMeanSpectra_1, xx_1] = ExtractDonuts_Function(NDVI_1, BW_plant_1, data_norm_t_1, h_1, Wavelengths_1, M1, name1, date1, plant_1, folder_data, donuts);


%% Plant B
plant_1 = 'B';

% Getting donut hyperspectral data from each of the two plants
folder_data = "Data Analysis/Donuts/"; %Have folder named: Images /  Mean Slopes / Mean Spectra / Pixel Values

[DonutPixelValues_1, DonutMeanSlopes_1, DonutAllSlopes_1, DonutMeanSpectra_1, xx_1] = ExtractDonuts_Function(NDVI_1, BW_plant_1, data_norm_t_1, h_1, Wavelengths_1, M1, name1, date1, plant_1, folder_data, donuts);


%% Plant C
plant_1 = 'C';

% Getting donut hyperspectral data from each of the two plants
folder_data = "Data Analysis/Donuts/"; %Have folder named: Images /  Mean Slopes / Mean Spectra / Pixel Values

[DonutPixelValues_1, DonutMeanSlopes_1, DonutAllSlopes_1, DonutMeanSpectra_1, xx_1] = ExtractDonuts_Function(NDVI_1, BW_plant_1, data_norm_t_1, h_1, Wavelengths_1, M1, name1, date1, plant_1, folder_data, donuts);

end


end

