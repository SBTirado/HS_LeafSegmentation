Extract pixel values, mean spectra, mean slope and pixel slope values for 10 leaf segments across a maize plant's longest leaf from top-down hyperspectral images. 


Usage:
[DonutPixelValues, DonutMeanSlopes, DonutAllSlopes, DonutMeanSpectra, xx] = ExtractDonuts_Function(NDVI, BW_plant, data_norm_t, h, Wavelengths, M, name, date, plant, folder_data);

Inputs:

NDVI, BW_plant, data_norm_t, h, Wavelengths = Output from from readHyperspecFile3.m function

M = matrix with center coordinates for plants A, B and C of each image file processed

name = orig file name to be used for naming new files (date, genotype, stress, etc.)

date = date of imaging in mmddyyyy format as named in folder and images

plant = A, B or C

folder_data = folder where all the donut output will be stored. The folders Imaging/ Mean Slopes/ Mean SPectra and Pixel Values must be contained in the folder specified by this path 


Outputs:

DonutPixelValues = cell matrix with all pixel values for each donut for given plant. 
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}
    
DonutMeanSlopes = Average slope of spectra for all pixels in given donut. 
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}
    
DonutAllSlopes = Slope of spectra for all pixels in given donut.  
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}
    
DonutMeanSpectra =  Average spectra for all pixels in given donut.  
    % Ex. donut # output will be in DonutAllSlopes{donut#,1}

xx = wavelength breaks for DonutMeanSlopes
 
