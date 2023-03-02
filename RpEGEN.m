%% RpEGEN
% Created by G. Burch Fisher beginning on 11/1/21
% Last updated 1/10/22

% REQUIREMENTS
% 1) From MATLAB:
    % Image Processing Toolbox
    % Curve Fitting Toolbox
    % Statistics and Machine Learning Toolbox

% 2) ReadImageJROI Toolbox - provided with RpEGEN folder but also available at:
    % https://github.com/DylanMuir/ReadImageJROI
    % Dylan Muir (2021). ReadImageJROI, GitHub. Retrieved November 1, 2021.

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

%------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES
% Directory locations for ROIs, IMGs, and where you want the outputs to go.
% NOTE: Directory addresses have different slash directions on PC vs Mac

ROIdir = '/Users/name/Desktop/JoVE_dataset/DMSO_4dpi/ROIs';
IMGdir = '/Users/name/Desktop/JoVE_dataset/DMSO_4dpi/TIFs 8-bit';
Outputdir = '/Users/name/Desktop/JoVE_dataset/Output';

% Output .mat file name
filename = 'DMSO_4dpi.mat';

% Location in the tif stack of the brightfield image
bf_tif_loc = 3;     % (e.g., 3rd tiff in the tiff stack here)

% NO NEED TO ALTER ANYTHING BELOW HERE UNLESS YOU ARE MODIFYING THE CODE
%------------------------------------------------------------------------------------
%% Check that user has the Image Processing Toolbox installed
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
  % User does not have the toolbox installed.
    sprintf('Image Processing Toolbox not installed')
    return;
else
    sprintf('Image Processing Toolbox installed')    
end

clear hasIPT ans

%% Import the ROI files created in FIJI

% ROI directory
cd (ROIdir)

% Creates a cell of the roi filenames used in ReadImageJROI.m
tempnames=dir('*.roi');
ROInames={tempnames.name};
clear tempnames

% Runs the matlab function ReadImageJROI on the ROIs
% Dylan Muir (2021). ReadImageJROI (https://github.com/DylanMuir/ReadImageJROI), GitHub. Retrieved November 1, 2021.

rois = ReadImageJROI(ROInames);

% Creates a data structure with the sample name (column 1) and polygon vertices
for x = 1:size(rois,2);
    data.ROIname(x,1) = {rois{1, x}.strName}; 
    data.ROIxy(x,1) = {rois{1, x}.mnCoordinates};
end 

clear x ROInames rois

%% Import the .tif image files (currently only the brightfield channel) 

% Image directory
cd (IMGdir)

% Creates a cell of the tif filenames 
TIFnames=dir('*.tif');

% Adds the image sample name and image data to the data structure 
for x = 1:size(TIFnames,1);
    X = imread(TIFnames(x).name,bf_tif_loc);
    samp_name = split(TIFnames(x).name, '.tif');
    data.IMGname(x,1) = samp_name(1);
    data.IMG(x,1) = {X} ;
end

clear X x samp_name TIFnames

%% Check that the file names between the ROI and IMG match
if ~isequal(data.ROIname,data.IMGname)
  % The names are not the same between the images and the ROIs
    sprintf('Mismatch between image filenames and ROI filenames')
    return;
else
    sprintf('Image and ROI names match')    
end

clear ans

%% Convert ROIs into binary BW masks
% Turn xy vertices into a BW image equal in size to the IMG data
for x = 1:size(data.ROIxy,1);
    X = poly2mask(data.ROIxy{x,1}(:,1),data.ROIxy{x,1}(:,2),size(data.IMG{x,1},1),size(data.IMG{x,1},2));
    data.ROIBW(x,1) = {X} ;
end

clear x X

%% Loop through the BW masks to create the centerline and distance data by calling centderfunc.m
for z = 1:size(data.ROIBW,1);
    [dist,xl,yl] = centderfunc(data.ROIBW{z, 1}); % Call the centderfunc.m file here
    xl = transpose(xl);
    yl = transpose(yl);
    dist = transpose(dist);
    centerline = table(xl,yl,dist);
    centerline.Properties.VariableNames = {'x','y','distance'};
    data.CLine(z,1) = {centerline};
end

clear z centerline dist xl yl

sprintf('6 of 10 Steps Done')
%% Create an index map based on distance from centerline pixels
for x = 1:size(data.ROIBW,1);
    szy = size(data.ROIBW{x,1},1);
    szx = size(data.ROIBW{x,1},2);
    centBW = zeros(szy,szx);
    centIND = sub2ind([szy szx],data.CLine{x, 1}.y,data.CLine{x, 1}.x);
    centBW(centIND) = 1;
    [D,idx] = bwdist(centBW);
    maskIDX = uint32(data.ROIBW{x,1}).*idx;

    % Calculates the angle of each centerline point from the midpoint
        % between the start and end points
    axst = data.CLine{x,1}.x(data.CLine{x,1}.distance == 0);
    ayst = data.CLine{x,1}.y(data.CLine{x,1}.distance == 0);
    axe = data.CLine{x,1}.x(data.CLine{x,1}.distance == max(data.CLine{x,1}.distance));
    aye = data.CLine{x,1}.y(data.CLine{x,1}.distance == max(data.CLine{x,1}.distance));
    midptx = axst-((axst - axe)./2);
    midpty = ayst-((ayst - aye)./2);
    % I round here to keep any subtle variations in the first few pixels
        % from creating negative angles 
    tilt = round(atand((axst-axe)/(ayst-aye)),1);
    for z = 1:size(data.CLine{x,1},1);
        data.CLine{x,1}.angle(z) = round(atand((data.CLine{x,1}.x(z)-midptx)./(data.CLine{x,1}.y(z)-midpty))-tilt,1);
        if data.CLine{x,1}.angle(z) < 0;
            data.CLine{x,1}.angle(z) = 180 + data.CLine{x,1}.angle(z);
        else
            continue
        end
    end
    id0 = find(data.CLine{x,1}.angle <= 0);
    id0 = id0(id0 > size(data.CLine{x,1},1)./2);
    data.CLine{x,1}.angle(id0) = 180;
    data.CLine{x,1}.angle = abs(data.CLine{x,1}.angle - 180);
    data.CLine{x,1}.ind = centIND;
    data.CIdx{x,1} = maskIDX;
    data.CLineBW{x,1} = centBW;
end

clear a szx szy x idx D centIND centBW maskIDX axe axst aye ayst midptx midpty tilt z id0

sprintf('7 of 10 Steps Done - Next step takes a bit longer...')
%% Extract data using the CIdx and store in the CLine table
data.RAW_data = table([],[], [], 'VariableNames',{'IMG_name' 'Angle' 'Raw_Data'});

for x = 1:size(data.IMG,1);
    for z = 1:size(data.CLine{x, 1},1);
        idc = find(data.CIdx{x,1} == data.CLine{x, 1}.ind(z));
        IMG_lin = reshape(data.IMG{x,1},[],1);
        IMG_data = single(IMG_lin(idc));
        data.CLine{x, 1}.num_pts(z) = numel(idc);   % of pixels contributing to the centerline pixel
        data.CLine{x, 1}.mean(z) = mean(IMG_data);
        data.CLine{x, 1}.stdev(z) = std(IMG_data);
        data.CLine{x, 1}.stderr(z) = std(IMG_data)/sqrt(numel(idc));
        data.CLine{x, 1}.median(z) = median(IMG_data);
        data.CLine{x, 1}.p25(z) = prctile(IMG_data,25);
        data.CLine{x, 1}.p75(z) = prctile(IMG_data,75);
        data.CLine{x, 1}.min(z) = min(IMG_data);
        data.CLine{x, 1}.max(z) = max(IMG_data);        
        a = cell(numel(idc),1); 
        a(:) = {data.IMGname{x, 1}};
        b = repmat(data.CLine{x, 1}.angle(z),numel(idc), 1);
        A = table(a, b, IMG_data, 'VariableNames',{'IMG_name' 'Angle' 'Raw_Data'}); 
        data.RAW_data = vertcat(data.RAW_data, A);
    end
    % Add in the smoothed mean and median to the table in data.CLine
    data.CLine{x, 1}.mean_sm10per = smooth(data.CLine{x, 1}.angle, data.CLine{x, 1}.mean, 0.1,'moving');
    data.CLine{x, 1}.median_sm10per = smooth(data.CLine{x, 1}.angle, data.CLine{x, 1}.median, 0.1,'moving');
end

clear z x idc IMG_data IMG_lin a b A 

sprintf('8 of 10 Steps Done')
%% Degree binned data for all the raw data
% 1-degree
a = cell(180,1); b=split(filename,'.'); a(:) = {b{1,1}};
data.BIN_1_deg_all = table(a,[zeros(180,1)],[zeros(180,1)],[zeros(180,1)],[zeros(180,1)],[zeros(180,1)],[zeros(180,1)],[zeros(180,1)],[zeros(180,1)],...
    'VariableNames',{'Group' 'Angle' 'Numel' 'Mean' 'STD' 'SE' 'Median' 'CI95U' 'CI95L'});

for x=1:180;
    idx = find(data.RAW_data.Angle >= x-1 & data.RAW_data.Angle < x);
    rd = sort(data.RAW_data.Raw_Data(idx));
    data.BIN_1_deg_all.Angle(x) = x-0.5;
    data.BIN_1_deg_all.Numel(x) = numel(rd);
    data.BIN_1_deg_all.Mean(x) = mean(rd);
    data.BIN_1_deg_all.STD(x) = std(rd);
    data.BIN_1_deg_all.SE(x) = std(rd)/sqrt(numel(rd));
    data.BIN_1_deg_all.Median(x) = median(rd);
    
    % Calculate the 95% CI for the median
    j = round(0.5*numel(rd) - 1.96*sqrt(0.5*numel(rd)*(1-0.5)));
    k = round(0.5*numel(rd) + 1.96*sqrt(0.5*numel(rd)*(1-0.5)));
    data.BIN_1_deg_all.CI95U(x) = rd(k);
    data.BIN_1_deg_all.CI95L(x) = rd(j);
end
    
clear x j k rd idx b a


% 5-degree
a = cell(36,1); b=split(filename,'.'); a(:) = {b{1,1}};
data.BIN_5_deg_all = table(a,[zeros(36,1)],[zeros(36,1)],[zeros(36,1)],[zeros(36,1)],[zeros(36,1)],[zeros(36,1)],[zeros(36,1)],[zeros(36,1)],...
    'VariableNames',{'Group' 'Angle' 'Numel' 'Mean' 'STD' 'SE' 'Median' 'CI95U' 'CI95L'});

for x=5:5:180;
    z = uint8(x/5);
    idx = find(data.RAW_data.Angle >= x-5 & data.RAW_data.Angle < x);
    rd = sort(data.RAW_data.Raw_Data(idx));
    data.BIN_5_deg_all.Angle(z) = x-2.5;
    data.BIN_5_deg_all.Numel(z) = numel(rd);
    data.BIN_5_deg_all.Mean(z) = mean(rd);
    data.BIN_5_deg_all.STD(z) = std(rd);
    data.BIN_5_deg_all.SE(z) = std(rd)/sqrt(numel(rd));
    data.BIN_5_deg_all.Median(z) = median(rd);
    
    % Calculate the 95% CI for the median
    j = round(0.5*numel(rd) - 1.96*sqrt(0.5*numel(rd)*(1-0.5)));
    k = round(0.5*numel(rd) + 1.96*sqrt(0.5*numel(rd)*(1-0.5)));
    data.BIN_5_deg_all.CI95U(z) = rd(k);
    data.BIN_5_deg_all.CI95L(z) = rd(j);
end
    
clear x j k rd idx b a 


% 10-degree
a = cell(18,1); b=split(filename,'.'); a(:) = {b{1,1}};
data.BIN_10_deg_all = table(a,[zeros(18,1)],[zeros(18,1)],[zeros(18,1)],[zeros(18,1)],[zeros(18,1)],[zeros(18,1)],[zeros(18,1)],[zeros(18,1)],...
    'VariableNames',{'Group' 'Angle' 'Numel' 'Mean' 'STD' 'SE' 'Median' 'CI95U' 'CI95L'});

for x=10:10:180;
    z = uint8(x/10);
    idx = find(data.RAW_data.Angle >= x-10 & data.RAW_data.Angle < x);
    rd = sort(data.RAW_data.Raw_Data(idx));
    data.BIN_10_deg_all.Angle(z) = x-5;
    data.BIN_10_deg_all.Numel(z) = numel(rd);
    data.BIN_10_deg_all.Mean(z) = mean(rd);
    data.BIN_10_deg_all.STD(z) = std(rd);
    data.BIN_10_deg_all.SE(z) = std(rd)/sqrt(numel(rd));
    data.BIN_10_deg_all.Median(z) = median(rd);
    
    % Calculate the 95% CI for the median
    j = round(0.5*numel(rd) - 1.96*sqrt(0.5*numel(rd)*(1-0.5)));
    k = round(0.5*numel(rd) + 1.96*sqrt(0.5*numel(rd)*(1-0.5)));
    data.BIN_10_deg_all.CI95U(z) = rd(k);
    data.BIN_10_deg_all.CI95L(z) = rd(j);
end
    
clear x j k rd idx b a


sprintf('9 of 10 Steps Done - Saving .MAT file...')
%% Save .mat and .xlsx files
% Output directory
cd (Outputdir)

% Convert name to the group name from filename
nm = strsplit(filename,'.');
D.(nm{1,1}) = data;

% Save the .mat file
save(filename, '-struct', 'D');

sprintf('File Saved. Saving QC plots to output folder... ')

clear D nm

%% Save and plot QC plots to ouput folder as PDFs
% Set up the figure
figure('Position',[0 0 2100 800], 'DefaultAxesFontSize',16, 'defaultTextFontName', 'Calibri');

% Run through each IMG and create and save a 3-panel figure to PDF
for x = 1:length(data.IMG);
    %This is the polyshape of the ROI
    pgon = polyshape(data.ROIxy{x,1}(:,1), data.ROIxy{x,1}(:,2));
    z = 1;
    
    % A - ROI overlaid on the IMG
    a = subplot(1, 3, z);
    imagesc(data.IMG{x,1})
    colormap(a,'gray')
    c1 = colorbar;
    c1.Visible = 'on';
    ylabel(c1,'Pixel Intensity (0-255)','FontSize',16,'Rotation',270);
    c1.Label.Position(1) = 5;
    daspect(a,[1 1 1])
    hold on
    plot(pgon, 'EdgeColor', 'r', 'FaceColor', 'none', 'Linewidth', 1)
    hold off

    % B - ROI with angular distance calculated for each centerline pixel
    b = subplot(1, 3, z+1);

    %The polyshape CLine values have to be corrected when using plot here because image Y axes are from the top down 
    pgon = polyshape(data.ROIxy{x,1}(:,1), size(data.IMG{x, 1},1) - data.ROIxy{x,1}(:,2));
    
    plot(pgon, 'EdgeColor', 'k', 'FaceColor', 'none', 'Linewidth', 2)
    hold on
    scatter(data.CLine{x,1}.x,size(data.IMG{x, 1},1) - data.CLine{x,1}.y, 50, data.CLine{x,1}.angle, 'fill')
    colormap(b, 'turbo')
    c2 = colorbar;
    c2.Visible = 'on';
    ylabel(c2,'Angular Distance (degrees)','FontSize',16,'Rotation',270);
    c2.Label.Position(1) = 5;
    ylim([0 size(data.IMG{x, 1},1)])
    xlim([0 size(data.IMG{x, 1},2)])
    axis square
    hold off

    % C - ROI with the median intensity value calculated for each
        % centerline pixel and indices
    c = subplot(1, 3, z+2);
    plot(pgon, 'EdgeColor', 'k', 'FaceColor', 'none', 'Linewidth', 2)
    hold on
    scatter(data.CLine{x,1}.x,size(data.IMG{x, 1},1) - data.CLine{x,1}.y, 50, data.CLine{x,1}.median, 'fill')
    colormap(c, 'turbo')
    c3 = colorbar;
    c3.Visible = 'on';
    ylabel(c3,'Median Intensity Value','FontSize',16,'Rotation',270);
    c3.Label.Position(1) = 5;
    ylim([0 size(data.IMG{x, 1},1)])
    xlim([0 size(data.IMG{x, 1},2)])
    axis square
    hold off

    f = gcf;
    name = strcat(data.IMGname{x,1},'_QC_plots.pdf')
    exportgraphics(f,name)

    clf
end

close

%% Clear all data in the workspace

clear all
