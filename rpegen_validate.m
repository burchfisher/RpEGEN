function manifest = rpegen_validate(config)
% Validate inputs without changing the published extraction calculations.
for field = {'ROIdir', 'IMGdir'}
    if ~isfolder(config.(field{1}))
        error('RpEGEN:Folder', 'Select an existing %s folder.', field{1});
    end
end
if ~isvarname(config.GroupName)
    error('RpEGEN:Group', 'Group name must be a MATLAB identifier, such as DMSO_4dpi.');
end
validateattributes(config.Channel, {'numeric'}, {'scalar', 'integer', 'positive'});
roiFiles = dir(fullfile(config.ROIdir, '*.roi'));
imageFiles = dir(fullfile(config.IMGdir, '*.tif'));
if isempty(roiFiles) || isempty(imageFiles)
    error('RpEGEN:Empty', 'The selected folders must contain .roi and .tif files.');
end
roiNames = erase(string({roiFiles.name})', '.roi');
imageNames = erase(string({imageFiles.name})', '.tif');
missingImages = setdiff(roiNames, imageNames);
missingROIs = setdiff(imageNames, roiNames);
if ~isempty(missingImages) || ~isempty(missingROIs)
    error('RpEGEN:Mismatch', 'Missing images: %s\nMissing ROIs: %s', ...
        strjoin(missingImages, ', '), strjoin(missingROIs, ', '));
end
if ~isequal(roiNames, imageNames)
    error('RpEGEN:Order', 'Image and ROI filenames must have the same case and ordering.');
end
manifest = table(imageNames, strings(numel(imageNames),1), ...
    strings(numel(imageNames),1), zeros(numel(imageNames),1), ...
    zeros(numel(imageNames),1), 'VariableNames', ...
    {'Sample', 'ImagePath', 'ROIPath', 'Height', 'Width'});
for k = 1:height(manifest)
    manifest.ImagePath(k) = fullfile(config.IMGdir, imageFiles(k).name);
    manifest.ROIPath(k) = fullfile(config.ROIdir, roiFiles(k).name);
    info = imfinfo(manifest.ImagePath(k));
    if numel(info) < config.Channel
        error('RpEGEN:Channel', '%s has %d TIFF pages; selected channel is %d.', ...
            imageFiles(k).name, numel(info), config.Channel);
    end
    image = imread(manifest.ImagePath(k), config.Channel);
    if ~isa(image, 'uint8') || ~ismatrix(image)
        error('RpEGEN:Image', '%s must contain an 8-bit grayscale brightfield page.', imageFiles(k).name);
    end
    roi = ReadImageJROI(char(manifest.ROIPath(k)));
    if ~strcmp(roi.strName, char(imageNames(k)))
        error('RpEGEN:ROIName', 'Internal ROI name "%s" does not match "%s".', roi.strName, imageNames(k));
    end
    if ~isfield(roi, 'mnCoordinates') || size(roi.mnCoordinates,1) < 3 || ...
            ~ismember(roi.strType, {'Polygon', 'Freehand', 'Traced'})
        error('RpEGEN:Polygon', '%s must be a closed polygon, freehand, or traced ROI.', roiFiles(k).name);
    end
    xy = roi.mnCoordinates;
    if any(~isfinite(xy(:))) || polyarea(xy(:,1), xy(:,2)) == 0
        error('RpEGEN:Polygon', '%s has an invalid or zero-area polygon.', roiFiles(k).name);
    end
    manifest.Height(k) = size(image,1);
    manifest.Width(k) = size(image,2);
    if any(xy(:,1) < 0 | xy(:,1) > size(image,2) | xy(:,2) < 0 | xy(:,2) > size(image,1))
        error('RpEGEN:Bounds', '%s extends beyond its image.', roiFiles(k).name);
    end
end
end
