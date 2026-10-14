function summary = rpegen_heatmap(raw, angleWidth, intensityWidth, normalization)
% Use identical histogram edges for every group in a comparison.
validateattributes(angleWidth, {'numeric'}, {'scalar','positive','integer'});
validateattributes(intensityWidth, {'numeric'}, {'scalar','positive','integer'});
if mod(180,angleWidth) ~= 0 || mod(255,intensityWidth) ~= 0
    error('RpEGEN:Bins','Bin widths must divide 180 degrees and 255 intensity units.');
end
if ~istable(raw) || ~all(ismember({'Angle','Raw_Data'},raw.Properties.VariableNames)) || isempty(raw)
    error('RpEGEN:RawData','Results must contain a nonempty RAW_data table with Angle and Raw_Data.');
end
angle = double(raw.Angle); intensity = double(raw.Raw_Data);
if ~iscolumn(angle) || ~iscolumn(intensity) || ...
        any(~isfinite(angle) | angle < 0 | angle > 180) || ...
        any(~isfinite(intensity) | intensity < 0 | intensity > 255)
    error('RpEGEN:RawData','Angles must be between 0 and 180 and intensities between 0 and 255, with no missing values.');
end
summary.AngleEdges = 0:angleWidth:180;
summary.IntensityEdges = 0:intensityWidth:255;
summary.Counts = histcounts2(angle,intensity,summary.AngleEdges,summary.IntensityEdges)';
summary.Values = summary.Counts;
switch normalization
    case 'Pixel count'
        summary.ColorLabel = 'Pixel count';
    case 'Percent per angle'
        summary.Values = 100*summary.Counts ./ max(sum(summary.Counts,1),1);
        summary.ColorLabel = 'Pixels per angular bin (%)';
    otherwise
        error('RpEGEN:Normalization','Unknown heatmap normalization.');
end
summary.PixelCount = height(raw);
end
