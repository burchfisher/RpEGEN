function file = preview_qc_centerline
% Standalone display comparison only; does not alter the GUI or saved data.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
loaded = load(fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'));
data = loaded.DMSO_4dpi; eye = 1;
center = data.CLine{eye}; polygon = data.ROIxy{eye};
polygon = [polygon;polygon(1,:)];
img = repmat(double(data.IMG{eye})/255,1,1,3);
[~,order] = sort(center.distance);
xy = [double(center.x(order)),double(center.y(order))];
values = double(center.median(order)); n = numel(values);
assert(all(isfinite(values)) && all(values>=0 & values<=255));
assert(all(vecnorm(diff(xy),2,2)<=sqrt(2)+1e-8),'Centerline has a gap.');
% Each value owns the path between neighboring midpoints, with no color blending.
midpoint = (xy(1:end-1,:)+xy(2:end,:))/2;
before = [xy(1,:);midpoint]; after = [midpoint;xy(end,:)];
colors = turbo(256); indices = min(256,max(1,floor(values/255*256)+1));
[~,q] = min(abs(double(center.angle)-90));
zoomX = double(center.x(q))+[-70 70]; zoomY = double(center.y(q))+[-70 70];
fig = figure('Visible','off','Color','w','Position',[100 100 1120 1050]);
cleanup = onCleanup(@() close(fig));
tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for panel = 1:4
    ax = nexttile(tiles); imagesc(ax,img); hold(ax,'on');
    plot(ax,polygon(:,1),polygon(:,2),'Color',[0.9 0.15 0.2],'LineWidth',1);
    proposed = mod(panel,2)==0;
    if proposed
        for k = 1:n
            plot(ax,[before(k,1),xy(k,1),after(k,1)], ...
                [before(k,2),xy(k,2),after(k,2)], ...
                'Color',colors(indices(k),:),'LineWidth',1.8);
        end
    else
        scatter(ax,center.x,center.y,10,center.median,'filled');
    end
    colormap(ax,colors); clim(ax,[0 255]);
    axis(ax,'image'); ax.YDir = 'reverse'; ax.XTick = []; ax.YTick = [];
    if panel<=2
        xlim(ax,[.5 size(img,2)+.5]); ylim(ax,[.5 size(img,1)+.5]);
        rectangle(ax,'Position',[zoomX(1),zoomY(1),diff(zoomX),diff(zoomY)], ...
            'EdgeColor',[.2 .2 .2],'LineStyle','--','LineWidth',1);
    else
        xlim(ax,zoomX); ylim(ax,zoomY);
    end
    names = {'Current: overlapping dots','Preview: colored centerline', ...
        'Current: same region enlarged','Preview: same region enlarged'};
    title(ax,names{panel},'FontSize',14);
end
bar = colorbar(ax); bar.Layout.Tile = 'east';
bar.Label.String = 'Median intensity (0-255)';
title(tiles,{'QC display comparison | DMSO 4 dpi, eye 1';data.IMGname{eye}}, ...
    'Interpreter','none','FontSize',12);
xlabel(tiles,{'Identical saved medians and color scale. Each segment keeps one point''s value; no intensity smoothing.'; ...
    'Dashed box marks the enlarged region. Display preview only; the GUI is unchanged.'},'FontSize',11);
folder = fullfile(root,'output','pdf'); if ~isfolder(folder), mkdir(folder); end
file = fullfile(folder,'qc_centerline_preview.pdf');
rpegen_export_figure(fig,file);
fprintf('Created standalone QC preview for %d unchanged centerline medians: %s\n',n,file);
end
