function rpegen_draw_qc(imageAxes,profileAxes,data,k)
% Keep on-screen and exported QC coordinates, overlays, and scales identical.
img = data.IMG{k}; xy = data.ROIxy{k}; center = data.CLine{k};
xy = [xy;xy(1,:)]; rgb = repmat(double(img)/255,1,1,3);
names = {'Image and ROI','Centerline Angle','Median Intensity'};
tags = {'QCImageAxes','QCAngleAxes','QCMedianAxes'};
for j = 1:3
    ax = imageAxes(j); cla(ax); imagesc(ax,rgb); hold(ax,'on');
    plot(ax,xy(:,1),xy(:,2),'Color',[0.9 0.15 0.2],'LineWidth',1,'Tag','QCPolygon');
    axis(ax,'image'); ax.YDir = 'reverse'; ax.XTick = []; ax.YTick = [];
    xlim(ax,[.5 size(img,2)+.5]); ylim(ax,[.5 size(img,1)+.5]);
    title(ax,names{j}); ax.Tag = tags{j};
end
plot(imageAxes(1),center.x,center.y,'Color',[0.1 0.8 0.8],'LineWidth',1,'Tag','QCCenterline');
scatter(imageAxes(2),center.x,center.y,10,center.angle,'filled','Tag','QCAnglePoints');
scatter(imageAxes(3),center.x,center.y,10,center.median,'filled','Tag','QCMedianPoints');
% All three panels reserve the same space for their color scale.
for j = 1:3
    ax = imageAxes(j);
    if j==1, colormap(ax,gray(256)); else, colormap(ax,turbo(256)); end
    upper = 255; if j==2, upper = 180; end
    clim(ax,[0 upper]); cb = colorbar(ax); cb.Ticks = [0 round(upper/2) upper];
    cb.TickLabels = {'0',sprintf('%g',round(upper/2)),sprintf('%g',upper)};
    hold(ax,'off');
end
[angles,order] = sort(center.angle);
plot(profileAxes,angles,center.median(order),'Color',[0.1 0.55 0.45],'Tag','QCIntensityProfile');
xlim(profileAxes,[0 180]); ylim(profileAxes,[0 255]); grid(profileAxes,'on');
xlabel(profileAxes,'Angular distance (degrees)'); ylabel(profileAxes,'Median intensity');
title(profileAxes,data.IMGname{k},'Interpreter','none'); profileAxes.Tag = 'QCProfileAxes';
end
