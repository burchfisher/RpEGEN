%% RpEGEN_PermPlot
% Created by G. Burch Fisher beginning on 11/21/21
% Last updated 1/10/22

% This is a script to run the permutation simulation that compares the median 
% values across the binned raw data of two datasets to prodiuce p-values. P-values 
% less than or equal to 0.05 indicates a statistically significant difference 
% between the medians of the two groups (e.g., IWR1_4dpi and DMSO_4dpi).
% In addition, this script contains the code to replicate the figures published 
% in the JoVE article from the provided dataset (see Leach et al, JoVE).

% NOTES
% The plots require the GRAMM toolbox in your Matlab path. It can be
    % downloaded for free at:
    % https://www.mathworks.com/matlabcentral/fileexchange/54465-gramm-complete-data-visualization-toolbox-ggplot2-r-like

% The code in this M file will require interaction and input to run properly. It 
    % has been created for our needs but can be easily modified by someone
    % with Matlab experience.

%% SECTION 1 - USER DEFINED VARIABLES
% Directory where the .mat files are located from RpEGEN.m
cd ('/Users/burch/Desktop/il34 RpEGEN/Output');

% Load all the .mat files and give them names (since data is the workspace structure for them all)
load('DMSO_4dpi.mat');  
load('IWR1_4dpi.mat');  
load('DMSO_9dpf.mat');  
load('IWR1_9dpf.mat');  

%% SECTION 2 - Run permutation simulation to get the p-value comparing the medians of each bin
%--------------------------------------------------------------------------
% ENTER DATA HERE
% Enter the name of the TWO DATASETS for comparison but keep the .RAW_data after the name
% Example: [Your_Group_Name].RAW_data 

% Dataset A and B for comparison
data_A = DMSO_4dpi.RAW_data;    
data_B = IWR1_4dpi.RAW_data;

%--------------------------------------------------------------------------
bin_sz = 1;                     % Default is 1 degree bins but you can change the number of degrees you want each bin to encompass here
pval = zeros(180/bin_sz,1);     % Intiating the variable that will store all of the p-values
reps = 2000;                   % # of permutations to run for each bin (bigger # = longer processing time but more statistically robust) 

for x=bin_sz:bin_sz:180;
    z=x/bin_sz;
    ida = find(data_A.Angle >= x-bin_sz & data_A.Angle < x);
    idb = find(data_B.Angle >= x-bin_sz & data_B.Angle < x);
    pval(z,1) = permsim(data_A.Raw_Data(ida),data_B.Raw_Data(idb),reps);
    perm_status = sprintf('%.1f%% complete',100*z/(180/bin_sz))
end

% This variable is used for plotting with pval in Fig 2B below
pval_x=bin_sz:bin_sz:180; pval_x=pval_x-(bin_sz/2);

clear x z ida idb reps data_A data_B perm_status

%% HEATMAP FIGURE - PLOT USING GRAMM (Figure 9 A-D from JoVE article) 
% Heatmap figure of all the data for each group (here a 1 x 4 plot)

% Need to have GRAMM toolbox in your Matlab path. You can download for free at 
% https://www.mathworks.com/matlabcentral/fileexchange/54465-gramm-complete-data-visualization-toolbox-ggplot2-r-like

% ENTER DATA HERE ---------------------------------------------------------
% Example: Replace groups in data and nms variables below. Can be any number 
    % of groups but more than 4 may start to compress each plot horizontally.  

data = [DMSO_9dpf; IWR1_9dpf; DMSO_4dpi; IWR1_4dpi];

%--------------------------------------------------------------------------
clear g

% To extract the group names into nms
for x = 1:length(data);
    nms{x,1} = data(x).BIN_1_deg_all.Group(1);
end

% Colormap Used
cmap = 'turbo';    

% # of columns and rows for the 2d histogram bins
bin_num = [36 51];  

% THE PLOTS
figure('Position',[100 100 500*numel(nms) 500])   %[Left Bottom Width Height]

% Create the GRAMMS (e.g., g(1,1) = ...)
for z = 1:numel(data);
    g(1,z)=gramm('x',data(z).RAW_data.Angle,'y',data(z).RAW_data.Raw_Data);
    g(1,z).stat_bin2d('nbins',bin_num,'geom','image');
    g(1,z).set_continuous_color('colormap',cmap);
    g(1,z).set_title(nms{z});
    g(1,z).no_legend();
end

% General settings for the plot
g.set_text_options('font','Calibri',...
    'base_size',20,...
    'label_scaling',1.2,...
    'legend_scaling',1.5,...
    'legend_title_scaling',1.5,...
    'facet_scaling',1,...
    'title_scaling',1.5);

g.axe_property('ylim',[0 255], 'xlim', [0 180], 'TickDir','out','XGrid','on',...
    'Ygrid','on','GridColor',[0.5 0.5 0.5], 'GridAlpha', 0.5, 'Linewidth', 1.5);

g.set_names('x','Angular Distance (0=dorsal, 180=ventral)','y','Pixel Intensity (0-255)');

% Draw the figure
g.draw();

% Export figure to a PDF
f = gcf;
exportgraphics(f,'Heatmaps_fig.pdf')

clear x1 x2 x3 x4 y1 y2 y3 y4 f g cmap bin_num ans data nms z

%% GROUP RESULTS AND P-VALUE FIGURE - PLOT USING GRAMM (Figure 10 A,B from JoVE article)
% Plot of A) group medians with 95% CI and B) p-values from permutation simulations (2 x 1 plot)

% Need to have GRAMM toolbox in your Matlab path. You can download for free at 
% https://www.mathworks.com/matlabcentral/fileexchange/54465-gramm-complete-data-visualization-toolbox-ggplot2-r-like

% This section reorganizes some of the data into a long table to be used in
    % Part A of the following plot and splines the data
%--------------------------------------------------------------------------
% ENTER DATA HERE
% Here you should enter all the groups you want to plot with the .BIN_5_deg_all on the end of each
% This can be more or less than the 4 examples belowand you can also change
    % to .BIN_1_deg_all or .BIN_10_deg_all to compare. 
    
data = [DMSO_4dpi.BIN_5_deg_all; IWR1_4dpi.BIN_5_deg_all; DMSO_9dpf.BIN_5_deg_all; IWR1_9dpf.BIN_5_deg_all];

% RGB triplets for coloring the 95% CI envelopes
polcol = [1 0 0; 0 1 0; 0 0.6 1; 1 0 1]; 
%--------------------------------------------------------------------------

grps = unique(data.Group);
spdata = table([],[],[],[],[],'VariableNames',{'Group' 'Angle' 'Median' 'CI95U' 'CI95L'});

% Data is splined from 2.5:1:177.5
for x = 1:length(grps);
    a = cell(176,1); a(:) = {grps{x,1}};
    b = transpose(2.5:1:177.5);
    idx = find(strcmp(grps{x,1},data.Group));
    med = spline(data.Angle(idx),data.Median(idx), b);
    ciu = spline(data.Angle(idx),data.CI95U(idx), b);
    cil = spline(data.Angle(idx),data.CI95L(idx), b);
    df = table(a,b,med,ciu,cil,'VariableNames',{'Group' 'Angle' 'Median' 'CI95U' 'CI95L'});
    spdata = [spdata; df];
end

clear a b idx med ciu cil df x ans 

%--------------------------------------------------------------------------
% FIGURE PART A - 5-degree binned median with 95% CI envelopes
%--------------------------------------------------------------------------
clear g

% Set figure location
figure('Position',[100 100 850 1600])

% Set the data source and tpe of plot
g(1,1)=gramm('x',spdata.Angle,'y',spdata.Median,'color',spdata.Group);
g(1,1).geom_line();

% Polygon for the central section that has not regenerated 
g(1,1).geom_polygon('x',{[41 142 142 41]},'y',{[0 0 255 255]},'color',[0 0.6 1],'alpha',0.05);

% Polygons for the median 95% confidence envelope for each group in spdata
for x=1:numel(grps);
    idx = find(strcmp(spdata.Group,grps{x})==1);   
    a = {[spdata.Angle(idx); flip(spdata.Angle(idx))]}; 
    b = {[spdata.CI95U(idx); flip(spdata.CI95L(idx))]};
    g(1,1).geom_polygon('x',a,'y',b,'color',polcol(x,:),'alpha',0.3);
end

% To move the legend into the figure
g(1,1).set_layout_options('legend_position',[0.8 0.885 0.1 0.1]) %We detach the legend from the plot and move it to the top right

% General settings for the plot
g(1,1).axe_property('ylim',[0 255], 'xlim', [0 180], 'TickDir','out','XGrid','on',...
    'Ygrid','on','GridColor',[0.5 0.5 0.5], 'GridAlpha', 0.1, 'Linewidth', 2);

g(1,1).set_text_options('font','Calibri',...
    'base_size',20,...
    'label_scaling',1.3,...
    'legend_scaling',1,...
    'legend_title_scaling',1,...
    'facet_scaling',1,...
    'title_scaling',1.5);

% X and Y label names
g(1,1).set_names('x','Angular Distance (0=dorsal, 180=ventral)','y','Median Pixel Intensity (0-255)');

clear z

%--------------------------------------------------------------------------
%FIGURE PART B - 1-degree binned p-values from permutation simulation
%--------------------------------------------------------------------------
% Set the data source and tpe of plot
g(2,1)=gramm('x',pval_x,'y',pval);
g(2,1).geom_bar('dodge',0,'width',bin_sz,'FaceColor',[0.7 0 1],'EdgeColor','k', 'Linewidth', 1);

% Polygon for the central section that has not regenerated 
g(2,1).geom_polygon('x',{[41 142 142 41]},'y',{[0 0 1 1]},'color',[0 0.6 1],'alpha',0.05);

% Horizontal line showing the 0.05 p-value equal to the 95% CI
g(2,1).geom_hline('yintercept',0.05,'style','k--','Linewidth', 3)

% General settings for the plot
g(2,1).axe_property('ylim',[0 1], 'xlim', [0 180], 'TickDir','out','XGrid','on',...
    'Ygrid','on','GridColor',[0.5 0.5 0.5], 'GridAlpha', 0.1, 'Linewidth', 2);

g(2,1).set_text_options('font','Calibri',...
    'base_size',20,...
    'label_scaling',1.3,...
    'legend_scaling',1,...
    'legend_title_scaling',1,...
    'facet_scaling',1,...
    'title_scaling',1.5);

% X and Y label names
g(2,1).set_names('x','Angular Distance (0=dorsal, 180=ventral)','y','p-value');

% Draw the figure
g.draw()

% Place 95% CI text label 
text(43,0.075,{'95% CI'},'Parent',g(2,1).facet_axes_handles(1),'FontName','Calibri','FontSize', 25);

% Export figure to a PDF
f = gcf;
exportgraphics(f,'Median_pvalue_5_deg_bin_fig.pdf')

clear a b g f ans idx grps x polcol
