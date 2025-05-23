clearvars
close all

fname = 'HST_Den_Install_Lat_1500x.mp4';
%fname = 'install_lateral_pullout.mp4';
vr = VideoReader(fname);

% initial settings for the identification of the boundaries (within this
% interval)
pile_start_y = 350;
pile_end_y = 1421;
  
pile_diameter=15.37; % mm, used to obtain the pixel 2 meters

% skipp and start end frame just raken from the pile pen main image
% tracking for same step to compare penetration
edge_limits = 10;
nf = vr.NumFrames;
fps=round(nf/59,0); %time of shortened video for fps=30
skippy = 10; % do the analysis every 10 frames


start_frame =20;
end_frame = 740;
% identifying the same positions used in the original analysis approximately 
% 5mm, 29, 81 an 68mm pen depth, these frames were checked against the original 
% main processed plug frames and are the same so can compare heave. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note by MP: I have no idea what this is, I think this was an edit by Thom %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frames2extractvert = [2*fps 7*fps 18*fps 23*fps];
%8 and 10fps,14fps in between
%frames2extractvert = [2*fps 7*fps 8*fps 10*fps 14*fps 18*fps 23*fps]
frames2extractvert2 = [2*fps:10:7*fps,7*fps:20:14*fps,14*fps:10:16*fps,14*fps:10:23*fps];
%frames2extractvert2 =[start_frame:skippy:end_frame]

frames2extractlat=[750:2:1100];
frames2extractlat2=[750, 850];

%%
% initialise the heave as zero-valued vector
max_heave = zeros(size(frames2extractvert2));
%Vert
for fidx = 1:length(frames2extractvert2)

actual_frame_idx = frames2extractvert2(fidx);

local_frame = read(vr,actual_frame_idx);


%binarize the entire domain
I = local_frame(pile_start_y:end,:,1);

% fill holes (dilate the current boundaries and fill the resulting
% boundaries)
E = edge(I);
E = imdilate(E,strel('disk',1));
E(end,:)=1;
E = imfill(E,"holes");
E(end,:)=0;

% figure(3)
% imshow(E)

%extract the biggest region (the soil)
rp = regionprops(E,'all');
[~,i]=sort([rp.Area]);
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
boundary = B{1};

ii = E;

% find the left and right limits of the pile
i1 = max([find(ii(1,:)==1,1,'first'),1])-edge_limits;
i2 = min([find(ii(1,:)==1,1,'last'),size(ii,2)])+edge_limits;
% find the left and right limits of the soil
i3 = bb(1)+edge_limits;
i4 = bb(1)+bb(3)-edge_limits;
px2mm = pile_diameter/(i2-i1-2*edge_limits);
if fidx == 1
static_i1 = i1;
static_i2 = i2;
% to plot 2.5D map of the top displacement
static_i3 = i3;
static_i4 = i4;

% the points of the mudline
vec_pos_x = linspace(static_i3,static_i4,500);
end

% i have no clue what this is? who made this edit?
pen=[1:vec_pos_x/100:100];
pile=[15.37:1.17:15.37]; % this doesnt even make sense?


%split the soil boundaries to include/exclude piles and edges
idx_boundary_pile = sum([bb(1)+boundary(:,2)>static_i1,bb(1)+boundary(:,2)<static_i2,boundary(:,1)<mean(boundary(:,1))],2)==3;
idx_boundary = (sum([~idx_boundary_pile,bb(1)+boundary(:,2)>static_i3,bb(1)+boundary(:,2)<static_i4,boundary(:,1)<mean(boundary(:,1))],2)==4);
    
idx_bd = find(idx_boundary);


figure(1)

imshow(mat2gray(I));
hold on
plot(bb(1)+boundary(idx_bd,2), bb(2)+boundary(idx_bd,1), 'b.', 'LineWidth', 1)




if fidx == 1
initial_pos = -min(-(bb(2)+boundary(idx_bd,1)));
end

try
vec_pos_y = interp1(bb(1)+boundary(idx_bd,2),-(bb(2)+boundary(idx_bd,1))+initial_pos,vec_pos_x);
vec_pos_y(and(vec_pos_x>static_i1,vec_pos_x<static_i2))=nan;
vec_pos_y=smoothdata(vec_pos_y,"movmean",40);
fig=figure(2)
fig.Position = [100 100 500 400];
hold on
fsz = 12;
set(gca, 'FontSize', fsz);
set(gcf, 'Color', 'w');
box off;
grid off;
% colormap = jet(length(frames2extractvert));
%colormap2 = jet(length(pen));

addpath('C:\Users\140007677.DUNDEE\OneDrive - University of Dundee\Desktop\Thomas Riccio - UoD\97. Conferences & Posters\YNU meet up');
% Generate a set of distinct colors using Brewer map
colormap2 = brewermap(8, 'Greys');
% Plot the lines with varying colors

scatter(vec_pos_x * px2mm - vec_pos_x(1) * px2mm, vec_pos_y * px2mm, 5, 'MarkerEdgeColor', colormap2(3,:), 'Marker', 's');
ylabel('Top surface contour $[mm]$', 'interpreter', 'latex', 'FontSize', fsz)
xlabel('Breadth $[mm]$', 'interpreter', 'latex', 'FontSize', fsz);
max_heave(fidx) =max(vec_pos_y * px2mm);
set(gca,'XAxisLocation','top')
hold on
% Assuming 'colormap2' is a colormap matrix

% Define the minimum and maximum values in 'pen' for color mapping
% min_pen = min(pen);
% max_pen = max(pen);
%  yyaxis right
% %  ylabel('Pile penetration $[mm]$', 'interpreter', 'latex', 'FontSize', fsz); 
% ylimit([100 0])
% % Define the minimum and maximum values in 'pen' for color mapping
% min_pen = min(pen);
% max_pen = max(pen);
% 
% % Scatter plot with different colors for each point based on 'colormap2'
% scatter(ones(1, length(pen)) * 57, -pen, 100, pen, 'filled');
% caxis([min_pen, max_pen]); % Set color limits based on 'pen' values
% colormap(colormap2); % Use 'colormap2' for colors
% 
% set(gca, 'YYDir', 'reverse'); % Reverse the Y-axis
% 
% % Create a color bar manually
% c = colorbar;
% c.Label.String = 'Penetration'; % Add a label to the color bar
% 
% % Draw the rectangle for penetration
% penetration_rect_x = [57 - 15.37/2, 57 + 15.37/2]; % Define the x-coordinates
% penetration_rect_y = [-max_pen, -min_pen]; % Use the min and max 'pen' values for y-coordinates
% rectangle('Position', [penetration_rect_x(1), penetration_rect_y(1), penetration_rect_x(2) - penetration_rect_x(1), penetration_rect_y(2) - penetration_rect_y(1)], 'EdgeColor', 'r');

% The rest of your code...

% The rest of your code...

% The rest of your code...


% Create a color bar manually
% c = colorbar;

% 
% % Customize the color bar limits and labels to match your data
% min_fidx = min(fidx);
% max_fidx = max(fidx);
% caxis([min_fidx, max_fidx]);
% 
% % Set color bar ticks and labels to match your colormap
% tick_values = linspace(min_fidx, max_fidx, 10); % You can adjust the number of ticks as needed
% tick_labels = arrayfun(@(x) sprintf('%.2f', x), tick_values, 'UniformOutput', false);
% c.Ticks = tick_values;
% c.TickLabels = tick_labels;
% 
% % Label the color bar
% c.Label.String = 'Your Colorbar Label';

% Adjust the color bar position if necessary
% c.Position = [x, y, width, height]; % Specify position in normalized coordinates

% Optional: If you want to change the color bar title position
% c.Label.Position = [x, y, z]; % Specify the position in data coordinates

% You can also customize other color bar properties as needed
%set(gca,'YDirection','reverse')
ylim([0 6.2])
 xlim([0 117]);
box on
hold on


end


end

for fidx = 1:length(frames2extractvert)

actual_frame_idx = frames2extractvert(fidx);

local_frame = read(vr,actual_frame_idx);


%binarize the entire domain
I = local_frame(pile_start_y:end,:,1);

E = edge(I);
E = imdilate(E,strel('disk',1));
E(end,:)=1;
E = imfill(E,"holes");
E(end,:)=0;

% figure(3)
% imshow(E)

%extract the biggest region (the soil)
rp = regionprops(E,'all');
[~,i]=sort([rp.Area]);
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
boundary = B{1};




ii = E;

% find the left and right limits of the pile
i1 = max([find(ii(1,:)==1,1,'first'),1])-edge_limits;
i2 = min([find(ii(1,:)==1,1,'last'),size(ii,2)])+edge_limits;
% find the left and right limits of the soil
i3 = bb(1)+edge_limits;
i4 = bb(1)+bb(3)-edge_limits;
px2mm = pile_diameter/(i2-i1-2*edge_limits);
if fidx == 1
static_i1 = i1;
static_i2 = i2;
% to plot 2.5D map of the top displacement
static_i3 = i3;
static_i4 = i4;

% the points of the mudline
vec_pos_x = linspace(static_i3,static_i4,500);
end

pen=[1:vec_pos_x/100:100]
pile=[15.37:1.17:15.37]


%split the soil boundaries to include/exclude piles and edges
idx_boundary_pile = sum([bb(1)+boundary(:,2)>static_i1,bb(1)+boundary(:,2)<static_i2,boundary(:,1)<mean(boundary(:,1))],2)==3;
idx_boundary = (sum([~idx_boundary_pile,bb(1)+boundary(:,2)>static_i3,bb(1)+boundary(:,2)<static_i4,boundary(:,1)<mean(boundary(:,1))],2)==4);
    
idx_bd = find(idx_boundary);


if fidx == 1
initial_pos = -min(-(bb(2)+boundary(idx_bd,1)));
end

try
vec_pos_y = interp1(bb(1)+boundary(idx_bd,2),-(bb(2)+boundary(idx_bd,1))+initial_pos,vec_pos_x);
vec_pos_y(and(vec_pos_x>static_i1,vec_pos_x<static_i2))=nan;
vec_pos_y=smoothdata(vec_pos_y,"movmean",40)
%vec_penetration(fidx)=myPileTip;
fig=figure(2)
fig.Position = [100 100 350 400];
hold on
fsz = 14;
set(gca, 'FontSize', fsz);
set(gcf, 'Color', 'w');
box off;
grid off;
% colormap = jet(length(frames2extractvert));
colormap3 = brewermap(9, 'Reds');

addpath('C:\Users\140007677.DUNDEE\OneDrive - University of Dundee\Desktop\Thomas Riccio - UoD\97. Conferences & Posters\YNU meet up');
% Generate a set of distinct colors using Brewer map
%colormap = brewermap(length(frames2extractvert2), 'Greys');
% Plot the lines with varying colors

scatter(vec_pos_x * px2mm - vec_pos_x(1) * px2mm, vec_pos_y * px2mm, 20, 'MarkerEdgeColor',colormap3(fidx+3,:), 'Marker', 's');
ylabel('Top surface contour $[mm]$', 'interpreter', 'latex', 'FontSize', fsz)
xlabel('Breadth $[mm]$', 'interpreter', 'latex', 'FontSize', fsz);

set(gca,'XAxisLocation','top')
hold on
ylim([0 6.2])
 xlim([0 117]);
box on
hold on


end


end
%%
print('heave_contour_4_position_examplenew','-dpng','-r500')
% 
% 
% figure(3)
% mat_data(isnan(mat_data))=nan;
% h = surf(mat_data*px2mm);
% h.LineStyle='none';
% view([0 0]);
% 
% % mat data is the actual interpolated data, with the columns being the 1000 points
% % defined in line 71 and rows being the change in elevation from the first
% % frame. this is all in pixels, so there is also px2mm to convert to
% % millimeters, not sure if the value is correct, you can check it with the
% % pile width vs pixel width maybe? the surface points on X is given by
% % vec_pos_x
% save('mat_upper_boundary.mat','mat_data','px2mm','vec_pos_x');

%%
%%%lateral
for fidx = 1:length(frames2extractlat)

actual_frame_idx = frames2extractlat(fidx);

local_frame = read(vr,actual_frame_idx);


%binarize the entire domain
I = local_frame(pile_start_y:end,:,1);

E = edge(I);
E = imdilate(E,strel('disk',1));
E(end,:)=1;
E = imfill(E,"holes");
E(end,:)=0;

% figure(3)
% imshow(E)

%extract the biggest region (the soil)
rp = regionprops(E,'all');
[~,i]=sort([rp.Area]);
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
boundary = B{1};




ii = E;

% find the left and right limits of the pile
i1 = max([find(ii(1,:)==1,1,'first'),1])-edge_limits;
i2 = min([find(ii(1,:)==1,1,'last'),size(ii,2)])+edge_limits;
% find the left and right limits of the soil
i3 = bb(1)+edge_limits;
i4 = bb(1)+bb(3)-edge_limits;
px2mm = pile_diameter/(i2-i1-2*edge_limits);
if fidx == 1
static_i1 = i1;
static_i2 = i2;
% to plot 2.5D map of the top displacement
static_i3 = i3;
static_i4 = i4;

% the points of the mudline
vec_pos_x = linspace(static_i3,static_i4,500);
end

pen=[1:vec_pos_x/100:100]
pile=[15.37:1.17:15.37]


%split the soil boundaries to include/exclude piles and edges
idx_boundary_pile = sum([bb(1)+boundary(:,2)>static_i1,bb(1)+boundary(:,2)<static_i2,boundary(:,1)<mean(boundary(:,1))],2)==3;
idx_boundary = (sum([~idx_boundary_pile,bb(1)+boundary(:,2)>static_i3,bb(1)+boundary(:,2)<static_i4,boundary(:,1)<mean(boundary(:,1))],2)==4);
    
idx_bd = find(idx_boundary);



figure(1)

imshow(mat2gray(I));
hold on
plot(bb(1)+boundary(idx_bd,2), bb(2)+boundary(idx_bd,1), 'b.', 'LineWidth', 1)




if fidx == 1
initial_pos = -min(-(bb(2)+boundary(idx_bd,1)));
end

try
vec_pos_y = interp1(bb(1)+boundary(idx_bd,2),-(bb(2)+boundary(idx_bd,1))+initial_pos,vec_pos_x);
vec_pos_y(and(vec_pos_x>static_i1,vec_pos_x<static_i2))=nan;
vec_pos_y=smoothdata(vec_pos_y,"movmean",40)
%vec_penetration(fidx)=myPileTip;
fig=figure(2)
fig.Position=[100 100 380 300];
hold on
fsz = 12;
set(gca, 'FontSize', fsz);
set(gcf, 'Color', 'w');
box off;
grid off;
% colormap = jet(length(frames2extractvert));
%colormap2 = jet(length(pen));

addpath('C:\Users\140007677.DUNDEE\OneDrive - University of Dundee\Desktop\Thomas Riccio - UoD\97. Conferences & Posters\YNU meet up');
% Generate a set of distinct colors using Brewer map
colormap2 = brewermap(8, 'Greys');
% Plot the lines with varying colors

scatter(vec_pos_x * px2mm - vec_pos_x(1) * px2mm, vec_pos_y * px2mm, 5, 'MarkerEdgeColor', colormap2(3,:), 'Marker', 's');
ylabel('Top surface contour $[mm]$', 'interpreter', 'latex', 'FontSize', fsz)
xlabel('Breadth $[mm]$', 'interpreter', 'latex', 'FontSize', fsz);

set(gca,'XAxisLocation','top')
hold on
% Assuming 'colormap2' is a colormap matrix

% Define the minimum and maximum values in 'pen' for color mapping
% min_pen = min(pen);
% max_pen = max(pen);
%  yyaxis right
% %  ylabel('Pile penetration $[mm]$', 'interpreter', 'latex', 'FontSize', fsz); 
% ylimit([100 0])
% % Define the minimum and maximum values in 'pen' for color mapping
% min_pen = min(pen);
% max_pen = max(pen);
% 
% % Scatter plot with different colors for each point based on 'colormap2'
% scatter(ones(1, length(pen)) * 57, -pen, 100, pen, 'filled');
% caxis([min_pen, max_pen]); % Set color limits based on 'pen' values
% colormap(colormap2); % Use 'colormap2' for colors
% 
% set(gca, 'YYDir', 'reverse'); % Reverse the Y-axis
% 
% % Create a color bar manually
% c = colorbar;
% c.Label.String = 'Penetration'; % Add a label to the color bar
% 
% % Draw the rectangle for penetration
% penetration_rect_x = [57 - 15.37/2, 57 + 15.37/2]; % Define the x-coordinates
% penetration_rect_y = [-max_pen, -min_pen]; % Use the min and max 'pen' values for y-coordinates
% rectangle('Position', [penetration_rect_x(1), penetration_rect_y(1), penetration_rect_x(2) - penetration_rect_x(1), penetration_rect_y(2) - penetration_rect_y(1)], 'EdgeColor', 'r');

% The rest of your code...

% The rest of your code...

% The rest of your code...


% Create a color bar manually
% c = colorbar;

% 
% % Customize the color bar limits and labels to match your data
% min_fidx = min(fidx);
% max_fidx = max(fidx);
% caxis([min_fidx, max_fidx]);
% 
% % Set color bar ticks and labels to match your colormap
% tick_values = linspace(min_fidx, max_fidx, 10); % You can adjust the number of ticks as needed
% tick_labels = arrayfun(@(x) sprintf('%.2f', x), tick_values, 'UniformOutput', false);
% c.Ticks = tick_values;
% c.TickLabels = tick_labels;
% 
% % Label the color bar
% c.Label.String = 'Your Colorbar Label';

% Adjust the color bar position if necessary
% c.Position = [x, y, width, height]; % Specify position in normalized coordinates

% Optional: If you want to change the color bar title position
% c.Label.Position = [x, y, z]; % Specify the position in data coordinates

% You can also customize other color bar properties as needed
%set(gca,'YDirection','reverse')
%ylim([0 1.6])
 xlim([0 117]);
box on
hold on


end


end

for fidx = 1:length(frames2extractlat2)

actual_frame_idx = frames2extractlat2(fidx);

local_frame = read(vr,actual_frame_idx);


%binarize the entire domain
I = local_frame(pile_start_y:end,:,1);

E = edge(I);
E = imdilate(E,strel('disk',1));
E(end,:)=1;
E = imfill(E,"holes");
E(end,:)=0;

% figure(3)
% imshow(E)

%extract the biggest region (the soil)
rp = regionprops(E,'all');
[~,i]=sort([rp.Area]);
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
boundary = B{1};




ii = E;

% find the left and right limits of the pile
i1 = max([find(ii(1,:)==1,1,'first'),1])-edge_limits;
i2 = min([find(ii(1,:)==1,1,'last'),size(ii,2)])+edge_limits;
% find the left and right limits of the soil
i3 = bb(1)+edge_limits;
i4 = bb(1)+bb(3)-edge_limits;
px2mm = pile_diameter/(i2-i1-2*edge_limits);
if fidx == 1
static_i1 = i1;
static_i2 = i2;
% to plot 2.5D map of the top displacement
static_i3 = i3;
static_i4 = i4;

% the points of the mudline
vec_pos_x = linspace(static_i3,static_i4,500);
end

pen=[1:vec_pos_x/100:100]
pile=[15.37:1.17:15.37]


%split the soil boundaries to include/exclude piles and edges
idx_boundary_pile = sum([bb(1)+boundary(:,2)>static_i1,bb(1)+boundary(:,2)<static_i2,boundary(:,1)<mean(boundary(:,1))],2)==3;
idx_boundary = (sum([~idx_boundary_pile,bb(1)+boundary(:,2)>static_i3,bb(1)+boundary(:,2)<static_i4,boundary(:,1)<mean(boundary(:,1))],2)==4);
    
idx_bd = find(idx_boundary);


if fidx == 1
initial_pos = -min(-(bb(2)+boundary(idx_bd,1)));
end

try
vec_pos_y = interp1(bb(1)+boundary(idx_bd,2),-(bb(2)+boundary(idx_bd,1))+initial_pos,vec_pos_x);
vec_pos_y(and(vec_pos_x>static_i1,vec_pos_x<static_i2))=nan;
vec_pos_y=smoothdata(vec_pos_y,"movmean",40)
%vec_penetration(fidx)=myPileTip;
fig=figure(2)
fig.Position = [100 100 380 300];
hold on
fsz = 12;
set(gca, 'FontSize', fsz);
set(gcf, 'Color', 'w');
box off;
grid off;
 colormap = jet(length(frames2extractvert));
colormap = brewermap(9, 'Reds');

addpath('C:\Users\140007677.DUNDEE\OneDrive - University of Dundee\Desktop\Thomas Riccio - UoD\97. Conferences & Posters\YNU meet up');
% Generate a set of distinct colors using Brewer map
%colormap = brewermap(length(frames2extractvert2), 'Greys');
% Plot the lines with varying colors

scatter(vec_pos_x * px2mm - vec_pos_x(1) * px2mm, vec_pos_y * px2mm, 20, 'MarkerEdgeColor',colormap(fidx+6,:), 'Marker', 's');
ylabel('Top surface contour $[mm]$', 'interpreter', 'latex', 'FontSize', fsz)
xlabel('Breadth $[mm]$', 'interpreter', 'latex', 'FontSize', fsz);

set(gca,'XAxisLocation','top')
hold on
%ylim([0 1.8])
 xlim([0 117]);
box off
hold on


end


end
print('Lateral heave_Example1 Figure','-dpng','-r500')



% mat data is the actual interpolated data, with the columns being the 1000 points
% defined in line 71 and rows being the change in elevation from the first
% frame. this is all in pixels, so there is also px2mm to convert to
% millimeters, not sure if the value is correct, you can check it with the
% pile width vs pixel width maybe? the surface points on X is given by
% vec_pos_x
