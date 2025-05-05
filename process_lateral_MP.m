clearvars
clear all
close all
addpath('../../..')

% name of the file to load/process
fname = 'install_lateral_pullout.mp4';

% filenames for output files
filename_out = strrep(fname,'.mp4','-lateral_postprocessed.gif');
filename_out2 = strrep(fname,'.mp4','-sides-lateral_postprocessed.gif');
filename_vid = strrep(fname,'.mp4','-lateral_postprocessed.mp4');



mkdir("frames_lateral");
% pile diameter to convert from pixel coordinates to real coordinates
pile_diameter = 0.01537; % 1 cm? idk

% if you want to see what is going on, set this to false
% otherwise, it will plot for each step, making the process slower
invisible = false;


% how much to cut from the perceived domain boundaries left and right
% this is because boundary lines can get blurry depending on xray intensity
edge_limits_domain = 30; %pixels

% how much to add to the perceived pile upper boundaries left and right
% this is used to cut a sub-domain portion around the pile to aid the
% identification of the boundaries
edge_limits_pile = 40; %pixels

% number of segments in a moving window average to smooth derivatives 
% (local gradients can be very noisy)
num_segments_smoothing = 50;

% load the video
vr = VideoReader(fname);


% data imports to plot force-displacement data curves on top, synced using
% the time / clock
% data = importdata("HST_dense_processed_130723.txt")
% FY = data(:,3);
% DY = data(:,2);
% data = importdata("HST_dense_Lat_processed_130723.txt")
% DLAT = data(:,1);
% FLAT = data(:,2);

% position of the pile tip at the end of penetration from the mudline
pile_tip_distance_at_eod = 888; %pixels



% load the first frame of the input video and select the upper portion of
% the region of interest (i.e. to exclude the frame bits seen in the video)
frame_ref = read(vr,1);

imshow(frame_ref);
title('Select upper limit of the sand (go slightly above)')
[~,start_pos_y] = ginput(1);

% calculate the framerate for the MP4 data. it says it has a framerate of
% 29.44 but it's not true
n_frames = vr.NumFrames;
actual_framerate = 10;%vr.NumFrames/29.440;

% do the analysis every X frames
skippy = 100;


%%
% from watching the input video, i only care about what happens between 31
% and 44 seconds
start_frame = 31*60*vr.FrameRate;
end_frame = 44*60*vr.FrameRate;

duration = 10;

% creating the output video
outputVideo = VideoWriter(filename_vid,'MPEG-4');
outputVideo.FrameRate = 60;
open(outputVideo)

close all
% start_frame = 14700
% end_frame = 16500;
% you can either select all frames from beginning to end, or create a
% vector of selected frames to loop through, same thing
vec_selected_frames = [14700:skippy:16500]; %[16500:17700];
start_frame = vec_selected_frames(1);
for fidx = vec_selected_frames%start_frame:skippy:end_frame

local_frame = read(vr,fidx);

disp(['Frame: ' num2str(fidx) '/' num2str(end_frame)])
%binarize the entire domain
warning off
I = local_frame(start_pos_y:end,:,1);
warning on
ii = I<255; % unsigned 8 bit integer -> 255 is the maximum value = white
I2 = double(I).*ii; % copy the image and keep everything except full white

%imshow(I)

%% extract the biggest region (the soil)
rp = regionprops(ii,'all');
[~,i]=sort([rp.Area]); %sort by area
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
% boundary is the boundary of both pile and domain
boundary = B{1}; 

% find the left and right limits of the pile
% the first line of the image is used -> select the edges based on the
% maximum value -> just outside the pile
i1 = max([find(ii(1,:)==1,1,'first')-edge_limits_pile/2,1]);
i2 = min([find(ii(1,:)==1,1,'last')+edge_limits_pile/2,size(ii,2)]);

% find the left and right limits of the soil by translating the bounding
% box of the soil 
i3 = bb(1)+edge_limits_domain;
i4 = bb(1)+bb(3)-edge_limits_domain;

% use the domain edges identified in the first frame for the rest of the
% analysis
if fidx == start_frame
% to plot 2.5D map of the top displacement
    static_i3 = i3;
    static_i4 = i4;
    % set pixel2meter as the ratio of pile diameter to what are identified
    % as the pile boundaries outside the domain
    % 1137 and 1016 correspond to the pile width in this specific case
    warning(['Pile width in pixel selected as: ' num2str(1137-1016) ' in line 139']);
    px2m = pile_diameter/(1137-1016);%(i2-i1)%-2*edge_limits_pile);
    ref_I = frame_ref(start_pos_y:end,:,1);

    % in the first step, get the max depth to search as the maximum tip
    % position measured
    max_depth_search = pile_tip_distance_at_eod; 
else
    % otherwise, check within 10% of the domain below the previous tip
    max_depth_search = myPileTip + length(frame_ref)/10;
end

% this just counts the number of times the boundary is crossed to identify them 
% three for the pile
idx_boundary_pile = sum([bb(1)+boundary(:,2)>i1,bb(1)+boundary(:,2)<i2,boundary(:,1)<mean(boundary(:,1))],2)==3;
% four for the domain
idx_boundary_domain = (sum([~idx_boundary_pile,bb(1)+boundary(:,2)>i3,bb(1)+boundary(:,2)<i4,boundary(:,1)<mean(boundary(:,1))],2)==4);

% pile center = the mean of the edges at the top of the domain
pile_center = round(mean([i2,i1]));

% find the position of the boundary based on the first y position in the frame
y_idx_mudline = boundary(find(idx_boundary_pile,1),1);
% translate the pile tip by the mudline offset
myPileTip = y_idx_mudline+pile_tip_distance_at_eod;

%identify the pile contour using edge analysis in a subspace beween the
%mudline, the tip and a slightly larger horizontal width than the pile
%edges at the mudline
rp = [];
% create the subdomain used to find the pile edges
I3 = mat2gray(I(y_idx_mudline:myPileTip,i1-edge_limits_pile:i2+edge_limits_pile));
%% carry out a Canny edge algorithm filtering on a smoothed image
I31 = I3(edge_limits_pile:end,1:size(I3,2)/2);
minThrs = 0.15;
maxThrs = 0.3;
% Gaussian smoothing is used to remove noise
EE = edge(imgaussfilt(mat2gray(I31)),'Canny',[minThrs maxThrs]);


rp_left.MajorAxisLength = 1;
rp_left.Orientation= 90;
rp_left.BoundingBox(4)=1;
join_size = 1;

% keep increasing the structural element size (strel) until a suitable
% sized domain is identified
% https://uk.mathworks.com/help/images/ref/strel.html
while rp_left.BoundingBox(4) * 1/sind(rp_left.Orientation)< size(I31,1)
SE = strel('square', join_size);
bw_left_edge = imdilate(EE,SE);

rp_left = regionprops(bw_left_edge,'all');
[~,i] = sort([rp_left.MajorAxisLength]);
rp_left = rp_left(i);
rp_left = rp_left(end);

join_size = join_size+1;
end
% 

%imshow(bw_left_edge)
pixel_list_left = rp_left.PixelList+[1 edge_limits_pile];
% just sorting the edges in the correct order
[~,i]=sort(pixel_list_left(:,1));%,'descend');
pixel_list_left=pixel_list_left(i,:);
[~,i]=unique(pixel_list_left(:,2));
pixel_list_left = pixel_list_left(i,:);


% repeat to identify the right edges
I32 = I3(edge_limits_pile:end,size(I3,2)/2:end);
EE = edge(imgaussfilt(mat2gray(I32)),'Canny',[minThrs maxThrs]);


rp_right.MajorAxisLength = 1;
join_size = 1;
rp_right.Orientation= 90;
rp_right.BoundingBox(4)=1;
while rp_right.BoundingBox(4) * abs(1/sind(rp_right.Orientation))< size(I32,1)
SE = strel('square', join_size);
bw_right_edge = imdilate(EE,SE);
rp_right = regionprops(bw_right_edge,'all');
[~,i] = sort([rp_right.MajorAxisLength]);
rp_right = rp_right(i);
rp_right = rp_right(end);
join_size = join_size+1;
end



pixel_list_right =rp_right.PixelList+[size(I3,2)/2-1 edge_limits_pile];

% remove doubles on the X direction
[~,i]=sort(pixel_list_right(:,1));
pixel_list_right=pixel_list_right(i,:);
[~,i]=unique(pixel_list_right(:,2));
pixel_list_right = pixel_list_right(i,:);

tt = (fidx-start_frame)/(end_frame-start_frame);
color = [tt,0,1-tt];
 figure(112)
 clf
 imshow(I3)
 hold on
 plot(pixel_list_left(:,1),pixel_list_left(:,2),'-','Color','y','LineWidth',0.5)
 plot(pixel_list_right(:,1),pixel_list_right(:,2),'-','Color','y','LineWidth',0.5)
 

 tmp_idx = find(fidx==vec_selected_frames,1);
tmp_vec_x = [183;182;173];
tmp_vec_x = linspace(183,183,length(vec_selected_frames));
x_pos_mudline = tmp_vec_x(tmp_idx);

% create a list of pixels for the edges
pixel_list_right = [x_pos_mudline 1;pixel_list_right];
start_x_poly = min(x_pos_mudline(:,1));

pixel_list_right2 = pixel_list_right;
pixel_list_right = [pixel_list_right(:,1) -pixel_list_right(:,2)];

% convert the list of pixels global scale (meters)
m_right_list = pixel_list_right * px2m;
m_left_list=pixel_list_left*px2m;

% translate to 0
centralized_right_list = m_right_list-min(m_right_list);
centralized_left_list=m_left_list-min(m_left_list);

initial_position = [i1,-y_idx_mudline]*px2m+min(m_right_list);

% fit third degree polynomials for the successive differentiation
poly_right = polyfit(centralized_right_list(:,2),centralized_right_list(:,1),3);
poly_left=polyfit(centralized_left_list(:,2),centralized_left_list(:,1),3);

% put the data for the current step in a structure
selected_data(tmp_idx).centralized_right_list = centralized_right_list;
selected_data(tmp_idx).initial_position = initial_position;
selected_data(tmp_idx).polynomial = poly_right;
selected_data(tmp_idx).frame = fidx;

 % FIT YOUR FUNCTION HERE
warning off % i dont want WARNING spammed because it's not well fitted
tmp_ply_l = polyfit(pixel_list_left(:,2),pixel_list_left(:,1),5);
tmp_ply_r = polyfit(pixel_list_right2(:,2),pixel_list_right2(:,1),5);
warning on

% create dummy variables for plotting
ymp_yl = linspace(min(pixel_list_left(:,2)),max(pixel_list_left(:,2)),100);
ymp_yr = linspace(min(pixel_list_right2(:,2)),max(pixel_list_right2(:,2)),100);
y_cr=centralized_right_list(:,2);
y_crl=centralized_left_list(:,2);


% plot the polynomials
plot(polyval(tmp_ply_l,ymp_yl),ymp_yl,'r--','LineWidth',1.5)
plot(polyval(tmp_ply_r,ymp_yr),ymp_yr,'r--','LineWidth',1.5)
%plot(centralized_right_list(:,1),centralized_right_list(:,2))
pixel_list_right = pixel_list_right2;
leftSideInside = flipud(([pixel_list_left(:,1) -pixel_list_left(:,2)]+[-pixel_list_left(1,1) myPileTip])*px2m);
rightSideInside = flipud(([pixel_list_right(:,1) -pixel_list_right(:,2)]+[-pixel_list_right(1,1) myPileTip])*px2m);


figure(1000)
hold on
plot(-centralized_left_list(:,1),-centralized_left_list(:,2))
plot(-polyval(poly_left,y_crl),-y_crl)


% FIT YOUR FUNCTION HERE
warning off % i dont want WARNING spammed because it's not well fitted
fitLeftInside = polyfit(leftSideInside(:,2),leftSideInside(:,1),5);
fitRightInside = polyfit(rightSideInside(:,2),rightSideInside(:,1),5);
warning on



inside_data(fidx).fitLeftInside=fitLeftInside;
inside_data(fidx).fitRightInside=fitRightInside;
inside_data(fidx).leftSideInside=leftSideInside;
inside_data(fidx).rightSideInside=rightSideInside;
inside_data(fidx).pixel_list_left=pixel_list_left;
inside_data(fidx).pixel_list_right=pixel_list_right;
%% THESE WERE ATTEMPTS TO FIND THE EDGES USING THE DIFF FROM INITIAL FRAME, WORKS OKAYISH
% imshow(mat2gray(I31))
% imshow()
%     ref_I = frame_ref(start_pos_y:end,:,1);
% 
% I_diff = I-ref_I;
% I3 = I_diff(y_idx_mudline:myPileTip+edge_limits_pile,i1-edge_limits_pile:i2+edge_limits_pile);
% I31 = I3(edge_limits_pile:end,1:size(I3,2)/2);
% I32 = I3(edge_limits_pile:end,size(I3,2)/2:end);
% 
% 
% 
% [Gx,Gy] = imgradientxy(imsharpen(I31),'prewitt');
% imshow(mat2gray(Gx))
% [Gmag,Gdir] = imgradient(Gx,Gy);
% XD = abs(Gdir)<30;
% imshow(XD)
% EE = edge(imsharpen(I31),'canny');
% imshow(EE)
% imshowpair(Gmag,Gdir,'montage')
% title('Gradient Magnitude (Left) and Gradient Direction (Right)')
% % define a threshold using canny edge detection blabla this can be refined
% % eventually?

% calculate the standard deviation of the intensity value within the pile,
% excluding the edges of the pile itself (i.e. the position along X at the
% mudline). 


%% THIS IS NOT USED
% find the inflection between side and bottom, keep the former
% gbdy = mat2gray(gradient(smooth(pixel_list_left(:,1))));
% tmp = gbdy<0.07;
% pixel_list_left = pixel_list_left(tmp,:);
% 
% 
% gbdy = mat2gray(-gradient(smooth(pixel_list_right(:,1))));
% tmp = gbdy<0.07;
% pixel_list_right = pixel_list_right(tmp,:);


%% find the pile boundaries outside the soil / above the mudline
rp_outside.Eccentricity = 1;
ii = 0;

while abs(rp_outside.Eccentricity)>0.5
outside_I = I(1:y_idx_mudline-ii,:);
tmp = outside_I<255;
tmp = imfill(tmp,'holes');

rp_outside = regionprops(tmp,'all');

[~,i]=sort([rp_outside.Area]);
rp_outside=rp_outside(i);
rp_outside=rp_outside(end);
ii=ii+1;
% figure(3)
% imshow(rp_outside.FilledImage)
% drawnow
end


boundaries_outside = bwboundaries(rp_outside.FilledImage);
boundaries_outside=boundaries_outside{1};

left_outside_boundary = boundaries_outside((boundaries_outside(:,2)<max(boundaries_outside(:,2))/10),:);
right_outside_boundary = boundaries_outside(boundaries_outside(:,2)>9*max(boundaries_outside(:,2))/10,:);


% remove doubles on the X direction
[~,i]=sort(left_outside_boundary(:,2));
left_outside_boundary=left_outside_boundary(i,:);
[~,i]=unique(left_outside_boundary(:,1));
left_outside_boundary = left_outside_boundary(i,:);

% remove doubles on the X direction
[~,i]=sort(right_outside_boundary(:,2),'descend');
right_outside_boundary=right_outside_boundary(i,:);
[~,i]=unique(right_outside_boundary(:,1));
right_outside_boundary = right_outside_boundary(i,:);


pixel_list_left_global= pixel_list_left + [i1-edge_limits_pile-1, y_idx_mudline];
pixel_list_right_global= pixel_list_right + [i1-edge_limits_pile-1, y_idx_mudline];
left_outside_boundary_global = left_outside_boundary + [0 i1+edge_limits_pile/2-2];
right_outside_boundary_global = right_outside_boundary + [0 i1+edge_limits_pile/2-2];




vec_pos_mudline(fidx) = y_idx_mudline*px2m;



left_boundary_global = [pixel_list_left_global;left_outside_boundary_global(:,2:-1:1)];
right_boundary_global = [pixel_list_right_global;right_outside_boundary_global(:,2:-1:1)];

% sort by y
[~,i]=sort(left_boundary_global(:,2),'descend');
left_boundary_global=left_boundary_global(i,:);
[~,i]=sort(right_boundary_global(:,2),'descend');
right_boundary_global=right_boundary_global(i,:);


figure(111)
clf
imshow(I);
hold on
plot(left_boundary_global(:,1),left_boundary_global(:,2),'r-','LineWidth',2.5)
plot(right_boundary_global(:,1),right_boundary_global(:,2),'r-','LineWidth',2.5)



% 
% figure(112)
% hold on
% plot(left_boundary_global(:,1)*px2m,-left_boundary_global(:,2)*px2m,'-','Color',color,'LineWidth',0.5)
% plot(right_boundary_global(:,1)*px2m,-right_boundary_global(:,2)*px2m,'-','Color',color,'LineWidth',0.5)
% axis equal
% ylim([-0.14 -0.03])


leftSide = flipud(([left_boundary_global(:,1) -left_boundary_global(:,2)]+[-left_boundary_global(1,1) myPileTip])*px2m);
rightSide = flipud(([right_boundary_global(:,1) -right_boundary_global(:,2)]+[-right_boundary_global(1,1) myPileTip])*px2m);


% figure(1)
% clf
% plot(leftSide(:,1),leftSide(:,2))
% axis equal
% hold on
% plot(rightSide(:,1),rightSide(:,2))

% FIT YOUR FUNCTION HERE
warning off % i dont want WARNING spammed because it's not well fitted
fitLeft = polyfit(leftSide(:,2),leftSide(:,1),5);
fitRight = polyfit(rightSide(:,2),rightSide(:,1),5);
warning on


%plot(polyval(fitLeft,linspace(0,0.12,100)),linspace(0,0.12,100))

% save the coefficients at the given step
fits_left(fidx,:)=fitLeft;
fits_right(fidx,:)=fitRight;

figure(111)
% write a spiffy giffy
frame = getframe(gcf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if fidx == start_frame
    imwrite(imind,cm,filename_out,'gif', 'Loopcount',inf,'DelayTime',1/((end_frame-start_frame)/skippy/duration));
else
    imwrite(imind,cm,filename_out,'gif','WriteMode','append','DelayTime',1/((end_frame-start_frame)/skippy/duration));
end

saveas(gcf,['frames_lateral/contours_' num2str(fidx) '.png']);





% write to video
if rem(fidx,30)==0
writeVideo(outputVideo,im)
end
figure(112)
% write a spiffy giffy
frame = getframe(gcf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if fidx == start_frame
    imwrite(imind,cm,filename_out2,'gif', 'Loopcount',inf,'DelayTime',1/((end_frame-start_frame)/skippy/duration));
else
    imwrite(imind,cm,filename_out2,'gif','WriteMode','append','DelayTime',1/((end_frame-start_frame)/skippy/duration));
end

saveas(gcf,['frames_lateral/sides_only_' num2str(fidx) '.png']);


end
% oh yeah this is to plot the evolution of the domain surface
close(outputVideo)
%%%
%close all
save('selected_data.mat','selected_data');
save('boundary_fits.mat','fits_left','fits_right','vec_pos_mudline','inside_data');
