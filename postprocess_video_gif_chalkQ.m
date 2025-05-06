clearvars
close all
fname = 'Stage1_zeroed_and_continuing_seenotes_to52.25mm.avi';
filename_out = strrep(fname,'.avi','_postprocessed.gif');
filename_vid = strrep(fname,'.avi','_postprocessed.mp4');
edge_limits = 30;
edge_limits_pile = 20;
num_segments_smoothing = 50;
vr = VideoReader(fname);
fig_visibility = 'off'; % to plot without visualising the fig
ifr_plr_data = importdata('avg_ifr_plr_chalkQ_rev3190325.txt').data %.data; % Assuming the file has columns: depth, IFR, PLR
depth_file_old = ifr_plr_data(:, 1);  % Depths from file
depth_file = -ifr_plr_data(:, 1);  %updated for new ct factor issue
IFR_file = ifr_plr_data(:, 3);    % IFR values from file
PLR_file = ifr_plr_data(:, 4);    % PLR values from file

filename_vid = 'CHALK_postprocessed.mp4';


data = importdata('chalk_Q_time_disp_load (1).txt'); %updated for new ct factor issue
data=data.data;
vec_time  = data(:,1);
vec_disp  = data(:,2)/0.930958095;
vec_force  = data(:,3);

idx_start_data = 14747;
idx_end_data = 97695;

outputVideo = VideoWriter(filename_vid,'MPEG-4');
outputVideo.FrameRate = 30;
open(outputVideo)


frame_ref = read(vr,1);
imshow(frame_ref);
title('Select upper limit of the sand (go slightly above)')
%[~,start_pos_y] = ginput(1);
start_pos_y = 711;

n_frames = vr.NumFrames;%vr.Duration/vr.FrameRate;
actual_framerate = vr.NumFrames/vr.Duration;
skippy = 10;
skippy_frame = 10;

frame_end = read(vr,n_frames);
imshow(frame_end);
title('Select lower limit of the pile')
lower_limit_pile = 1421;
%[~,lower_limit_pile] = ginput(1);


start_frame = actual_framerate*0+1;
end_frame = n_frames;% actual_framerate*22;
%%
close all
    h = figure(4);
          set(h,'visible',fig_visibility);
 

for fidx = start_frame:skippy_frame:end_frame
    disp(['Frame:' num2str(fidx) '/' num2str(end_frame)])
local_frame = read(vr,fidx);


%binarize the entire domain
I = local_frame(start_pos_y:end,:,1);
ii = I<125;%255;
I_diff = I - frame_ref(start_pos_y:end,:,1);
% figure(101)
% imshow(mat2gray(I_diff));
I2 = double(I).*ii;

%extract the biggest region (the soil)
rp = regionprops(ii,'all');
[~,i]=sort([rp.Area]);
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
boundary = B{1};



% find the left and right limits of the pile
i1 = 939;%max([find(ii(1,:)==1,1,'first')-edge_limits_pile,1]);
i2 = 1087;%min([find(ii(1,:)==1,1,'last')+edge_limits_pile,size(ii,2)]);
% find the left and right limits of the soil
i3 = 517;%bb(1)+edge_limits;
i4 = 1501;%bb(1)+bb(3)-edge_limits;

if fidx == 1

% to plot 2.5D map of the top displacement
static_i3 = i3;
static_i4 = i4;
end


%split the soil boundaries to include/exclude piles and edges
idx_boundary_pile = sum([bb(1)+boundary(:,2)>i1,bb(1)+boundary(:,2)<i2,boundary(:,1)<mean(boundary(:,1))],2)==3;
idx_boundary = (sum([~idx_boundary_pile,bb(1)+boundary(:,2)>i3,bb(1)+boundary(:,2)<i4,boundary(:,1)<mean(boundary(:,1))],2)==4);
    


% extract the pile plug
pile_center = round(mean([i2,i1]));

vec_val_pile_center = I(:,pile_center);
vec_smooth_val_pile_center = smoothdata(vec_val_pile_center);
vec_g_smooth_pile = gradient(vec_smooth_val_pile_center);
min_gradient_smooth = find(vec_g_smooth_pile==min(vec_g_smooth_pile),1);
y_idx_bottom_pile = find(vec_g_smooth_pile(min_gradient_smooth:end)>0,1)+min_gradient_smooth;
y_idx_top_plug = find(I(:,pile_center)<I(1,pile_center),1);
y_idx_bottom_plug = 140;%boundary(find(idx_boundary_pile,1),1);
plug_thickness = y_idx_bottom_plug-y_idx_top_plug;
% the tip is found as the peak in the gradient of value along the pile depth
tmp = mat2gray(double(I(y_idx_bottom_plug-10:end,i1:i2)));

tmp2 = gradient(smooth(mean(tmp,2)));
tmp3 = find(min(tmp2)==tmp2);
tmp4 = find(tmp2(tmp3:end)==max(tmp2(tmp3:end)),1)+tmp3-4;

if fidx == 1
myPileTip=y_idx_bottom_plug;
elseif abs(tmp4+y_idx_bottom_plug-10-myPileTip)<30
myPileTip = tmp4+y_idx_bottom_plug-10;
end

tmp = (double(I(1:y_idx_bottom_plug-25,i1:i2)));
tmp1 = smoothdata(mean(tmp,2));
tmp2 = gradient(smooth(mean(tmp,2)));
if fidx == 1
myPlug = y_idx_bottom_plug-20;
elseif abs(find(tmp2==min(tmp2))-myPlug)<40
myPlug = find(tmp2==min(tmp2));%else
%
end

vec_avg_value = mat2gray(mean((double(I(y_idx_bottom_plug+edge_limits*3:end,i1+edge_limits/2:i2-edge_limits/2))),2));
img_plot_plug = mat2gray(((double(I(y_idx_bottom_plug+edge_limits*3:end,i1:i2)))));

%myPileTip = y_idx_bottom_pile;
gv = (mat2gray(smoothdata(gradient(vec_avg_value),'movmean',length(vec_avg_value)/num_segments_smoothing)));
y_idx_plug_inside = [];
thrs = 0.2;
while and(length(y_idx_plug_inside)==0, sum(vec_avg_value) >0)
[~,y_idx_plug_inside] = findpeaks(-gv,'MinPeakProminence',thrs);
thrs = thrs-thrs*0.1;
end

%y_idx_plug_inside = start_gv+find(gv(start_gv:myPileTip-y_idx_bottom_plug-edge_limits)==min(gv(start_gv:myPileTip-y_idx_bottom_plug-edge_limits)),1);


%if 1
y_idx_plug_inside=y_idx_plug_inside(1);





idx_start_data = 14747;
idx_end_data = 97695;
range_data = (idx_end_data-idx_start_data);
max_disp_data = vec_disp(idx_end_data);
min_disp_data = vec_disp(idx_start_data);
start_soil_y = 126;
curr_disp_xct = (myPileTip - start_soil_y)/(lower_limit_pile - start_pos_y);
myStep = round(curr_disp_xct*range_data+idx_start_data);

px2mm = 0.089972/ 0.9002;

pilePen(fidx) = (myPileTip - start_soil_y) * px2mm;
plugHeight(fidx) = (myPlug-y_idx_bottom_plug+20)*px2mm;
pilePen2=pilePen;
plugHeight2=plugHeight;
plugHeight2(pilePen==0)=[];
pilePen2(pilePen==0)=[];

% Remove duplicate depth values for interpolation
[depth_file_unique, unique_idx] = unique(depth_file);

% Also remove corresponding IFR and PLR values at those indices
IFR_file_unique = IFR_file(unique_idx);
PLR_file_unique = PLR_file(unique_idx);

% Now interpolate using the cleaned data
IFR_interp = interp1(depth_file_unique, IFR_file_unique, pilePen2, 'linear', 'extrap');
PLR_interp = interp1(depth_file_unique, PLR_file_unique, pilePen2, 'linear', 'extrap');


%%

map2 = brewermap(4, 'Set1');
% Within the loop:
if rem(fidx-1, skippy) == 0


    clf
    set(gcf, 'Position', get(0, 'Screensize'))

    tiledlayout(1,4)

    nexttile([1 1])
    imshow(mat2gray(I));
    hold on
    plot([i1 i2]-4, [myPlug myPlug], 'r-.', 'LineWidth', 2)
    plot([i1 i2]-4, [myPileTip myPileTip], 'LineStyle', '--', 'Color', map2(:,3), 'LineWidth', 2)
    box off
    xlim([i3-30 i4+30])
    set(gca, 'FontSize', 14)

    nexttile % Second column
    plot((-vec_force(idx_start_data:myStep)), (vec_disp(idx_start_data:myStep)), 'LineStyle', '--', 'Color', map2(:,3), 'LineWidth', 3)
    ylim([-100 0]*1.05);
    yticks([ -105 -90 -75 -60 -45 -30 -15 0])
    xlim([min(-vec_force) max(-vec_force)]*1.1);
    xlabel('Pile insertion force, {\it Q} [kN]', 'FontSize', 14);
    ylabel('Pile insertion below the mudline, {\it z} [mm]', 'FontSize', 14);
    set(gca, 'FontSize', 14)
     box off
    nexttile % Third column
    plot(-plugHeight2, -pilePen2, 'r-.', 'LineWidth', 3);
    ylim([-100 0]*1.05);
     yticks([ -105 -90 -75 -60 -45 -30 -15 0])
     box off
    xlim([-20 20]);
    xlabel('Plug height above the mudline, {\it H_p}  [mm]', 'FontSize', 14);
    set(gca, 'FontSize', 14)

    set(gcf, 'color', 'w');
         nexttile
        hold on
        plot(IFR_interp, -pilePen2, '-k', 'LineWidth', 2, 'DisplayName', 'IFR')
        plot(PLR_interp, -pilePen2, '--k', 'LineWidth', 2, 'DisplayName', 'PLR')
            ylim([-100 0]*1.05);
         yticks([ -105 -90 -75 -60 -45 -30 -15 0])
        xlim([0 1.5]);
             box off
        %set(gca, 'YDir', 'reverse')  % Depth should increase downwards
        xlabel('IFR or PLR [-]', 'FontSize', 14)
        %ylabel('Pile insertion below the mudline, z [mm]', 'FontSize', 14)
        %legend('IFR', 'PLR', 'Location', 'northeast')
        set(gca, 'FontSize', 14)
        set(gcf, 'color', 'w');
    
    % Capture the frame and convert to indexed image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write the first frame
    if fidx == start_frame
        imwrite(imind, cm, 'force_disp_install_3.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        % Append subsequent frames
        imwrite(imind, cm, 'force_disp_install_3.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
end
%%

fname = 'Stage1_zeroed_and_continuing_seenotes_to103.648mm.avi';
edge_limits = 30;
edge_limits_pile = 20;
num_segments_smoothing = 50;
vr = VideoReader(fname);


idx_start_data = 103487;
idx_end_data = 97695;
time_pile_level = 16*60*actual_framerate;
time_pile_restart = (18*60)*actual_framerate;
n_frames1=1920;%length(pilePen)+skippy_frame-1;
myPileTip = pilePen2(n_frames1/skippy_frame)/px2mm+start_soil_y;

%start_pos_y = 711;

n_frames = vr.NumFrames;%vr.Duration/vr.FrameRate;
actual_framerate = vr.NumFrames/vr.Duration;


frame_end = read(vr,n_frames);
imshow(frame_end);
title('Select lower limit of the pile')
lower_limit_pile = 1957;
plug_below = 1334;
%[~,lower_limit_pile] = ginput(1);


start_frame = actual_framerate*0+1;
end_frame = n_frames;% actual_framerate*22;

close all
frame_ref = read(vr,1);

clear myDistance

for fidx = start_frame:skippy_frame:end_frame
    disp(['Frame:' num2str(fidx) '/' num2str(end_frame)])

local_frame = read(vr,fidx);


%binarize the entire domain
I = local_frame(start_pos_y:end,:,1);
ii = I<125;%255;
I_diff = I - frame_ref(start_pos_y:end,:,1);
% figure(101)
% imshow(mat2gray(I));
I2 = double(I).*ii;

%extract the biggest region (the soil)
rp = regionprops(ii,'all');
[~,i]=sort([rp.Area]);
rp = rp(i);
rp = rp(end);

%extract its bundaries
[B,L] = bwboundaries(rp.FilledImage);
bb = rp.BoundingBox;
boundary = B{1};



% find the left and right limits of the pile
i1 = 939;%max([find(ii(1,:)==1,1,'first')-edge_limits_pile,1]);
i2 = 1087;%min([find(ii(1,:)==1,1,'last')+edge_limits_pile,size(ii,2)]);
% find the left and right limits of the soil
i3 = 517;%bb(1)+edge_limits;
i4 = 1501;%bb(1)+bb(3)-edge_limits;


% the tip is found as the peak in the gradient of value along the pile depth
tmp = mat2gray(double(I(y_idx_bottom_plug-10:end,i1:i2)));

tmp2 = gradient(smooth(mean(tmp,2)));
tmp3 = find(min(tmp2)==tmp2);
tmp4 = find(tmp2(tmp3:end)==max(tmp2(tmp3:end)),1)+tmp3-4;

 if abs(tmp4+y_idx_bottom_plug-myPileTip)<30
myPileTip = tmp4+y_idx_bottom_plug;
 end

if fidx > time_pile_level && fidx < time_pile_restart

    myPlug = myPileTip-myDistance;

elseif fidx < time_pile_level

tmp = (double(I(1:y_idx_bottom_plug,i1:i2)));
tmp2 = gradient(smooth(mean(tmp,2)));

if abs(find(tmp2==min(tmp2),1)-myPlug)<40
myPlug = find(tmp2==min(tmp2));%else
myDistance = myPileTip-myPlug;
end


elseif fidx > time_pile_restart

tmp = mat2gray(double(I(y_idx_bottom_plug:myPileTip,i1+20:i2-20)));
tmp2 = gradient(smooth(mean(tmp,2)));
tmp3 = find(tmp2>0,1);
tmp4 = tmp2(tmp3:end);
tmp5 = find(tmp4==min(tmp4),1);
tmp6 = tmp5 + tmp3 + y_idx_bottom_plug;
if ~isnan(myPlug)
if abs(myPlug-tmp6)<50
myPlug = tmp6;
end
else
myPlug = 28 + y_idx_bottom_plug;
end

end
% else
% 
% if ~exist('myDistance')
%     myDistance = myPileTip+myPlug-y_idx_bottom_plug; 
% end
% 
% tmp = mat2gray(double(I(y_idx_bottom_plug:myPileTip,i1+20:i2-20)));
% tmp1 = smooth(mean(tmp,2));
% tmp2 = gradient(tmp1);
% if abs(find(tmp2==min(tmp2))-myPlug+y_idx_bottom_plug-1)<50
% myPlug = myPileTip-myDistance;%find(tmp2==min(tmp2))+y_idx_bottom_plug-1;%else
% 
% end




vec_avg_value = mat2gray(mean((double(I(y_idx_bottom_plug+edge_limits*3:end,i1+edge_limits/2:i2-edge_limits/2))),2));
img_plot_plug = mat2gray(((double(I(y_idx_bottom_plug+edge_limits*3:end,i1:i2)))));

%myPileTip = y_idx_bottom_pile;
gv = (mat2gray(smoothdata(gradient(vec_avg_value),'movmean',length(vec_avg_value)/num_segments_smoothing)));
y_idx_plug_inside = [];
thrs = 0.2;
while and(length(y_idx_plug_inside)==0, sum(vec_avg_value) >0)
[~,y_idx_plug_inside] = findpeaks(-gv,'MinPeakProminence',thrs);
thrs = thrs-thrs*0.1;
end

%y_idx_plug_inside = start_gv+find(gv(start_gv:myPileTip-y_idx_bottom_plug-edge_limits)==min(gv(start_gv:myPileTip-y_idx_bottom_plug-edge_limits)),1);


%if 1
y_idx_plug_inside=y_idx_plug_inside(1);




idx_start_data = 14747;
idx_end_data = length(vec_force);
range_data = (idx_end_data-idx_start_data);
max_disp_data = vec_disp(idx_end_data);
min_disp_data = vec_disp(idx_start_data);
start_soil_y = 126;
curr_disp_xct = (myPileTip - start_soil_y)/(lower_limit_pile - start_soil_y);
myStep = round(curr_disp_xct*range_data+idx_start_data);

px2mm = 0.089972/ 0.9002;


pilePen_tmp = (myPileTip - start_soil_y) * px2mm;
myStep = find(pilePen_tmp<-vec_disp,1);


pilePen(end+1) = (myPileTip - start_soil_y) * px2mm;
plugHeight(end+1) = (myPlug-y_idx_bottom_plug+20)*px2mm;
pilePen2=pilePen;
plugHeight2=plugHeight;
plugHeight2(pilePen==0)=[];
pilePen2(pilePen==0)=[];

% Interpolate IFR and PLR values to match the depth of pilePen2
% Now interpolate using the cleaned data
IFR_interp2 = interp1(depth_file_unique, IFR_file_unique, pilePen2, 'linear', 'extrap');
PLR_interp2 = interp1(depth_file_unique, PLR_file_unique, pilePen2, 'linear', 'extrap');
if  rem(fidx-1,skippy)==0

h = figure(4);
    set(h,'visible',fig_visibility);

clf
set(gcf, 'Position', get(0, 'Screensize'))
tiledlayout(1,4) 


    nexttile([1 1])
    imshow(mat2gray(I));
    hold on
    plot([i1 i2]-4, [myPlug myPlug], 'r-.', 'LineWidth', 2)
    plot([i1 i2]-4, [myPileTip myPileTip], 'LineStyle', '--', 'Color', map2(:,3), 'LineWidth', 2)
 box off
    xlim([i3-30 i4+30])
    set(gca, 'FontSize', 14)

    nexttile % Second column
    plot((-vec_force(idx_start_data:myStep)), (vec_disp(idx_start_data:myStep)), 'LineStyle', '--', 'Color', map2(:,3), 'LineWidth', 3)
    ylim([-100 0]*1.05);
    yticks([ -105 -90 -75 -60 -45 -30 -15 0])
    xlim([min(-vec_force) max(-vec_force)]*1.1);
    xlabel('Pile insertion force, {\it Q} [kN]', 'FontSize', 24);
    ylabel('Pile insertion below the mudline, {\it z} [mm]', 'FontSize', 14);
    set(gca, 'FontSize',14)
 box off
    nexttile % Third column
    plot(-plugHeight2, -pilePen2, 'r-.', 'LineWidth', 3);
    ylim([-100 0]*1.05);
    yticks([ -105 -90 -75 -60 -45 -30 -15 0])
 box off
    xlim([-20 20]);
    xlabel('Plug height above the mudline, {\it H_p}  [mm]', 'FontSize', 14);
    set(gca, 'FontSize', 14)

    set(gcf, 'color', 'w');

    % Third tile: IFR and PLR vs Depth (interpolated to match pilePen2)
        nexttile
        hold on
        plot(IFR_interp2, -pilePen2, '-k', 'LineWidth', 2, 'DisplayName', 'IFR')
        plot(PLR_interp2, -pilePen2, '--k', 'LineWidth', 2, 'DisplayName', 'PLR')
            ylim([-100 0]*1.05);
         yticks([ -105 -90 -75 -60 -45 -30 -15 0])
        xlim([0 1.5]);
        %set(gca, 'YDir', 'reverse')  % Depth should increase downwards
        xlabel('IFR or PLR [-]', 'FontSize', 14)
       % ylabel('Pile insertion below the mudline, z [mm]', 'FontSize', 14)
       % legend('IFR', 'PLR', 'Location', 'northeast')
        set(gca, 'FontSize', 14)
        set(gcf, 'color', 'w');
         box off
% frame = getframe(gcf);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% 
% writeVideo(outputVideo,im)

% 
% frame = getframe(gcf);
% im = frame2im(frame);
%    writeVideo(outputVideo,im)
saveas(gcf,['chalk_in/' num2str(fidx) '.png'])
[imind,cm] = rgb2ind(im,256);

    imwrite(imind,cm,'force_disp_install_3.gif','gif','WriteMode','append','DelayTime',0.1);
end

end
%%
% Final Static Figure
% h = figure(5);
% set(h, 'visible', fig_visibility);
%  
% % Clear figure and set fullscreen
% clf
figure1 = figure;
set(figure1, 'Position', [100, 100, 860, 400]);
 hold on
% Second Tile: Pile Insertion Force vs. Depth
subplot(1,3,1)
 hold on
plot(-vec_force(idx_start_data:myStep), vec_disp(idx_start_data:myStep), '--', 'Color', map2(:, 3), 'LineWidth', 3);
ylim([-110 0] * 1.05);
yticks([-105 -90 -75 -60 -45 -30 -15 0]);
xlim([min(-vec_force) max(-vec_force)] * 1.1);
xlabel('$F_z$ [kN]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$z$ [mm]', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
box off;
 
% Third Tile: Plug Height above Mudline
subplot(1,3,2)
 hold on
plot(-plugHeight2, -pilePen2, 'r-.', 'LineWidth', 3);
ylim([-110 0] * 1.05);
yticks([-105 -90 -75 -60 -45 -30 -15 0]);
xlim([-20 20]);
xlabel('$H_p$ [mm]', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
box off;
 
subplot(1,3,3)
 hold on
plot(IFR_interp2, -pilePen2, '-k', 'LineWidth', 2, 'DisplayName', 'IFR');
plot(PLR_interp2, -pilePen2, '--k', 'LineWidth', 2, 'DisplayName', 'PLR');
%plot(1.09 * PLR_interp2 - 0.22, -pilePen2, ':r', 'LineWidth', 2, 'DisplayName', 'IFR = PLR * 1.09 - 0.22');
 
% Adding Paik & Lee 2003 Model
plot(PLR_interp2 * 1.09 -0.22, -pilePen2, ':b', 'LineWidth', 2, 'DisplayName', 'Paik \& Lee 2003');
 
ylim([-110 0] * 1.05);
yticks([-105 -90 -75 -60 -45 -30 -15 0]);
xlim([0 1.5]);
xlabel('IFR or PLR [-]', 'FontSize', 12, 'Interpreter', 'latex');
legend('IFR', 'PLR', "Paik" +newline+ "& Lee (2003)", 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
set(gca, 'FontSize', 12);
box off;
 %%
set(gcf, 'color', 'w');
 
% ---- Save data to CSV file ----
% Prepare data for saving with two tables
output_data = table();
output_data2 = table();
 
% Second Tile Data: Pile Insertion Force vs. Depth
output_data.PileInsertionForce_Q_kN = -vec_force(idx_start_data:myStep);
output_data.PileInsertionDepth_z_mm = vec_disp(idx_start_data:myStep);
 
% Third and Fourth Tile Data: Plug Height, Pile Penetration, IFR, PLR, IFR Prediction
output_data2.PlugHeight_Hp_mm = -plugHeight2;
output_data2.PilePenetrationDepth_z_mm = -pilePen2;
output_data2.IFR = IFR_interp2;
output_data2.PLR = PLR_interp2;
output_data2.IFR_Predicted = 1.09 * PLR_interp2 - 0.22;
output_data2.Depth_mm = -pilePen2;
 
% Write tables to CSV
writetable(output_data, 'figure_data_output.csv');
writetable(output_data2, 'figure_data_forIFR_PLR_etc.csv');
disp('Data has been saved to "figure_data_output.csv" and "figure_data_forIFR_PLR_etc.csv".');
 
% Save figure as PNG
print('Fig_ICE_PICK_ChalkQ_IFR_Continuous_rev2200325', '-dpng', '-r500');
disp('Figure saved as "Fig_ICE_PICK_ChalkQ_IFR_Continuous_rev2200325.png"');
%%
close(outputVideo)

IR = plugHeight2./plugHeight2;
ff = interp1(1:myStep-idx_start_data+1,vec_force(idx_start_data:myStep),1:length(IR));
ff2 = downsample(vec_force(idx_start_data:myStep),ceil((myStep-idx_start_data+1)/length(IR)))
plot(-IR,-ff2)
xlabel('Plug Height / Penetration ')
ylabel('Force [kN]')
saveas(gcf,'plug_vs_force.png');

figure
IR = plugHeight2./pilePen2;
ff = interp1(1:myStep-idx_start_data+1,vec_force(idx_start_data:myStep),1:length(IR));
ff2 = downsample(vec_force(idx_start_data:myStep),ceil((myStep-idx_start_data+1)/length(IR)))
IR2 = diff(smooth(plugHeight2))./diff(smooth(pilePen2));
plot(-IR2,-ff2(2:end))
xlabel('\partial Plug Height / \partial Penetration ')
ylabel('Force [kN]')
saveas(gcf,'plug_ratio_vs_force.png');
figure
SIR = smooth(IR);
DSIR = gradient(SIR);
plot(-DSIR,-ff2)
xlabel('\partial (Plug Height / Penetration) ')
ylabel('Force [kN]')
saveas(gcf,'delta_plug_ratio_vs_force.png');
plugpos=plugHeight2'
tippos=pilePen2'

