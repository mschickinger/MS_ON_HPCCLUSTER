% What's new in version 10b?

% Gaussian fitting procedure is (or rather was already in v10a)
% adapted for less data communication overhead -> supposedly more efficient
% working order is changed to movie -> frame -> spots
% mapping of fitted positions is drawn out of parallel loop

% v10b -> batch job assignment for parallel computation on clusters 

% Also includes drift (x,y) correction routine (optional) and new
% intensity tracer.

%% startup
clc, clear all, close all
path0 = cd;
run('/nfs/matlabuser/matthiasschickinger/MATLAB/my_prefs.m')


%% choose colors
rgb={'red','green','blue'};
[colors,ok]=listdlg('PromptString', 'Select two colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
while ne(length(colors),2) && ok>0
    [colors,ok]=listdlg('PromptString', 'Select _TWO_ colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
end

channel = cell(2,1);
channel{1} = rgb{colors(1)};
channel{2} = rgb{colors(2)};

[chb,ok]=listdlg('PromptString', 'Which one is surface-bound?',...
                'ListString', channel, 'SelectionMode', 'single',...
                'OKString', 'Confirm');           
channel_bound=rgb{chb};

%% LOAD STACK OF MOVIES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_ch1 = pickFirstFitsFiles(pname, channel{1}); 
files_ch2 = pickFirstFitsFiles(pname, channel{2});

N_movie = length(files_ch1);
if length(files_ch1) ~= length(files_ch2)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)

%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', ['Sequence ' channel{1} ':'], ['Sequence ' channel{2} ':'],... % sample options
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Minimal length [frames]:',...
    'Average over first N_frames:'};
input_default = {'2', '-1', '1', '1', '4', '3', '20', '100'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence_ch1 = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence_ch1(1,i) =1;
    end
end
sequence_ch2 = zeros(1, size(tmp{4},2));
for i=1:size(tmp{4},2)
    if(tmp{4}(i) == '1')
        sequence_ch2(1,i) =1;
    end
end
r_find = str2double(tmp(5)); % radius used to find spots
r_integrate = str2double(tmp(6)); % radius used for integration of intensities
min_length = str2double(tmp(7)); % minimal number of found spots in a trace
N_frames = str2double(tmp(8)); % in get_h_min, average over first N_frames is used

%% generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    ch2{i} = movie(pname, files_ch2{i}, first, last, sequence_ch2); % pname, fname, first, last, sequence
end

%%
button = questdlg(['Map positions ' channel{1} ' ON ' channel{2} ' and vice versa?'],'Mapping','Yes','No','No');
mapping = strcmp(button, 'Yes');

button = questdlg('Perform drift correction?','Drift correction','Yes','No','Yes');
drift_cor = strcmp(button, 'Yes');
%%
if mapping
    [mapping_file_1TO2, mapping_dir]=uigetfile(data_dir,['Choose the ' channel{1} '2' channel{2} ' mapping file:']);
    map1TO2 =load([mapping_dir mapping_file_1TO2], 'tform');
    tform_1TO2 = map1TO2.tform;
    display(['loaded ' channel{1} ' TO ' channel{2} ' mapping file: ' mapping_dir mapping_file_1TO2]);
    
    [mapping_file_2TO1]=uigetfile(mapping_dir,['Choose the ' channel{2} '2' channel{1} ' mapping file:']);
    map2TO1=load([mapping_dir mapping_file_2TO1]);
    tform_2TO1 = map2TO1.tform; %['tform_' channel{2} 'ON' channel{1}];
    display(['loaded ' channel{2} ' TO ' channel{1} ' mapping file: ' mapping_dir mapping_file_2TO1]);
end

%% compute average images
avg_img = cell(N_movie, 4);

for i=1:N_movie
    avg_img{i, 1} = ch1{i}.average_image(N_frames);
    avg_img{i, 2} = ch2{i}.average_image(N_frames);
    avg_img{i, 3} = ch1{i}.average_image_last(N_frames); % for fitting threshold assignment, drift correction
    avg_img{i, 4} = ch2{i}.average_image_last(N_frames); % for fitting threshold assignment, drift correction
end

%% get threshold and find peaks from first N_frames

peaks_raw = zeros(0,5);
all_positions = cell(N_movie, 2);

if mapping
pos1on2 = cell(N_movie, 1);
pos2on1 = cell(N_movie, 1);
end

for i=1:N_movie 
    [h_min, pos_ch1] = ch1{i}.get_h_min(r_find, N_frames);
    [h_min, pos_ch2] = ch2{i}.get_h_min(r_find, N_frames);
    all_positions{i,1} = pos_ch1;
    all_positions{i,2} = pos_ch2;
    
    if mapping
        pos1on2{i} = transformPointsInverse(tform_2TO1, pos_ch1(:,1:2));  %%this takes coords in ch1 and transforms them to coords in ch2
        pos2on1{i} = transformPointsInverse(tform_1TO2, pos_ch2(:,1:2));  %%this takes coords in ch2 and transforms them to coords in ch1
    end
    
    % map peaks
    trace_map = map_traces(pos_ch1(:,1:2), pos_ch2(:,1:2), pos_ch2(:,1:2), r_find*2)+1; %map the traces from average positions

    tmp = zeros(size(trace_map,1),5);
    
    % combine pairs
    for j=1:size(trace_map,1)
        tmp(j,:) = [pos_ch1(trace_map(j,1), 1:2)+1 pos_ch2(trace_map(j,2), 1:2)+1 i]; %x_1 y_1 x_2 y_2 frame
    end
    
    peaks_raw = [peaks_raw; tmp];
end

N_peaks_raw = size(peaks_raw,1);
display(['You have ' num2str(N_peaks_raw) ' pairs'])

%% drift correction
if drift_cor
display('Perform drift correction now or hit enter to proceed')
pause
end

%% Fit psf to spots
s_x = 2.5;
s_y = 2.5;
w_fit = 8;

ch1_fit_raw = zeros(N_peaks_raw, 7); 
ch1_fit_err_raw = zeros(N_peaks_raw, 7);
ch2_fit_raw = zeros(N_peaks_raw, 7); 
ch2_fit_err_raw = zeros(N_peaks_raw, 7);

h = waitbar(0,'Fitting spots.. please wait');

for i=1:N_peaks_raw 
    
    % channel 1
    x1 = round(peaks_raw(i,1));
    y1 = round(peaks_raw(i,2));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x1, y1, s_x, w_fit, avg_img{peaks_raw(i, 5),1});
    ch1_fit_raw(i,:) = c;
    ch1_fit_err_raw(i,:) = c_err;

    % channel 2
    x2 = round(peaks_raw(i,3));
    y2 = round(peaks_raw(i,4));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x2, y2, s_x, w_fit, avg_img{peaks_raw(i, 5),2});
    ch2_fit_raw(i,:) = c;
    ch2_fit_err_raw(i,:) = c_err;
    
    waitbar( i/N_peaks_raw , h, ['Fitting spot... ' num2str(i) ' of ' num2str(N_peaks_raw) ' done']) % update waitbar
end

close(h)

%% SORT OUT: remove spots where ratio of width is not close to 1 and which are too large

criteria = ones(N_peaks_raw,2);

if chb == 1
criteria(:,1:2) = filter_spots(ch1_fit_raw(:,3:4), [0.9 10/9], 2);
elseif chb == 2
criteria(:,1:2) = filter_spots(ch2_fit_raw(:,3:4), [0.9 10/9], 2);
end

accepted = [criteria(:,1) & criteria(:,2)];

%remove not-accepted spots
ch1_fit = ch1_fit_raw(accepted==1, :);
ch1_fit_err = ch1_fit_err_raw(accepted==1, :);
ch2_fit = ch2_fit_raw(accepted==1, :);
ch2_fit_err = ch2_fit_err_raw(accepted==1, :);
peaks = peaks_raw(accepted==1, :);
peaks = [ch1_fit(:,1:2) ch2_fit(:,1:2) peaks(:,5)]; % use fitted poistions for further analysis


plot_discarded = strcmp(questdlg('Plot discarded spots?','Plot discarded','Yes','No','No'), 'Yes');
if plot_discarded
    path_out_discarded = [path_out filesep 'discarded'];
    mkdir(path_out_discarded)
end

plot_accepted = strcmp(questdlg('Plot accepted spots?','Plot accepted','Yes','No','No'), 'Yes');
if plot_accepted
    path_out_accepted= [path_out filesep 'accepted'];
    mkdir(path_out_accepted)
end

close all
fig_dim =1*[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
colormap gray
w_plot = 10;


if plot_discarded
    display('Plotting discarded spots...')
    if chb == 1
    for i=1:N_peaks_raw
        if ~accepted(i) % discarded spot
            message = {['Sigma ratio ' channel_bound ' OK: ' num2str(ch1_fit_raw(i,3)./ch1_fit_raw(i,4))],...
                ['Spotsize ' channel_bound ' OK: ' num2str(sqrt(ch1_fit_raw(i,3).^2+ch1_fit_raw(i,4).^2))]};
            if criteria(i,1)==0
                message{1} = ['Sigma ratio ' channel_bound ' BAD: ' num2str(ch1_fit_raw(i,3)./ch1_fit_raw(i,4))];
            end
            if criteria(i,2)==0
                message{2} = ['Spotsize ' channel_bound 'BAD: ' num2str(sqrt(ch1_fit_raw(i,3).^2+ch1_fit_raw(i,4).^2))];
            end
           
            x_1 = ch1_fit_raw(i,1);
            y_1 = ch1_fit_raw(i,2);

            plot_subframe(avg_img{peaks_raw(i, 5), 1}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'g.')
            ellipse(ch1_fit_raw(i,3), ch1_fit_raw(i,4), -ch1_fit_raw(i,5), x_1, y_1, channel{1});
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel_bound ' channel'], message{1},message{2}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_' num2str(i) '.tif'])
        end 
    end
    else
    for i=1:N_peaks_raw
        if ~accepted(i) % discarded spot
            message = {['Sigma ratio ' channel_bound ' OK: ' num2str(ch2_fit_raw(i,3)./ch2_fit_raw(i,4))],...
                ['Spotsize ' channel_bound ' OK: ' num2str(sqrt(ch2_fit_raw(i,3).^2+ch2_fit_raw(i,4).^2))]};
            if criteria(i,1)==0
                message{1} = ['Sigma ratio ' channel_bound ' BAD: ' num2str(ch2_fit_raw(i,3)./ch2_fit_raw(i,4))];
            end
            if criteria(i,2)==0
                message{2} = ['Spotsize ' channel_bound 'BAD: ' num2str(sqrt(ch2_fit_raw(i,3).^2+ch2_fit_raw(i,4).^2))];
            end
           
            x_1 = ch2_fit_raw(i,1);
            y_1 = ch2_fit_raw(i,2);

            plot_subframe(avg_img{peaks_raw(i, 5), 2}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'g.')
            ellipse(ch2_fit_raw(i,3), ch2_fit_raw(i,4), -ch2_fit_raw(i,5), x_1, y_1, channel{2});
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel_bound ' channel'], message{1},message{2}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_' num2str(i) '.tif'])
        end 
    end
    end
end

if plot_accepted
    display('Plotting accepted spots...')
    if chb == 1
    for i=1:N_peaks_raw
        if  accepted(i) % discarded spot
            message = {['Sigma ratio ' channel_bound ' OK: ' num2str(ch1_fit_raw(i,3)./ch1_fit_raw(i,4))],...
                ['Spotsize ' channel_bound ' OK: ' num2str(sqrt(ch1_fit_raw(i,3).^2+ch1_fit_raw(i,4).^2))]};
           
            x_1 = ch1_fit_raw(i,1);
            y_1 = ch1_fit_raw(i,2);

            plot_subframe(avg_img{peaks_raw(i, 5), 1}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'g.')
            ellipse(ch1_fit_raw(i,3), ch1_fit_raw(i,4), -ch1_fit_raw(i,5), x_1, y_1, channel{1});
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel_bound ' channel'], message{1},message{2}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_accepted filesep 'Accepted_' num2str(i) '.tif'])
        end 
    end
    else
    for i=1:N_peaks_raw
        if  accepted(i) % discarded spot
            message = {['Sigma ratio ' channel_bound ' OK: ' num2str(ch2_fit_raw(i,3)./ch2_fit_raw(i,4))],...
                ['Spotsize ' channel_bound ' OK: ' num2str(sqrt(ch2_fit_raw(i,3).^2+ch2_fit_raw(i,4).^2))]};
           
            x_1 = ch2_fit_raw(i,1);
            y_1 = ch2_fit_raw(i,2);

            plot_subframe(avg_img{peaks_raw(i, 5), 2}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'g.')
            ellipse(ch2_fit_raw(i,3), ch2_fit_raw(i,4), -ch2_fit_raw(i,5), x_1, y_1, channel{2});
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel_bound ' channel'], message{1},message{2}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_accepted filesep 'Accepted_' num2str(i) '.tif'])
        end 
    end
    end
end

display(['Accepted ' num2str(sum(accepted)) ' spots.'])
display(['Discarded ' num2str(sum(~accepted)) ' spots.'])
close all
N_peaks = size(peaks,1);

%% Get intensity traces 'itraces' plus median filtered itraces
display('Getting intensity traces... please wait')
tic
merged_itraces = cell(N_movie,5);
iEndval = cell(N_movie,2);
iEndval_sorted = cell(N_movie,2);
avg_iEndval = zeros(N_movie,2);
avg_iEndval_tenth = zeros(N_movie,2);
for i=1:N_movie 
    %get fluorescence intensity traces from position    
    ch1_itraces_full = ch1{i}.int_spots_in_frames(1:length(ch1{i}.frames), peaks(peaks(:,5)==i,1:2), r_integrate);
    display(['Tracing ' channel{1} ' channel in movie #' num2str(i) ' done'])
    toc %
    ch2_itraces_full = ch2{i}.int_spots_in_frames(1:length(ch2{i}.frames), peaks(peaks(:,5)==i,3:4), r_integrate);
    display(['Tracing ' channel{2} ' channel in movie #' num2str(i) ' done'])
    toc %
    movnumber = cell(size(ch1_itraces_full));
    movnumber(:) = {i};
    merged_itraces{i,1} = ch1_itraces_full;
    merged_itraces{i,2} = ch2_itraces_full;
    merged_itraces{i,3} = movnumber;
    
    %add median filtered itraces and average intensity values (over 100 frames) at the end
    tmp = cell(length(ch1_itraces_full),2);
    iEndval{i,1} = zeros(size(ch1_itraces_full,1),1);
    iEndval{i,2} = zeros(size(ch2_itraces_full,1),1);
    for j=1:length(ch1_itraces_full)
    tmp{j,1} = medfilt1(ch1_itraces_full{j}(:,4),20);
    tmp{j,2} = medfilt1(ch2_itraces_full{j}(:,4),20);
    iEndval{i,1}(j) = mean(ch1_itraces_full{j}(end-100:end,4));
    iEndval{i,2}(j) = mean(ch2_itraces_full{j}(end-100:end,4));
    end
    merged_itraces{i,4} = tmp;
    for ch = 1:2
    avg_iEndval(i,ch) = mean(iEndval{i,ch});
    [tmp_val, tmp_spot] = sort(iEndval{i,ch});
    iEndval_sorted{i,ch} = [tmp_spot, tmp_val];
    avg_iEndval_tenth(i,ch) = mean(iEndval_sorted{i,ch}(1:ceil(length(iEndval_sorted{i,ch})/10),2));
    end
end
display('itraces complete')
toc %
%% Determine fitting parameters
fit_cutoff = cell(N_movie,2);

%%
fc = figure('PaperPositionMode', 'manual', 'PaperUnits', 'centimeters',...
'PaperPosition', [0 0 25 9]);

for m = 1:N_movie
    for ch = 1:2
    fit_cutoff{m,ch} = zeros(size(merged_itraces{m,ch},1),1);
    def_fc = iEndval_sorted{m,ch}(1,2);
    for j = 1:size(merged_itraces{m,ch},1) %cycle through spots, in iEndval ascending order
        
        act_spotnum = iEndval_sorted{m,ch}(j,1);
        x_0 = round(merged_itraces{m,ch}{act_spotnum}(1,2));
        y_0 = round(merged_itraces{m,ch}{act_spotnum}(1,3));
        x_0_end = round(x_0 + ch1{m}.drift(end,1));
        y_0_end = round(y_0 + ch2{m}.drift(end,2));
        
        subplot('Position', [0.05,0.55,0.65,0.4])
        hold off
        plot (1:size(iEndval_sorted{m,ch},1),iEndval_sorted{m,ch}(:,2), [channel{ch} '.'], 'MarkerSize', 3);
        hold on
        plot (j,iEndval_sorted{m,ch}(j,2), 'ko', 'MarkerSize', 4);
        plot (1:length(iEndval_sorted{m,ch}),ones(1,length(iEndval_sorted{m,ch}))*1.5*avg_iEndval_tenth(m,ch), '-b');
        plot (1:length(iEndval_sorted{m,ch}),ones(1,length(iEndval_sorted{m,ch}))*1.5*avg_iEndval(m,ch), '-k');
        title('Average end value distribution (ascending order)')
        
        subplot('Position', [0.05,0.05,0.65,0.4])
        hold off
        plot(merged_itraces{m,ch}{act_spotnum}(:,1), ...
            merged_itraces{m,ch}{act_spotnum}(:,4),...
            ['-' channel{ch}(1)], 'LineWidth', 0.5)
        hold on
        plot(merged_itraces{m,ch}{act_spotnum}(:,1), ...
            merged_itraces{m,4}{act_spotnum,ch},...
            '-k', 'LineWidth', 0.25)
        title('Intensity trace')
        plot(merged_itraces{m,ch}{act_spotnum}(:,1), ones(1,size(merged_itraces{m,ch}{act_spotnum},1)).*def_fc, ...
            'color', [1 1 1].*0.5, 'LineStyle', ':')
        
        
        subplot('Position', [0.75,0.55,0.25,0.4])
        hold off
        plot_subframe(avg_img{m,ch}, x_0, y_0, 6)
        hold on
        ellipse(r_integrate, r_integrate, 0, x_0, y_0, channel{ch})
        title('Averaged over first 100 frames');
        axis square
        set(gca, 'YDir', 'normal')
        
        subplot('Position', [0.75,0.05,0.25,0.4])
        hold off
        plot_subframe(avg_img{m,ch+2}, x_0_end, y_0_end, 6)
        hold on
        ellipse(r_integrate, r_integrate, 0, x_0_end, y_0_end, channel{ch})
        title('Averaged over last 100 frames');
        axis square
        set(gca, 'YDir', 'normal')
        
        tmp = inputdlg(['Cutoff intensity for movie #' num2str(m) ', ' channel{ch} ' channel, spot #' ...
            num2str(act_spotnum)], 'Cutoff intensity', 1, {num2str(def_fc)});

        fit_cutoff{m,ch}(act_spotnum) = str2double(tmp);
        def_fc = fit_cutoff{m,ch}(j);
    end
    end
end
%%
%pos_in_frame: cell of arrays that for each frame gives starting fit
%coordinates for all spots in respective channel. If both parameters
%return zero, spot is not fitted in that frame.

pos_in_frame = cell(N_movie,2);
for m = 1:N_movie
    % channel 1
    pos_in_frame{m,1} = cell(size(ch1{m}.frames,2),1);
    for j = 1:size(pos_in_frame{m,1},1)
        pos_in_frame{m,1}{j} = zeros(size(merged_itraces{m,1},1),2);
        for s=1:size(pos_in_frame{m,1}{j},1)
        pos_in_frame{m,1}{j}(s,1:2) = (merged_itraces{m,4}{s,1}(j)>=fit_cutoff{m,1}(s))*(merged_itraces{m,1}{s}(j,2:3)); % x_0, y_0 remain zero if intensity is below threshold
        end
    end

    % channel 2
    pos_in_frame{m,2} = cell(size(ch2{m}.frames,2),1);
    for j = 1:size(pos_in_frame{m,2},1)
        pos_in_frame{m,2}{j} = zeros(size(merged_itraces{m,2},1),2);
        for s=1:size(pos_in_frame{m,2}{j},1)
        pos_in_frame{m,2}{j}(s,1:2) = (merged_itraces{m,4}{s,2}(j)>=fit_cutoff{m,2}(s))*(merged_itraces{m,2}{s}(j,2:3)); % x_0, y_0 remain zero if intensity is below threshold
        end
    end
end
   
%% FIT GAUSSIAN TO DATA % MAP 1 On 2 and 2 ON 1 and store
fit_result = cell(N_movie,1);

for m=1:N_movie %loop through movies
    fit_result{m} = cell (size(merged_itraces{m,1},1),2);
    for c=1:2
    for i=1:size(fit_result{m},1)
    fit_result{m}{i,c} = zeros(size(merged_itraces{m,c}{i},1),3);
    for j=1:size(merged_itraces{m,c}{i},1)
    fit_result{m}{i,c}(j,1:3) = merged_itraces{m,c}{i}(j,1:3)*...
        (merged_itraces{m,c}{i}(j,4)>fit_cutoff{m,c}(i));
    end
    
    fit_result{m}{i,c} = [fit_result{m}{i,c} zeros(size(fit_result{m}{i,c},1),9)];
    end
    end
    all_fit_result1 = cell(size(ch1{m}.frames,2),1);
    for i=1:size(all_fit_result1,1)
        all_fit_result1{i} = zeros(size(pos_in_frame{m,1}{i},1),12);
    end
    all_fit_result2 = cell(size(ch2{m}.frames,2),1);
    for i=1:size(all_fit_result2,1)
        all_fit_result2{i} = zeros(size(pos_in_frame{m,2}{i},1),12);
    end
end

%% par_fit
mycluster=parcluster('SharedCluster');
fit_job = batch(mycluster, @par_fit_v1, 1, { ch1, ch2, pos_in_frame, all_fit_result1, all_fit_result2, channel, N_movie, path_out} ...
    ,'CaptureDiary',true, 'CurrentDirectory', '.', 'Pool', 63 ...
    ,'AdditionalPaths', {[matlab_dir filesep 'TOOLBOX_GENERAL'], [matlab_dir filesep 'TOOLBOX_MOVIE'], [matlab_dir filesep 'FM_applications']});

%% retrieve data... load par_fit_result.mat file into fit_result.
cd(path_out)
load('par_fit_result')

%%
for m =1:size(ch1,1)
for s = 1:size(fit_result{m},1)
fit_result{m}{s,1}(:,6:12) = par_fit_result{m}{s,1}(:,6:12);
fit_result{m}{s,2}(:,6:12) = par_fit_result{m}{s,2}(:,6:12);
end
end


%% map fitted coordinates
if mapping
    for m = 1:size(fit_result,1)
    for i = 1:size(fit_result{m},1)
    for j = 1:size(fit_result{m}{i,1},1)
        if fit_result{m}{i,1}(j,1)~=0
        fit_result{m}{i,1}(j,4:5) = transformPointsInverse(tform_2TO1, fit_result{m}{i,1}(j,6:7));  %%this takes coords in ch1 and transfroms them to coords in ch2
        end
    end
    for j = 1:size(fit_result{m}{i,2},1)
        if fit_result{m}{i,2}(j,1)~=0
        fit_result{m}{i,2}(j,4:5) = transformPointsInverse(tform_1TO2, fit_result{m}{i,2}(j,6:7));  %%this takes coords in ch2 and transfroms them to coords in ch1
        end
    end
    end
    end
end

% save data
cd(path_out)
save -v7.3 'all_data.mat'

disp('Done')

