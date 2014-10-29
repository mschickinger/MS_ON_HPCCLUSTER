% What's new in version 10a?

% Preparative steps for batch job assignment: gaussian fitting procedure is
% adapted for less data communication overhead -> supposedly more efficient
% working order is changed to movie -> frame -> spots
% mapping of fitted positions is drawn out of parallel loop

% Also includes drift (x,y) correction routine (optional) and new
% intensity tracer.
% -> ALSO INCLUDE IN v10b!!!

%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs.m')


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
        pos1on2{i} = transformPointsInverse(tform_2TO1, ch1_pos);  %%this takes coords in ch1 and transforms them to coords in ch2
        pos2on1{i} = transformPointsInverse(tform_1TO2, ch2_pos);  %%this takes coords in ch2 and transforms them to coords in ch1
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
    
% find spot pairs
pairs_for_drift_cor = cell(N_movie,2);

for i=1:N_movie
    
    [h_min_end, pos_ch1_end] = ch1{i}.get_threshold(r_find, avg_img{i,3});
    [h_min_end, pos_ch2_end] = ch2{i}.get_threshold(r_find, avg_img{i,4});
    
    pos_ch1_start = peaks_raw(peaks_raw(:,5)==i,1:2);
    pos_ch2_start = peaks_raw(peaks_raw(:,5)==i,3:4);
    
    % map peaks
    trace_map_drift_ch1 = map_traces(pos_ch1_start(:,1:2), pos_ch1_end(:,1:2), pos_ch1_end(:,1:2), r_find*2)+1;
    trace_map_drift_ch2 = map_traces(pos_ch2_start(:,1:2), pos_ch2_end(:,1:2), pos_ch2_end(:,1:2), r_find*2)+1;

    tmp_ch1 = zeros(size(trace_map_drift_ch1,1),4);    
    tmp_ch2 = zeros(size(trace_map_drift_ch2,1),4);
    
    % combine pairs
    for j=1:size(trace_map_drift_ch1,1)
        tmp_ch1(j,:) = [pos_ch1_start(trace_map_drift_ch1(j,1),1:2)+1 pos_ch1_end(trace_map_drift_ch1(j,2),1:2)+1];
    end
    for j=1:size(trace_map_drift_ch2,1)
        tmp_ch2(j,:) = [pos_ch2_start(trace_map_drift_ch2(j,1),1:2)+1 pos_ch2_end(trace_map_drift_ch2(j,2),1:2)+1];
    end   
    pairs_for_drift_cor{i,1} = [tmp_ch1 zeros(size(tmp_ch1,1),4)];
    pairs_for_drift_cor{i,2} = [tmp_ch2 zeros(size(tmp_ch2,1),4)];
end

% fit psf to spots

s_x = 2.5;
s_y = 2.5;
w_fit = 5;
h = waitbar(0,'Fitting spots for drift correction.. please wait');
for i=1:N_movie
    for ch = 1:2
        for s = 1:size(pairs_for_drift_cor{i,ch},1)
            x1 = round(pairs_for_drift_cor{i,ch}(s,1));
            y1 = round(pairs_for_drift_cor{i,ch}(s,2));
            x2 = round(pairs_for_drift_cor{i,ch}(s,3));
            y2 = round(pairs_for_drift_cor{i,ch}(s,4));
            [c1] = fit_gauss2d_mainaxis_bg(x1, y1, s_x, w_fit, avg_img{i,ch});
            [c2] = fit_gauss2d_mainaxis_bg(x2, y2, s_x, w_fit, avg_img{i,ch+2});
            pairs_for_drift_cor{i,ch}(s,5:6) = c1(1:2);
            pairs_for_drift_cor{i,ch}(s,7:8) = c2(1:2);
        end
    end
waitbar(i/N_movie , h, ['Fitting spots from movie... ' num2str(i) ' of ' num2str(N_movie) ' done'])
end
close(h);

for i = 1:N_movie
    ch1{i}.drift(1) = mean(pairs_for_drift_cor{i,1}(:,7)-pairs_for_drift_cor{i,1}(:,5))/(ch1{i}.mov_length-N_frames);
    ch1{i}.drift(2) = mean(pairs_for_drift_cor{i,1}(:,8)-pairs_for_drift_cor{i,1}(:,6))/(ch1{i}.mov_length-N_frames);
    ch2{i}.drift(1) = mean(pairs_for_drift_cor{i,2}(:,7)-pairs_for_drift_cor{i,2}(:,5))/(ch2{i}.mov_length-N_frames);
    ch2{i}.drift(2) = mean(pairs_for_drift_cor{i,2}(:,8)-pairs_for_drift_cor{i,2}(:,6))/(ch2{i}.mov_length-N_frames);
end
end %if-loop

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
    ch1_itraces_full = ch1{i}.int_spots_in_frames(ch1{i}.frames, peaks(peaks(:,5)==i,1:2), r_integrate);
    display(['Tracing ' channel{1} ' channel in movie #' num2str(i) ' done'])
    toc %
    ch2_itraces_full = ch2{i}.int_spots_in_frames(ch2{i}.frames, peaks(peaks(:,5)==i,3:4), r_integrate);
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

fc = figure('PaperPositionMode', 'manual', 'PaperUnits', 'centimeters',...
'PaperPosition', [0 0 25 9]);

for m = 1:N_movie
    for ch = 1:2
    fit_cutoff{m,ch} = zeros(size(merged_itraces{m,ch},1),1);
    def_fc = iEndval_sorted{m,ch}(1,2);
    for j = 1:size(merged_itraces{m,ch},1) %cycle through spots, in iEndval ascending order
        
        act_spotnum = iEndval_sorted{m,ch}(j,1);
        x_0 = merged_itraces{m,ch}{act_spotnum}(1,2);
        y_0 = merged_itraces{m,ch}{act_spotnum}(1,3);
        x_0_end = round(x_0 + (ch==1)*ch1{m}.mov_length*ch1{m}.drift(1) + (ch==2)*ch2{m}.mov_length*ch2{m}.drift(1));
        y_0_end = round(y_0 + (ch==1)*ch1{m}.mov_length*ch1{m}.drift(2) + (ch==2)*ch2{m}.mov_length*ch2{m}.drift(2));
        
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
    
    parfor_progress(length(ch1{m}.frames));
end

for m=1:N_movie %loop through movies
parfor p=1:length(ch1{m}.frames)
    display(['Fitting ' channel{1} ' spots in frame #' num2str(ch1{m}.frames(p)) ' in movie ' num2str(m) ' / ' num2str(N_movie)])
    tmp1 = ch1{m}.fit_sym_psfs_to_frame(ch1{m}.frames(p), pos_in_frame{m,1}{p}, 2); %fit all spots that occur in this frame, sigma = 2

    for j=1:size(tmp1,1)
        all_fit_result1{p}(j,1) = tmp1(j,1);
        all_fit_result1{p}(j,2) = tmp1(j,2);
        all_fit_result1{p}(j,3) = tmp1(j,3);
        all_fit_result1{p}(j,4) = tmp1(j,4);
        all_fit_result1{p}(j,5) = tmp1(j,5);
        all_fit_result1{p}(j,6) = tmp1(j,6);
    end
    parfor_progress;
end
msgbox(['Finished fitting all ' channel{1} ' spots from movie ' num2str(m)])
%write fit loop data to fit_result cell
for j=1:size(merged_itraces{m,1},1)
    for i=1:size(all_fit_result1,1)
        fit_result{m}{j,1}(i,6:8) = all_fit_result1{i}(j,1:3);
        fit_result{m}{j,1}(i,10:12) = all_fit_result1{i}(j,4:6);
    end
    fit_result{m}{j,1}(:,9) = fit_result{m}{j,1}(:,8);
end

parfor_progress(0);
parfor_progress(length(ch2{m}.frames));

parfor p = 1:length(ch2{m}.frames);
    display(['Fitting ' channel{2} ' spots in frame #' num2str(ch2{m}.frames(p)) ' in movie ' num2str(m) ' / ' num2str(N_movie)])
    tmp2 = ch2{m}.fit_sym_psfs_to_frame(ch2{m}.frames(p), pos_in_frame{m,2}{p}, 2); %fit spot in each frame, sigma = 2

    for j=1:size(tmp2,1);
        all_fit_result2{p}(j,1) = tmp2(j,1);
        all_fit_result2{p}(j,2) = tmp2(j,2);
        all_fit_result2{p}(j,3) = tmp2(j,3);
        all_fit_result2{p}(j,4) = tmp2(j,4);
        all_fit_result2{p}(j,5) = tmp2(j,5);
        all_fit_result2{p}(j,6) = tmp2(j,6);
    end
    parfor_progress;
end
msgbox(['Finished fitting all spots from movie ' num2str(m)])
%write fit loop data to fit_result cell
for j=1:size(merged_itraces{m,2},1)
    for i=1:size(all_fit_result2,1)
        fit_result{m}{j,2}(i,6:8) = all_fit_result2{i}(j,1:3);
        fit_result{m}{j,2}(i,10:12) = all_fit_result2{i}(j,4:6);
    end
    fit_result{m}{j,2}(:,9) = fit_result{m}{j,2}(:,8);
end

parfor_progress(0);

% map fitted coordinates on other channel
if mapping
    for i = 1:size(merged_itraces{m,1},1)
    for j = 1:size(fit_result{m}{i,1},1)
    fit_result{m}{i,1}(j,4) = transformPointsInverse(tform_2TO1, fit_result{m}{i,1}(j,6));  %%this takes coords in ch1 and transfroms them to coords in ch2
    fit_result{m}{i,1}(j,5) = transformPointsInverse(tform_2TO1, fit_result{m}{i,1}(j,7));  %%this takes coords in ch1 and transfroms them to coords in ch2
    end
    for j = 1:size(fit_result{m}{i,2},1)
    fit_result{m}{i,2}(j,4) = transformPointsInverse(tform_1TO2, fit_result{m}{i,2}(j,6));  %%this takes coords in ch2 and transfroms them to coords in ch1
    fit_result{m}{i,2}(j,5) = transformPointsInverse(tform_1TO2, fit_result{m}{i,2}(j,7));  %%this takes coords in ch2 and transfroms them to coords in ch1
    end
    end
end

msgbox(['Finished fitting all spots from movie ' num2str(m)])

% save data
save([path_out filesep 'data_until_movie' num2str(m) '.mat']);
end


%% Plot data
delta = 5;

for m=1:N_movie
for i=1:length(fit_result{m})

% traces for delta_r and sigma over frame number
cf = figure('Visible', 'Off', 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters',...
'PaperPosition', [1 1 32 18]);

    x0_ch1 = merged_itraces{m,1}{i}(1,2);
    y0_ch1 = merged_itraces{m,1}{i}(1,3);
    xy_ch1 = fit_result{m}{i,1}(fit_result{m}{i,1}(:,1)>0,6:7);
    xy_mean_ch1 = [mean(xy_ch1(abs(xy_ch1(:,1)-x0_ch1)<=5,1))...
        mean(xy_ch1(abs(xy_ch1(:,2)-y0_ch1)<=5,2))];
    d_ch1 = sqrt((xy_ch1(:,1)-xy_mean_ch1(1)).^2 + (xy_ch1(:,2)-xy_mean_ch1(2)).^2);
    
    sigma_ch1 = sqrt(fit_result{m}{i,1}(fit_result{m}{i,1}(:,1)>0,8).^2  + fit_result{m}{i,1}(fit_result{m}{i,1}(:,1)>0,9).^2);
    frame_ch1 = fit_result{m}{i,1}(fit_result{m}{i,1}(:,1)>0,1);
    chisquared_ch1 = fit_result{m}{i,1}(fit_result{m}{i,1}(:,1)>0,12);
    meanchisquared_ch1=mean(chisquared_ch1);
    
    x0_ch2 = merged_itraces{m,2}{i}(1,2);
    y0_ch2 = merged_itraces{m,2}{i}(1,3);    
    xy_ch2 = fit_result{m}{i,2}(fit_result{m}{i,2}(:,1)>0,6:7);
    xy_mean_ch2 = [mean(xy_ch2(abs(xy_ch2(:,1)-x0_ch2)<=5,1))...
        mean(xy_ch2(abs(xy_ch2(:,2)-y0_ch2)<=5,2))];
    d_ch2 = sqrt((xy_ch2(:,1)-xy_mean_ch2(1)).^2 + (xy_ch2(:,2)-xy_mean_ch2(2)).^2);
    
    sigma_ch2 = sqrt(fit_result{m}{i,2}(fit_result{m}{i,2}(:,1)>0,8).^2  + fit_result{m}{i,2}(fit_result{m}{i,2}(:,1)>0,9).^2);
    frame_ch2 = fit_result{m}{i,2}(fit_result{m}{i,2}(:,1)>0,1);
    chisquared_ch2 = fit_result{m}{i,2}(fit_result{m}{i,2}(:,1)>0,12);
    meanchisquared_ch2 = mean(chisquared_ch2);
    
    center = xy_mean_ch1;
    if size(xy_ch1,1)>0 && size(xy_ch2,1)>0
    center = ((size(xy_ch1,1)>0)*xy_mean_ch1+(size(xy_ch2,1)>0)*xy_mean_ch2)./2;    
    else if size(xy_ch1,1)==0 && size(xy_ch2,1)>0
            center = xy_mean_ch2;
        else
            center = [delta delta];
        end
    end
if length(frame_ch1)>0
    subplot(3, 2, 1)
    plot(frame_ch1, d_ch1, ['.-' channel{1}(1)], 'Markersize', 5, 'LineWidth', 0.5)
    title(['Distance from average position. Spot: ' num2str(i) ' / movie #' num2str(m)])
    xlabel('Frame')
    ylabel('$\sqrt{{\Delta x}^2 +{\Delta y}^2}$', 'Interpreter', 'latex')
    xlim([1 frame_ch1(end)])
    ylim([0 5])

    subplot(3, 2, 3)
    plot(frame_ch1, sigma_ch1, ['.-' channel{1}(1)], 'Markersize', 5, 'LineWidth', 0.5)
    title(['Width of peak. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel('$\sqrt{{\sigma_x}^2+{\sigma_y}^2}$', 'Interpreter', 'latex')
    xlim([1 frame_ch1(end)])
    ylim([0 5])
    
    subplot(3, 2, 5)
    %bar(frame_green, chi2green)
    plot(frame_ch1, chisquared_ch1, ['.-' channel{1}(1)], 'Markersize', 5, 'LineWidth', 0.5)
    title(['Chi^2 of ' rgb{colors(1)} ' fit. Mean value: ' sprintf('%.2s',meanchisquared_ch1) '. Spot: ' num2str(i)])    
    xlim([1 frame_ch1(end)])
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(1)} '}$'], 'Interpreter', 'latex')
else
    subplot(3, 2, 1)
    cla
    subplot(3, 2, 3)
    cla
    subplot(3, 2, 5)
    cla
end

if length(frame_ch2)>0
    subplot(3, 2, 2)
    plot(frame_ch2, d_ch2, ['.-' channel{2}(1)], 'Markersize', 5, 'LineWidth', 0.5)
    xlabel('Frame')
    ylabel('$\sqrt{{\Delta x}^2 +{\Delta y}^2}$', 'Interpreter', 'latex')
    xlim([1 frame_ch2(end)])
    ylim([0 5])
   
    subplot(3, 2, 4)
    plot(frame_ch2, sigma_ch2, ['.-' channel{2}(1)], 'Markersize', 5, 'LineWidth', 0.5)
    xlabel('Frame')
    ylabel('$\sqrt{{\sigma_x}^2+{\sigma_y}^2}$', 'Interpreter', 'latex')
    xlim([1 frame_ch2(end)])
    ylim([0 5])
      
    subplot(3, 2, 6)
    plot(frame_ch2, chisquared_ch2, ['.-' channel{2}(1)], 'Markersize', 5, 'LineWidth', 0.5)
    title(['Chi^2 of ' rgb{colors(2)} ' fit. Mean value: ' sprintf('%.2s',meanchisquared_ch2) '. Spot: ' num2str(i)])    
    xlim([1 frame_ch2(end)])
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(2)} '}$'], 'Interpreter', 'latex')
else
    subplot(3, 2, 2)
    cla
    subplot(3, 2, 4)
    cla
    subplot(3, 2, 6)
    cla
end
    
    print(cf, '-dpng', [path_out filesep 'trace_m' num2str(m) '_s' num2str(i) '.png' ])
    
close all
    
% traces for position (dots unconnected)
cf = figure('Visible', 'Off');    
hold off
    
    plot(xy_ch1(:,1), xy_ch1(:,2), [channel{1}(1) 'o']), hold on
    plot(xy_ch2(:,1), xy_ch2(:,2), [channel{2}(1) 'x']), hold off
    axis([floor(center(1)-delta)  ceil(center(1)+delta) floor(center(2)-delta)  ceil(center(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([channel{1} ' and ' channel{1} ' fitted positions (unmapped). Spot: ' num2str(i) ' / movie #' num2str(m)]);
    hold off;
    axis([floor(center(1)-delta)  ceil(center(1)+delta) floor(center(2)-delta)  ceil(center(2)+delta)]);
    print(cf, '-dpng', [path_out filesep 'position_m' num2str(m) '_s' num2str(i) '.png' ])
    
% traces for position (dots connected)
  
hold off
    
    xy_ch1 = [xy_ch1(d_ch1<=delta,1:2)];
    xy_ch2 = [xy_ch2(d_ch2<=delta,1:2)];
    plot(xy_ch1(:,1), xy_ch1(:,2), [channel{1}(1) '-o']), hold on
    plot(xy_ch2(:,1), xy_ch2(:,2), [channel{2}(1) '-x']), hold off
    axis([floor(center(1)-delta)  ceil(center(1)+delta) floor(center(2)-delta)  ceil(center(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([channel{1} ' and ' channel{1} ' fitted positions (unmapped). Spot: ' num2str(i) ' / movie #' num2str(m)]);
    hold off;
    axis([floor(center(1)-delta)  ceil(center(1)+delta) floor(center(2)-delta)  ceil(center(2)+delta)]);
    print(cf, '-dpng', [path_out filesep 'position2_m' num2str(m) '_s' num2str(i) '.png' ])

close all
end
end

close all

% traces for position (dots unconnected, red dots mapped rONg)

if mapping
delta = 5;
cf = figure(1);
for i=1:size(merged_traces,1)

    xy_ch1 = fit_result{m}{i,1}(fit_result{m}{i,1}(:,1)>0,6:7);
    xy_mean_ch1 = [mean(xy_ch1((xy_ch1(:,1)-merged_itraces{m,1}{i}(1,2))<=5,1))...
    mean(xy_ch1((xy_ch1(:,2)-merged_itraces{m,1}{i}(1,3))<=5,2))];
    
    xy_2on1 = fit_result{m}{i,2}(fit_result{m}{i,2}(:,1)>0,4:5);
    
    plot(xy_ch1(:,1), xy_ch1(:,2), [channel{1}(1) 'o']), hold on
    plot(xy_2on1(:,1), xy_2on1(:,2), [channel{2}(1) 'x']), hold off
    axis([floor(xy_mean_ch1(1)-delta)  ceil(xy_mean_ch1(1)+delta) floor(xy_mean_ch1(2)-delta)  ceil(xy_mean_ch1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (rONg mapped). Spot: ' num2str(i) ' / movie #' num2str(m)]);
    hold off;
    axis([floor(xy_mean_ch1(1)-delta)  ceil(xy_mean_ch1(1)+delta) floor(xy_mean_ch1(2)-delta)  ceil(xy_mean_ch1(2)+delta)]);
    print(cf, '-dpng', [path_out filesep 'position_mapped_rONg_m' num2str(m) '_s' num2str(i) '.png' ])
end
end

close all


%% Plot positions of all found spots on channel 1 average image

cf=figure(1);
for i=1:N_movie
    imagesc(avg_img{i,1}),colormap gray, colorbar, axis image;
    hold on
    plot(all_positions{i,1}(:,1),all_positions{i,1}(:,2),[channel{1}(1) 'o']);
    plot(all_positions{i,2}(:,1),all_positions{i,2}(:,2),[channel{2}(1) 'x']);
    print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{1} '_' channel{2} '.eps'])
    hold off
end

% Plot positions of all found spots on channel 2 average image
cf=figure(1);
for i=1:N_movie
    imagesc(avg_img{i,2}),colormap gray, colorbar, axis image;
    hold on
    plot(all_positions{i,1}(:,1),all_positions{i,1}(:,2),[channel{1}(1) 'x']);
    plot(all_positions{i,2}(:,1),all_positions{i,2}(:,2),[channel{2}(1) 'o']);
    print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{2} '_' channel{1} '.eps'])
    hold off
end

% Plot positions taking into account bead mapping

if mapping
    
    cf=figure(1);
    for i=1:N_movie
        imagesc(avg_img{i,1}),colormap gray, colorbar, axis image;
        hold on
        plot(all_positions{i,1}(:,1),all_positions{i,1}(:,2),[channel{1} 'o']);
        plot(pos2on1{i}(:,1),pos2on1{i}(:,2),[channel{2} 'x']);

        print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{2} 'TO' channel{1} '.eps'])
        hold off
    end
    
    for i=1:N_movie
        imagesc(avg_img{i,2}),colormap gray, colorbar, axis image;
        hold on
        plot(pos1on2{i}(:,1),pos1on2{i}(:,2),[channel{1} 'x']);
        plot(all_positions{i,2}(:,1),all_positions{i,2}(:,2),[channel{2} 'o']);

        print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{1} 'TO' channel{2} '.eps'])
        hold off
    end
    
end

% End of program
disp('Done')
