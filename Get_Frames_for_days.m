
% load relevant data
animal = 'AP021';
verbose = true; % this turns on/off progress display in command line
recordings = plab.find_recordings(animal,[],'visual_conditioning_prob_right');


for i=1:length({recordings.day})
    % grab  the recording day and time for each day
    rec_day= recordings(i).day;
    rec_time= char(recordings(i).recording);

    if size(rec_time) >1 % check if there is more than 1 recording, if so then take the second trial
        rec_time= rec_time(2,:)
    end

    load_parts.mousecam= true; % to not load the widefield data as well
    ap.load_recording;
    timelite
    bonsai_workflow

    % Bonsai can change how it saves data, and loading scripts can specify how
    % that data is loaded and parsed. In this demo, ap.load_recording loads
    % data from Bonsai as 'trial_events':
    % trial_events
    % mousecam_fn;
    % mousecam_times;
    % 
    mousecam_vr = VideoReader(mousecam_fn);


    % read mousecam
    mousecam_video= read(mousecam_vr,1);
    % find the center indecies
    center_point_indecies= [round(size(mousecam_video,1)*2/3),round(size(mousecam_video,2)*0.56)];
    mouscam_video_right_quadrant= mousecam_video(1:center_point_indecies(1),center_point_indecies(2):size(mousecam_video,2));

    % Plot the video frame
    subplot(3,4, i);
    imagesc(mouscam_video_right_quadrant);
    colormap gray;

    % figure;imagesc(mouscam_video_right_quadrant)

    title(['Recording Day: ', rec_day], 'FontSize', 8);

end
% title(['Frame for Day ', rec_day, ' for Animal ', animal]);

% Get the current axes
h = gca;

% Get the children objects within the axes (assuming the heatmap is the only object)
children = get(h, 'Children');

% Check if children is an array (if there are multiple objects)
if numel(children) == 1
    % If there's only one child, access its CData and XData properties
    c = get(children, 'CData');
    x = get(children, 'XData');
else
    % If there are multiple children, find the heatmap object and access its CData and XData properties
    for i = 1:numel(children)
        if isprop(children(i), 'CData') && isprop(children(i), 'XData')
            c = get(children(i), 'CData');
            x = get(children(i), 'XData');
            break;
        end
    end
end

figure;plot(x,c(1,:),LineWidth=2)
hold on;
plot(x,c(5,:),LineWidth=2)
plot(x,c(9,:),LineWidth=2)
xline(0,'-')
xlabel('Time (sec)');
ylabel('Fraction Moving');
title('Change in movement across days')

h=gca;
[h.XAxis.FontSize, h.YAxis.FontSize] = deal(20);

%format the plots
h=gca;
h.Title.FontSize=15;


h1 = gca;
children1 = get(h, 'Children');


% Filter out Line objects
line_objects = children(strcmpi(get(children, 'Type'), 'line'));

% Now you can work with the line_objects, for example, get their XData and YData
for i = 1:numel(line_objects)
    x1_data = get(line_objects(i), 'XData');
    y1_data = get(line_objects(i), 'YData');

    % Do something with x_data and y_data
end

figure;plot(x1_data(1:end-3),log(y1_data(1:end-3)))





% Extract data from the first plot
h1 = gca;   % Get the current axes for the first plot
line_objects1 = findobj(h1, 'Type', 'line');  % Find all line objects in the first plot

if ~isempty(line_objects1)
    x1_data = line_objects1(2).XData;
    y1_data = line_objects1(2).YData;
else
    error('No line objects found in the first plot.');
end

figure;plot(x1_data,y1_data,LineWidth=2)
ylim([0,10]);
xlim([0,14]);
xlabel('Days');
ylabel('Median Reaction Time');
title('Median Reaction Time Change Across Days')



%format the plots
h=gca;
[h.XAxis.FontSize, h.YAxis.FontSize] = deal(20);
h.Title.FontSize=15;