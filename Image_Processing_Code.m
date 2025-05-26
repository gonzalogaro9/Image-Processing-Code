clc; clear; close all;

%% Mode Selection
option = questdlg('What do you want to process?', 'Input Selection', 'Single Image', 'Entire Folder', 'Cancel', 'Single Image');

if strcmp(option, 'Single Image')
    [file, folder] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif'}, 'Select an image');
    if file == 0, disp('Selection canceled.'); return; end
    images = {fullfile(folder, file)};
elseif strcmp(option, 'Entire Folder')
    folder = uigetdir('', 'Select a folder');
    if folder == 0, disp('Selection canceled.'); return; end
    files = dir(fullfile(folder, '*.png'));
    images = fullfile(folder, {files.name});
else
    disp('Operation canceled.'); return;
end

%% Scale
scale_options = {'50', '100', '200', '500', '1000'};
scale_selection = listdlg('PromptString', 'Select scale in microns:', 'SelectionMode', 'single', 'ListString', scale_options);
if isempty(scale_selection), disp('Selection canceled.'); return; end
scale_microns = str2double(scale_options{scale_selection});

%% Initialize Excel variables
excel_file = 'Particle_Analysis.xlsx';
if exist(excel_file, 'file')
    delete(excel_file); 
end

% Structure to store total data
total_data = struct();

%% Initialize total variables
total_ball_areas = []; total_fiber_areas = []; total_chip_areas = [];
total_fiber_lengths = []; total_ball_circularities = [];
total_chip_eccentricities = []; total_particle_count = 0;

%% Colors for plots
colors = {'b', 'r', 'g', 'm', 'c', 'y', [0.25 0.5 0.9]};

%% Pre-processing
for i = 1:length(images)
    img_path = images{i};
    [~, img_name, ext] = fileparts(img_path);
    full_name = [img_name, ext];
    
    %% Show original image
    img = imread(img_path);
    figure('Name','Original'); imshow(img); title('Original Image', 'Interpreter','none');
    text(10,size(img, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');

    %% Grayscale
    if size(img,3)==3
        img_gray = rgb2gray(img);
    else
        img_gray = img;
    end
    figure('Name','Grises'); imshow(img_gray); title('Gray Scale', 'Interpreter','none');
    text(10,size(img_gray, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');

    %% Median filter
    img_median = medfilt2(img_gray, [5 5]);
    figure('Name','Median'); imshow(img_median); title('Median Filter', 'Interpreter','none');
    text(10,size(img_median, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');

    %% Binarization for Scale
    img_binary = img_gray > 250;
    img_binary = bwareaopen(img_binary, 100);
    figure('Name','Binarization'); imshow(img_binary); title('Binarization for the Scale', 'Interpreter','none');
    text(10,size(img_binary, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');
   
    %% Scale Detection
    scale_stats = regionprops(img_binary, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid');
    scale_length_px = NaN;

    for s = 1:length(scale_stats)
        orientation = abs(scale_stats(s).Orientation);
        cy = scale_stats(s).Centroid(2);
        
        if orientation < 10 && cy > size(img_binary,1)*0.8 && scale_stats(s).MajorAxisLength > 30
            scale_length_px = scale_stats(s).MajorAxisLength;
            
            % Coordinates to highlight the scale
            bbox = scale_stats(s).BoundingBox;
            cx = scale_stats(s).Centroid(1);
    
            % Show image with highlighted scale
            figure, imshow(img_binary); hold on;
            rectangle('Position', bbox, 'EdgeColor', 'cyan', 'LineWidth', 2);
            text(cx, cy+30, sprintf('%.2f px', scale_length_px),'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold','HorizontalAlignment', 'center');
            hold off;
    
            break; 
        end
    end
    
    if isnan(scale_length_px)
        warning('No scale detected.');
        return;
    end
    
    conversion_factor = scale_microns / scale_length_px;
   
    %% Remove Scale
    img_no_scale = img_gray;
    img_no_scale(img_binary) = median(img_gray(:));
    figure('Name','Non-Scale'); imshow(img_no_scale); title('Scale Eliminated', 'Interpreter','none');
    text(10,size(img_no_scale, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');

    %% Gaussian Filter
    filtered = imgaussfilt(img_no_scale, 1.5);
    figure('Name','Gaussian'); imshow(filtered); title('Gaussian Filter', 'Interpreter','none');
    text(10,size(filtered, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');

   %% Canny + Dilation + Fill
    edges = edge(filtered, 'Canny', [0.2 0.4]);
    dilated = imdilate(edges, strel('disk', 5));
    opened = imopen(dilated, strel('disk', 5));
    closed = imclose(opened, strel('disk', 10));
    closed = imfill(closed, 'holes');
    closed = bwareaopen(closed, 500);
    [B, L] = bwboundaries(imresize(closed,2,"nearest"));
    figure('Name','Particles'); imshow(closed); title('Particles Detected', 'Interpreter','none');
    text(10,size(closed, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');

    %% Properties
    stats = regionprops(L, 'Area','Perimeter','MajorAxisLength','MinorAxisLength','Eccentricity','Centroid');
    if isempty(stats), continue; end
    A_px2 = [stats.Area]; 
    P_px = [stats.Perimeter];
    major_length_px = [stats.MajorAxisLength]; 
    minor_length_px = [stats.MinorAxisLength]; 
    Ecc = [stats.Eccentricity];
    Centroid = reshape([stats.Centroid], 2, [])';

    A_um2 = A_px2 * (conversion_factor)^2;
    P_um = P_px * conversion_factor; 
    major_length_um = major_length_px * conversion_factor;
    minor_length_um = minor_length_px * conversion_factor;
    circ = (4*pi*A_um2) ./ (P_um.^2); 
    aspect_ratio = major_length_px ./ minor_length_px;
    A_um2 = A_um2/4;
    P_um = P_um/2;
    major_length_um = major_length_um/2;
    minor_length_um = minor_length_um/2;
    Centroid = Centroid/2;
    
    %% Classification
    balls = find(circ >= 0.80);
    fibers = find(aspect_ratio >= 1.8 & Ecc >= 0.80 & A_um2 >= 19000 & major_length_um > 50 & circ < 0.3);
    chips = find(~ismember(1:length(stats), union(balls, fibers)) & (circ < 0.6) & (Ecc >= 0.4 & Ecc < 0.8) & (P_um ./ sqrt(A_um2) > 4.5));

    %% Prepare Data for Excel
    num_particles = length(stats);
    image_data = cell(num_particles + 1, 10);
    
    image_data(1,:) = {'ID', 'Type', 'Area (µm²)', 'Circularity/Length/Eccentricity','Major Length (µm)', 'X Position', 'Y Position', 'Image Name', 'Scale (µm)', 'Scale (px)'};
    
    for j = 1:num_particles
        if ismember(j, balls)
            type = 'Ball';
            measurement = circ(j);
        elseif ismember(j, fibers)
            type = 'Fiber';
            measurement = major_length_um(j);
        else
            type = 'Chip';
            measurement = Ecc(j);
        end
        
        image_data{j+1,1} = j;
        image_data{j+1,2} = type;
        image_data{j+1,3} = A_um2(j);
        image_data{j+1,4} = measurement;
        image_data{j+1,5} = major_length_um(j);
        image_data{j+1,6} = Centroid(j,1);
        image_data{j+1,7} = Centroid(j,2);
        image_data{j+1,8} = full_name;  
        image_data{j+1,9} = scale_microns;
        image_data{j+1,10} = scale_length_px;
    end
    
    sheet_name = img_name;
    if length(sheet_name) > 31 
        sheet_name = sheet_name(1:31);
    end
    
    xlswrite(excel_file, image_data, sheet_name);
    
    % Save image summary
    total_data(i).image_name = img_name;
    total_data(i).num_balls = length(balls);
    total_data(i).num_fibers = length(fibers);
    total_data(i).num_chips = length(chips);
    total_data(i).ball_areas = A_um2(balls);
    total_data(i).fiber_areas = A_um2(fibers);
    total_data(i).chip_areas = A_um2(chips);
    total_data(i).fiber_lengths = major_length_um(fibers);
    total_data(i).ball_circularities = circ(balls);
    total_data(i).chip_eccentricities = Ecc(chips);
    
    %% Show Classification with Interactive Tooltips
    figure('Name','Classification');
    imshow(filtered); hold on;
    [B,~] = bwboundaries(closed,'noholes');
    
    h_balls = []; h_fibers = []; h_chips = [];
    particle_data = cell(length(stats), 1);
    h_particles = gobjects(length(stats), 1);
    ab=[];
    af=[];
    av=[];
    cb=[];
    lf=[];
    ev=[];
    
    for j = 1:length(stats)
        boundary = B{j};
        if ismember(j, balls)
            h = plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 1.5);
            h_balls = [h_balls; h];
            type = 'Ball';
            measurement = sprintf('Circularity: %.2f', circ(j));
            ab(end+1)=A_um2(j);
            cb(end+1)=circ(j);
        elseif ismember(j, fibers)
            h = plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
            h_fibers = [h_fibers; h];
            type = 'Fiber';
            measurement = sprintf('Length: %.2f µm', major_length_um(j));
            af(end+1)=A_um2(j);
            lf(end+1)=major_length_um(j);
        else
            h = plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1.5);
            h_chips = [h_chips; h];
            type = 'Chip';
            measurement = sprintf('Eccentricity: %.2f', Ecc(j));
            av(end+1)=A_um2(j);
            ev(end+1)=Ecc(j);
        end
        h_particles(j) = h;
        abb=mean(ab);
        aff=mean(af);
        avv=mean(av);
        cbb=mean(cb);
        lff=mean(lf);
        evv=mean(ev);
        
        particle_data{j} = struct('ID', j, 'Type', type, 'Area', sprintf('%.2f µm²', A_um2(j)), 'Measurement', measurement, 'Position', sprintf('(%.1f, %.1f)', Centroid(j,1), Centroid(j,2)));
    end
    
    title('Classification (hover over particles)', 'Interpreter','none');
    text(10,size(filtered, 1) - 10,sprintf('Path: %s | Scale: %d µm', img_path, scale_microns),'Color','w','FontSize',8,'Interpreter','none','BackgroundColor','k');
    
    set([h_balls; h_fibers; h_chips], 'HitTest', 'on', 'PickableParts', 'visible');
    set(gcf, 'WindowButtonMotionFcn', @(src,evt) showTooltip(src, evt, h_particles, particle_data));
    
    hold off;

    %% Individual plots
    tit = @(txt) sprintf('%s\n%s | Scale: %d µm', txt, img_path, scale_microns);

    figure('Name','1 - Distribution');
    bar(categorical({'Balls','Fibers','Chips'}),100 * [length(balls), length(fibers), length(chips)] / length(stats),'FaceColor', colors{1});
    title(tit('Percentage distribution'), 'Interpreter','none'); ylabel('%');

    figure('Name','2 - Fiber area');
    if ~isempty(fibers)
        histogram(A_um2(fibers),'FaceColor',colors{2});
        title(tit('Fiber area'), 'Interpreter','none');
        xlabel('Area (µm²)'); ylabel('Count');
    end

    figure('Name','3 - Chip area');
    if ~isempty(chips)
        histogram(A_um2(chips),'FaceColor',colors{3});
        title(tit('Chip area'), 'Interpreter','none');
        xlabel('Area (µm²)'); ylabel('Count');
    end

    figure('Name','4 - Ball area');
    if ~isempty(balls)
        histogram(A_um2(balls),'FaceColor',colors{4});
        title(tit('Ball area'), 'Interpreter','none');
        xlabel('Area (µm²)'); ylabel('Count');
    end

    figure('Name','5 - Fiber length');
    if ~isempty(fibers)
        histogram(major_length_um(fibers),'FaceColor',colors{5});
        title(tit('Fiber length'), 'Interpreter','none');
        xlabel('Length (µm)'); ylabel('Count');
    end

    figure('Name','6 - Ball circularity');
    if ~isempty(balls)
        histogram(circ(balls),'FaceColor',colors{6});
        title(tit('Ball circularity'), 'Interpreter','none');
        xlabel('Circularity'); ylabel('Count');
    end

    figure('Name','7 - Chip eccentricity');
    if ~isempty(chips)
        histogram(Ecc(chips),'FaceColor',colors{7});
        title(tit('Chip eccentricity'), 'Interpreter','none');
        xlabel('Eccentricity'); ylabel('Count');
    end

    %% Total accumulation
    if strcmp(option, 'Entire Folder')
        total_ball_areas = [total_ball_areas, A_um2(balls)];
        total_fiber_areas = [total_fiber_areas, A_um2(fibers)];
        total_chip_areas = [total_chip_areas, A_um2(chips)];
        total_fiber_lengths = [total_fiber_lengths, major_length_um(fibers)];
        total_ball_circularities = [total_ball_circularities, circ(balls)];
        total_chip_eccentricities = [total_chip_eccentricities, Ecc(chips)];
        total_particle_count = total_particle_count + length(stats);
    end
end

%% Create summary sheet if entire folder
if strcmp(option, 'Entire Folder')
    summary_data = {'Image', 'Balls', 'Fibers', 'Chips', 'Total Particles','Mean Ball Area (µm²)', 'Mean Fiber Area (µm²)', 'Mean Chip Area (µm²)','Scale (µm)', 'Scale (px)'};
    for i = 1:length(total_data)
        dt = total_data(i);
        total_particles = dt.num_balls + dt.num_fibers + dt.num_chips;
        mean_ball_area = mean(dt.ball_areas);
        mean_fiber_area = mean(dt.fiber_areas);
        mean_chip_area = mean(dt.chip_areas);
        
        summary_data{i+1,1} = dt.image_name;
        summary_data{i+1,2} = dt.num_balls;
        summary_data{i+1,3} = dt.num_fibers;
        summary_data{i+1,4} = dt.num_chips;
        summary_data{i+1,5} = total_particles;
        summary_data{i+1,6} = mean_ball_area;
        summary_data{i+1,7} = mean_fiber_area;
        summary_data{i+1,8} = mean_chip_area;
        summary_data{i+1,9} = scale_microns;
        summary_data{i+1,10} = scale_length_px;
    end
 
    xlswrite(excel_file, summary_data, 'Summary');
end

%% Total Plots
if strcmp(option, 'Entire Folder')
    n1 = numel(total_ball_areas);
    n2 = numel(total_fiber_areas);
    n3 = numel(total_chip_areas);
    total = total_particle_count;

    figure('Name','Total - Distribution');
    bar(categorical({'Balls','Fibers','Chips'}), 100*[n1,n2,n3]/total, 'FaceColor', colors{1});
    title('Total particles by type (%)', 'Interpreter','none'); ylabel('%');

    figure('Name','Total - Fiber area');
    histogram(total_fiber_areas, 'FaceColor', colors{2});
    title('Total fiber area', 'Interpreter','none'); xlabel('Area (µm²)'); ylabel('Count');

    figure('Name','Total - Chip area');
    histogram(total_chip_areas, 'FaceColor', colors{3});
    title('Total chip area', 'Interpreter','none'); xlabel('Area (µm²)'); ylabel('Count');

    figure('Name','Total - Ball area');
    histogram(total_ball_areas, 'FaceColor', colors{4});
    title('Total ball area', 'Interpreter','none'); xlabel('Area (µm²)'); ylabel('Count');
    
    figure('Name','Total - Fiber length');
    histogram(total_fiber_lengths, 'FaceColor', colors{5});
    title('Total fiber length', 'Interpreter','none'); xlabel('Length (µm)'); ylabel('Count');

    figure('Name','Total - Ball circularity');
    histogram(total_ball_circularities, 'FaceColor', colors{6});
    title('Total ball circularity', 'Interpreter','none'); xlabel('Circularity'); ylabel('Count');

    figure('Name','Total - Chip eccentricity');
    histogram(total_chip_eccentricities, 'FaceColor', colors{7});
    title('Total chip eccentricity', 'Interpreter','none'); xlabel('Eccentricity'); ylabel('Count');
end

disp('Analysis complete and data exported to Excel.');

%% Tooltip Function
function showTooltip(src, ~, h_particles, particle_data)
    persistent tooltip;

    cp = get(gca, 'CurrentPoint');
    x = cp(1,1); y = cp(1,2);

    found = false;

    for i = 1:length(h_particles)
        h = h_particles(i);
        if strcmp(get(h, 'Visible'), 'on')
            xdata = get(h, 'XData');
            ydata = get(h, 'YData');
            if any(abs(x - xdata) < 5 & abs(y - ydata) < 5)
                found = true;
                data = particle_data{i};
                str = sprintf('ID: %d\nType: %s\nArea: %s\n%s\nPosition: %s', data.ID, data.Type, data.Area, data.Measurement, data.Position);

                if ~isempty(tooltip) && all(ishghandle(tooltip))
                    delete(tooltip);
                end

                tooltip = text(x, y, str, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 5, 'FontSize', 8, 'Tag', 'tooltip');
                return;
            end
        end
    end

    if ~found && ~isempty(tooltip) && all(ishghandle(tooltip))
        delete(tooltip);
        tooltip = [];
    end
end