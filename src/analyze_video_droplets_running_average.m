function results = analyze_video_droplets_running_average()
%% Droplet Freezing Video Analysis
%
% MATLAB workflow for analyzing droplet freezing experiments recorded on a
% cold plate.
%
% The script:
%   1. Loads a video file.
%   2. Selects the time interval to analyze.
%   3. Selects the cold-plate region of interest (ROI).
%   4. Detects droplets using Cellpose.
%   5. Allows manual correction of droplet detections.
%   6. Builds a background mask excluding droplets.
%   7. Normalizes each frame using the cold-plate background.
%   8. Computes running-average frames over fixed time windows.
%   9. Extracts droplet intensity and contrast metrics.
%  10. Saves tables, plots, and detection images.
%
% Outputs:
%   results.dropletInfoTable
%   results.runningTable
%   results.summaryTable
%
% Author: Jorge H. Melillo

close all;
clc;

%% =================== User parameters ===================

cellposeModel = 'cyto';          % Cellpose model used for droplet detection
defaultWindowSec = 2;            % Default running-average window [s]

minAreaRadiusFactor = 0.25;      % Minimum radius factor for area filtering
maxAreaRadiusFactor = 1.50;      % Maximum radius factor for area filtering

minDiameterFactor = 0.50;        % Minimum equivalent diameter factor
maxDiameterFactor = 1.60;        % Maximum equivalent diameter factor

brightPercentile = 90;           % Bright-pixel threshold percentile

%% =================== Select video ===================

[videoName, videoPath] = uigetfile( ...
    {'*.mp4;*.avi;*.mov;*.m4v', 'Video files'}, ...
    'Select video');

if isequal(videoName, 0)
    error('No video selected.');
end

videoFile = fullfile(videoPath, videoName);

%% =================== Select output folder ===================

outputFolder = uigetdir(videoPath, 'Select folder to save results');

if isequal(outputFolder, 0)
    error('No output folder selected.');
end

%% =================== Read video information ===================

videoObj = VideoReader(videoFile);

fps = videoObj.FrameRate;
durationSec = videoObj.Duration;
nFramesApprox = floor(durationSec * fps);

fprintf('Video: %s\n', videoFile);
fprintf('Duration: %.3f s\n', durationSec);
fprintf('FPS: %.3f\n', fps);
fprintf('Approximate number of frames: %d\n', nFramesApprox);

%% =================== Choose analysis interval ===================

answerPortion = inputdlg( ...
    {'Start time (s):', 'End time (s):', 'Running-average window (s):'}, ...
    'Video portion', ...
    [1 50; 1 50; 1 50], ...
    {'0', sprintf('%.3f', durationSec), num2str(defaultWindowSec)});

if isempty(answerPortion)
    error('No video portion selected.');
end

tStart = str2double(answerPortion{1});
tEnd = str2double(answerPortion{2});
windowSec = str2double(answerPortion{3});

if isnan(tStart) || isnan(tEnd) || isnan(windowSec)
    error('Invalid numeric input.');
end

if tStart < 0 || tEnd <= tStart || tEnd > durationSec
    error('Invalid start/end times.');
end

if windowSec <= 0
    error('Running-average window must be larger than zero.');
end

%% =================== Read first analysis frame ===================

videoObj.CurrentTime = tStart;

if ~hasFrame(videoObj)
    error('Could not read first frame at selected start time.');
end

firstFrame = readFrame(videoObj);
firstGray = convertToGray(firstFrame);

%% =================== Select cold-plate ROI ===================

figROI = figure('Color', 'w', 'Name', 'ROI selection');
imshow(firstFrame);
title('Draw rectangle around the cold-plate region');

roiHandle = drawrectangle('Color', 'b');
wait(roiHandle);

coldPlateROI = round(roiHandle.Position);

hold on;
rectangle('Position', coldPlateROI, 'EdgeColor', 'b', 'LineWidth', 1.5);
hold off;

uiwait(msgbox('Cold-plate ROI selected. Press OK to continue.', ...
    'Continue', 'modal'));

close(figROI);

%% =================== Crop first cold-plate patch ===================

coldPlatePatch = imcrop(firstFrame, coldPlateROI);

%% =================== Estimate droplet diameter ===================

figDiameter = figure('Color', 'w', 'Name', 'Select droplet diameter');
imshow(coldPlatePatch);
title('Draw a line across one typical droplet diameter');

lineHandle = drawline('Color', 'g');
wait(lineHandle);

linePosition = lineHandle.Position;
estimatedDropletDiameter_px = sqrt(sum((linePosition(2,:) - linePosition(1,:)).^2));

close(figDiameter);

fprintf('Estimated droplet diameter: %.2f px\n', estimatedDropletDiameter_px);

%% =================== Cellpose droplet segmentation ===================

grayPatch = im2gray(coldPlatePatch);
invertedGrayPatch = imcomplement(grayPatch);

cellposeObject = cellpose(Model=cellposeModel);

labels = segmentCells2D(cellposeObject, invertedGrayPatch, ...
    ImageCellDiameter=estimatedDropletDiameter_px);

stats = regionprops(labels, ...
    'Area', ...
    'Centroid', ...
    'EquivDiameter', ...
    'BoundingBox');

if isempty(stats)
    warning('No droplets detected by Cellpose.');
end

%% =================== Filter automatic detections ===================

nDetected = numel(stats);

minArea = pi * (minAreaRadiusFactor * estimatedDropletDiameter_px)^2;
maxArea = pi * (maxAreaRadiusFactor * estimatedDropletDiameter_px)^2;

keepDetection = false(nDetected, 1);

for k = 1:nDetected

    objectArea = stats(k).Area;
    objectDiameter = stats(k).EquivDiameter;

    if objectArea >= minArea && objectArea <= maxArea && ...
       objectDiameter >= minDiameterFactor * estimatedDropletDiameter_px && ...
       objectDiameter <= maxDiameterFactor * estimatedDropletDiameter_px

        keepDetection(k) = true;
    end
end

stats = stats(keepDetection);
nDroplets = numel(stats);

%% =================== Build automatic detection table ===================

DetectionType = repmat('Auto', nDroplets, 1);
DropletID = (1:nDroplets)';

CenterX_px = zeros(nDroplets, 1);
CenterY_px = zeros(nDroplets, 1);
Area_px2 = zeros(nDroplets, 1);
EquivDiameter_px = zeros(nDroplets, 1);

BoundingBoxX_px = zeros(nDroplets, 1);
BoundingBoxY_px = zeros(nDroplets, 1);
BoundingBoxW_px = zeros(nDroplets, 1);
BoundingBoxH_px = zeros(nDroplets, 1);

for k = 1:nDroplets

    centroid = stats(k).Centroid;
    boundingBox = stats(k).BoundingBox;

    CenterX_px(k) = centroid(1) + coldPlateROI(1) - 1;
    CenterY_px(k) = centroid(2) + coldPlateROI(2) - 1;

    BoundingBoxX_px(k) = boundingBox(1) + coldPlateROI(1) - 1;
    BoundingBoxY_px(k) = boundingBox(2) + coldPlateROI(2) - 1;
    BoundingBoxW_px(k) = boundingBox(3);
    BoundingBoxH_px(k) = boundingBox(4);

    Area_px2(k) = stats(k).Area;
    EquivDiameter_px(k) = stats(k).EquivDiameter;
end

%% =================== Display automatic detections ===================

figAuto = figure('Color', 'w', 'Name', 'Automatic detections');
imshow(firstFrame);
hold on;

rectangle('Position', coldPlateROI, 'EdgeColor', 'b', 'LineWidth', 1.5);

for k = 1:nDroplets

    rectangle('Position', ...
        [BoundingBoxX_px(k), BoundingBoxY_px(k), BoundingBoxW_px(k), BoundingBoxH_px(k)], ...
        'EdgeColor', 'y', ...
        'LineWidth', 1.2);

    text(CenterX_px(k), CenterY_px(k), sprintf('%d', DropletID(k)), ...
        'Color', 'y', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

title(sprintf('Automatic detections: %d droplets', nDroplets));
hold off;

uiwait(msgbox('Check automatic detections, then press OK to continue.', ...
    'Continue', 'modal'));

close(figAuto);

%% =================== Optionally remove false detections ===================

deleteAnswer = questdlg('Is there any false droplet to delete?', ...
    'Delete false droplets', ...
    'Yes', 'No', 'No');

if strcmp(deleteAnswer, 'Yes')

    figRemove = figure('Color', 'w', 'Name', 'Remove false detections');
    imshow(firstFrame);
    hold on;

    rectangle('Position', coldPlateROI, 'EdgeColor', 'b', 'LineWidth', 1.5);

    for k = 1:nDroplets

        rectangle('Position', ...
            [BoundingBoxX_px(k), BoundingBoxY_px(k), BoundingBoxW_px(k), BoundingBoxH_px(k)], ...
            'EdgeColor', 'y', ...
            'LineWidth', 1.2);

        text(CenterX_px(k), CenterY_px(k), sprintf('%d', DropletID(k)), ...
            'Color', 'y', ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
    end

    title({'Automatic detections', ...
           'Check IDs and enter in Command Window which ones to remove'});
    hold off;

    removeIDs = input('Enter IDs to REMOVE as a vector, e.g. [2 5 8]: ');

    if isempty(removeIDs)
        keepManual = true(nDroplets, 1);
    else
        keepManual = ~ismember(DropletID, removeIDs);
    end

    DetectionType = DetectionType(keepManual);
    CenterX_px = CenterX_px(keepManual);
    CenterY_px = CenterY_px(keepManual);
    Area_px2 = Area_px2(keepManual);
    EquivDiameter_px = EquivDiameter_px(keepManual);
    BoundingBoxX_px = BoundingBoxX_px(keepManual);
    BoundingBoxY_px = BoundingBoxY_px(keepManual);
    BoundingBoxW_px = BoundingBoxW_px(keepManual);
    BoundingBoxH_px = BoundingBoxH_px(keepManual);

    close(figRemove);
end

%% =================== Optionally add missing droplets ===================

addAnswer = questdlg('Is there any droplet that was not detected?', ...
    'Add missing droplets', ...
    'Yes', 'No', 'No');

if strcmp(addAnswer, 'Yes')

    nCurrent = numel(CenterX_px);

    figAdd = figure('Color', 'w', 'Name', 'Add missing droplets');
    imshow(firstFrame);
    hold on;

    rectangle('Position', coldPlateROI, 'EdgeColor', 'b', 'LineWidth', 1.5);

    for k = 1:nCurrent

        rectangle('Position', ...
            [BoundingBoxX_px(k), BoundingBoxY_px(k), BoundingBoxW_px(k), BoundingBoxH_px(k)], ...
            'EdgeColor', 'y', ...
            'LineWidth', 1.2);

        text(CenterX_px(k), CenterY_px(k), sprintf('%d', k), ...
            'Color', 'y', ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
    end

    title({'Draw circles on missing droplets', ...
           'After each circle double-click to confirm. Press No when finished.'});

    manualCenters = [];
    manualRadii = [];

    answer = 'Yes';

    while strcmpi(answer, 'Yes')

        circleHandle = drawcircle('Color', 'g');
        wait(circleHandle);

        manualCenters = [manualCenters; circleHandle.Center]; %#ok<AGROW>
        manualRadii = [manualRadii; circleHandle.Radius]; %#ok<AGROW>

        viscircles(circleHandle.Center, circleHandle.Radius, ...
            'Color', 'g', ...
            'LineWidth', 1);

        answer = questdlg('Add another missing droplet?', ...
            'Manual addition', ...
            'Yes', 'No', 'No');
    end

    close(figAdd);

    nManual = size(manualCenters, 1);

    for k = 1:nManual

        cx = manualCenters(k, 1);
        cy = manualCenters(k, 2);
        r = manualRadii(k);

        DetectionType(end+1, 1) = "Manual";
        CenterX_px(end+1, 1) = cx;
        CenterY_px(end+1, 1) = cy;
        Area_px2(end+1, 1) = pi * r^2;
        EquivDiameter_px(end+1, 1) = 2 * r;
        BoundingBoxX_px(end+1, 1) = cx - r;
        BoundingBoxY_px(end+1, 1) = cy - r;
        BoundingBoxW_px(end+1, 1) = 2 * r;
        BoundingBoxH_px(end+1, 1) = 2 * r;
    end
end

%% =================== Rebuild IDs after corrections ===================

nFinal = numel(CenterX_px);
DropletID = (1:nFinal)';

%% =================== Save corrected detection image ===================

figDet = figure('Visible', 'off', 'Color', 'w', ...
    'Name', 'Corrected droplets - first frame');

imshow(firstFrame);
hold on;

rectangle('Position', coldPlateROI, 'EdgeColor', 'b', 'LineWidth', 1.5);

for k = 1:nFinal

    if DetectionType(k) == "Auto"
        boxColor = 'g';
        textColor = 'y';
    else
        boxColor = 'c';
        textColor = 'c';
    end

    rectangle('Position', ...
        [BoundingBoxX_px(k), BoundingBoxY_px(k), BoundingBoxW_px(k), BoundingBoxH_px(k)], ...
        'EdgeColor', boxColor, ...
        'LineWidth', 1.2);

    text(CenterX_px(k), CenterY_px(k), sprintf('%d', DropletID(k)), ...
        'Color', textColor, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

title(sprintf('Detected droplets on first frame | N = %d', nFinal));
hold off;

%% =================== Droplet and background masks ===================

dropletMaskFull = false(size(firstGray));

[imageHeight, imageWidth] = size(firstGray);
[xx, yy] = meshgrid(1:imageWidth, 1:imageHeight);

for k = 1:nFinal

    cx = CenterX_px(k);
    cy = CenterY_px(k);

    radius_px = 0.5 * EquivDiameter_px(k);

    if isnan(radius_px) || radius_px <= 0
        radius_px = 0.25 * (BoundingBoxW_px(k) + BoundingBoxH_px(k));
    end

    dropletMaskFull = dropletMaskFull | ...
        ((xx - cx).^2 + (yy - cy).^2 <= radius_px.^2);
end

coldPlateMaskFull = false(size(firstGray));

x1 = coldPlateROI(1);
y1 = coldPlateROI(2);
w1 = coldPlateROI(3);
h1 = coldPlateROI(4);

coldPlateMaskFull(y1:y1+h1-1, x1:x1+w1-1) = true;

backgroundMask = coldPlateMaskFull & ~dropletMaskFull;

backgroundFirst = double(firstGray(backgroundMask));

if isempty(backgroundFirst)
    error('Background mask is empty in first frame.');
end

backgroundMeanFirst = mean(backgroundFirst, 'omitnan');
backgroundStdFirst = std(backgroundFirst, 'omitnan');

if backgroundStdFirst == 0 || isnan(backgroundStdFirst)
    error('Background standard deviation in first frame is zero or invalid.');
end

%% =================== Running-average setup ===================

framesPerWindow = max(1, round(windowSec * fps));

startFrame = max(1, floor(tStart * fps) + 1);
endFrame = min(nFramesApprox, floor(tEnd * fps));

if endFrame < startFrame
    error('Selected time range contains no frames.');
end

frameIndexList = startFrame:endFrame;

nFramesSelected = numel(frameIndexList);
nWindows = ceil(nFramesSelected / framesPerWindow);

fprintf('Selected frames: %d to %d (%d frames)\n', ...
    startFrame, endFrame, nFramesSelected);

fprintf('Frames per averaging window: %d\n', framesPerWindow);
fprintf('Number of averaging windows: %d\n', nWindows);

%% =================== Initialize running metrics ===================

rowsTotal = nFinal * nWindows;

DropletID_col = zeros(rowsTotal, 1);
DetectionType_col = strings(rowsTotal, 1);
WindowIndex_col = zeros(rowsTotal, 1);
FrameStart_col = zeros(rowsTotal, 1);
FrameEnd_col = zeros(rowsTotal, 1);
TimeStart_s_col = zeros(rowsTotal, 1);
TimeCenter_s_col = zeros(rowsTotal, 1);
TimeEnd_s_col = zeros(rowsTotal, 1);

MeanGray_col = zeros(rowsTotal, 1);
ContrastStd_col = zeros(rowsTotal, 1);
BrightMean_col = zeros(rowsTotal, 1);
BrightFrac_col = zeros(rowsTotal, 1);

BGMean_raw_col = zeros(rowsTotal, 1);
BGStd_raw_col = zeros(rowsTotal, 1);
BGMean_norm_col = zeros(rowsTotal, 1);
BGStd_norm_col = zeros(rowsTotal, 1);

rowCounter = 0;

%% =================== Loop over averaging windows ===================

for windowIndex = 1:nWindows

    idxA = (windowIndex - 1) * framesPerWindow + 1;
    idxB = min(windowIndex * framesPerWindow, nFramesSelected);

    frameWindow = frameIndexList(idxA:idxB);
    nFramesInWindow = numel(frameWindow);

    accumulatedNormalizedFrame = zeros(size(firstGray), 'double');

    bgMeanRawValues = zeros(nFramesInWindow, 1);
    bgStdRawValues = zeros(nFramesInWindow, 1);
    bgMeanNormValues = zeros(nFramesInWindow, 1);
    bgStdNormValues = zeros(nFramesInWindow, 1);

    for j = 1:nFramesInWindow

        frameNumber = frameWindow(j);

        videoObj.CurrentTime = (frameNumber - 1) / fps;

        if ~hasFrame(videoObj)
            error('Could not read frame %d.', frameNumber);
        end

        currentFrame = readFrame(videoObj);
        currentGray = convertToGray(currentFrame);

        backgroundCurrent = double(currentGray(backgroundMask));

        if isempty(backgroundCurrent)
            error('Background mask is empty in frame %d.', frameNumber);
        end

        bgMean = mean(backgroundCurrent, 'omitnan');
        bgStd = std(backgroundCurrent, 'omitnan');

        if bgStd == 0 || isnan(bgStd)
            bgStd = 1;
        end

        normalizedGray = double(currentGray);
        normalizedGray = ...
            (normalizedGray - bgMean) * (backgroundStdFirst / bgStd) + backgroundMeanFirst;

        normalizedGray = max(0, min(255, normalizedGray));

        normalizedBackground = normalizedGray(backgroundMask);

        bgMeanRawValues(j) = bgMean;
        bgStdRawValues(j) = bgStd;
        bgMeanNormValues(j) = mean(normalizedBackground, 'omitnan');
        bgStdNormValues(j) = std(normalizedBackground, 'omitnan');

        accumulatedNormalizedFrame = accumulatedNormalizedFrame + normalizedGray;
    end

    averageFrame = accumulatedNormalizedFrame / nFramesInWindow;

    frameStart = frameWindow(1);
    frameEnd = frameWindow(end);

    timeStart = (frameStart - 1) / fps;
    timeEnd = (frameEnd - 1) / fps;
    timeCenter = 0.5 * (timeStart + timeEnd);

    for k = 1:nFinal

        rowCounter = rowCounter + 1;

        cx = CenterX_px(k);
        cy = CenterY_px(k);

        radius_px = 0.5 * EquivDiameter_px(k);

        if isnan(radius_px) || radius_px <= 0
            radius_px = 0.25 * (BoundingBoxW_px(k) + BoundingBoxH_px(k));
        end

        dropletMask = (xx - cx).^2 + (yy - cy).^2 <= radius_px.^2;

        dropletPixels = double(averageFrame(dropletMask));

        meanGray = mean(dropletPixels, 'omitnan');
        contrastStd = std(dropletPixels, 'omitnan');

        brightThreshold = prctile(dropletPixels, brightPercentile);
        brightPixels = dropletPixels(dropletPixels >= brightThreshold);

        if isempty(brightPixels)
            brightMean = NaN;
        else
            brightMean = mean(brightPixels, 'omitnan');
        end

        brightFraction = sum(dropletPixels >= brightThreshold) / numel(dropletPixels);

        DropletID_col(rowCounter) = DropletID(k);
        DetectionType_col(rowCounter) = DetectionType(k);
        WindowIndex_col(rowCounter) = windowIndex;
        FrameStart_col(rowCounter) = frameStart;
        FrameEnd_col(rowCounter) = frameEnd;
        TimeStart_s_col(rowCounter) = timeStart;
        TimeCenter_s_col(rowCounter) = timeCenter;
        TimeEnd_s_col(rowCounter) = timeEnd;

        MeanGray_col(rowCounter) = meanGray;
        ContrastStd_col(rowCounter) = contrastStd;
        BrightMean_col(rowCounter) = brightMean;
        BrightFrac_col(rowCounter) = brightFraction;

        BGMean_raw_col(rowCounter) = mean(bgMeanRawValues, 'omitnan');
        BGStd_raw_col(rowCounter) = mean(bgStdRawValues, 'omitnan');
        BGMean_norm_col(rowCounter) = mean(bgMeanNormValues, 'omitnan');
        BGStd_norm_col(rowCounter) = mean(bgStdNormValues, 'omitnan');
    end
end

%% =================== Build droplet information table ===================

dropletInfoTable = table( ...
    DropletID, ...
    DetectionType, ...
    CenterX_px, ...
    CenterY_px, ...
    Area_px2, ...
    EquivDiameter_px, ...
    BoundingBoxX_px, ...
    BoundingBoxY_px, ...
    BoundingBoxW_px, ...
    BoundingBoxH_px, ...
    repmat(string(videoName), nFinal, 1), ...
    repmat(tStart, nFinal, 1), ...
    repmat(tEnd, nFinal, 1), ...
    'VariableNames', { ...
    'DropletID', ...
    'DetectionType', ...
    'CenterX_px', ...
    'CenterY_px', ...
    'Area_px2', ...
    'EquivDiameter_px', ...
    'BoundingBoxX_px', ...
    'BoundingBoxY_px', ...
    'BoundingBoxW_px', ...
    'BoundingBoxH_px', ...
    'VideoFile', ...
    'AnalysisStart_s', ...
    'AnalysisEnd_s'});

%% =================== Build running metrics table ===================

runningTable = table( ...
    repmat(string(videoName), rowsTotal, 1), ...
    DropletID_col, ...
    DetectionType_col, ...
    WindowIndex_col, ...
    FrameStart_col, ...
    FrameEnd_col, ...
    TimeStart_s_col, ...
    TimeCenter_s_col, ...
    TimeEnd_s_col, ...
    MeanGray_col, ...
    ContrastStd_col, ...
    BrightMean_col, ...
    BrightFrac_col, ...
    BGMean_raw_col, ...
    BGStd_raw_col, ...
    BGMean_norm_col, ...
    BGStd_norm_col, ...
    'VariableNames', { ...
    'VideoFile', ...
    'DropletID', ...
    'DetectionType', ...
    'WindowIndex', ...
    'FrameStart', ...
    'FrameEnd', ...
    'TimeStart_s', ...
    'TimeCenter_s', ...
    'TimeEnd_s', ...
    'MeanGray', ...
    'ContrastStd', ...
    'BrightMean', ...
    'BrightFrac', ...
    'BGMean_raw', ...
    'BGStd_raw', ...
    'BGMean_norm', ...
    'BGStd_norm'});

%% =================== Build summary table ===================

summaryTable = table( ...
    string(videoName), ...
    tStart, ...
    tEnd, ...
    fps, ...
    windowSec, ...
    framesPerWindow, ...
    nFramesSelected, ...
    nWindows, ...
    nFinal, ...
    backgroundMeanFirst, ...
    backgroundStdFirst, ...
    'VariableNames', { ...
    'VideoFile', ...
    'AnalysisStart_s', ...
    'AnalysisEnd_s', ...
    'FPS', ...
    'Window_s', ...
    'FramesPerWindow', ...
    'FramesAnalyzed', ...
    'NumWindows', ...
    'NumDroplets', ...
    'BGMean_first', ...
    'BGStd_first'});

disp(dropletInfoTable);
disp(summaryTable);

%% =================== Plot BrightMean vs time ===================

figBright = figure('Visible', 'off', 'Color', 'w', ...
    'Name', 'BrightMean vs averaged time');

hold on;

for k = 1:nFinal
    idx = runningTable.DropletID == DropletID(k);

    plot(runningTable.TimeCenter_s(idx), ...
         runningTable.BrightMean(idx), ...
         '-o', ...
         'LineWidth', 1.2);
end

xlabel('Averaged time (s)');
ylabel('BrightMean');
title('BrightMean of each droplet vs averaged time');
grid on;
legend(compose('Droplet %d', DropletID), 'Location', 'eastoutside');
hold off;

%% =================== Plot ContrastStd vs time ===================

figContrast = figure('Visible', 'off', 'Color', 'w', ...
    'Name', 'ContrastStd vs averaged time');

hold on;

for k = 1:nFinal
    idx = runningTable.DropletID == DropletID(k);

    plot(runningTable.TimeCenter_s(idx), ...
         runningTable.ContrastStd(idx), ...
         '-o', ...
         'LineWidth', 1.2);
end

xlabel('Averaged time (s)');
ylabel('ContrastStd');
title('ContrastStd of each droplet vs averaged time');
grid on;
legend(compose('Droplet %d', DropletID), 'Location', 'eastoutside');
hold off;

%% =================== Save outputs ===================

[~, baseVideo, ~] = fileparts(videoName);

outDetectionImage = fullfile(outputFolder, ...
    [baseVideo '_first_frame_detections.png']);

outDropletInfoCSV = fullfile(outputFolder, ...
    [baseVideo '_droplet_info.csv']);

outRunningMetricsCSV = fullfile(outputFolder, ...
    [baseVideo '_running_average_metrics.csv']);

outSummaryCSV = fullfile(outputFolder, ...
    [baseVideo '_summary.csv']);

outBrightPlot = fullfile(outputFolder, ...
    [baseVideo '_brightmean_vs_time.png']);

outContrastPlot = fullfile(outputFolder, ...
    [baseVideo '_contrast_vs_time.png']);

outMAT = fullfile(outputFolder, ...
    [baseVideo '_running_average_results.mat']);

exportgraphics(figDet, outDetectionImage, 'Resolution', 300);
exportgraphics(figBright, outBrightPlot, 'Resolution', 300);
exportgraphics(figContrast, outContrastPlot, 'Resolution', 300);

writetable(dropletInfoTable, outDropletInfoCSV);
writetable(runningTable, outRunningMetricsCSV);
writetable(summaryTable, outSummaryCSV);

save(outMAT, ...
    'dropletInfoTable', ...
    'runningTable', ...
    'summaryTable', ...
    'coldPlateROI', ...
    'backgroundMask', ...
    'estimatedDropletDiameter_px', ...
    'videoFile', ...
    'tStart', ...
    'tEnd', ...
    'windowSec');

close(figDet);
close(figBright);
close(figContrast);

fprintf('\nProcessing finished.\n');
fprintf('Droplets detected: %d\n', nFinal);
fprintf('Running windows: %d\n', nWindows);
fprintf('Detection image saved to:\n%s\n', outDetectionImage);
fprintf('Droplet info table saved to:\n%s\n', outDropletInfoCSV);
fprintf('Running metrics table saved to:\n%s\n', outRunningMetricsCSV);
fprintf('Summary table saved to:\n%s\n', outSummaryCSV);
fprintf('BrightMean plot saved to:\n%s\n', outBrightPlot);
fprintf('Contrast plot saved to:\n%s\n', outContrastPlot);

%% =================== Return results ===================

results.dropletInfoTable = dropletInfoTable;
results.runningTable = runningTable;
results.summaryTable = summaryTable;

end

%% =================== Local function ===================

function grayImage = convertToGray(inputImage)

    if size(inputImage, 3) == 3
        grayImage = rgb2gray(inputImage);
    else
        grayImage = inputImage;
    end
end