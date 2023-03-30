%% Paul_OPT_4By4.m (version 3.0)
%Author: Paul Ernst
%Date Created: 3/17/2022
%Date of Last Update: 3/18/2022
%Update History:
%PE 7/20/2022 - Changed for OPT
%PE 3/18/2022 - Added Error plot
%PE 3/17/2022 - Created
%--------------------------------------
%Purpose: Creates 4x4 panel figure of eddy characteristics, comparing different statistics
%Inputs: Eddy characteristics from CE/AE_traj from Win/Sum files
%Outputs: 4x4 panel figure of: SumCE/SumAE/WintWin/WintSum vs. Gen#, Eddy#, Radius, Amplitude
%--------------------------------------
%% Important dynamic inputs - Alter these fields
latlonbounds = [30, -10, 80, 40]; % [N, S, E, W] lat long boundaries
depthlevel = 2;
titlestring = '4by4Plot';
%Date Classification into winter and summer
winMonths = [1, 2, 11, 12];
sumMonths = [5, 6, 7, 8, 9];

%Load trajectories
load('/EXTRACTION/EDDY_TRAJECTORIES/AE_Filtered_Trajectories.mat')
load('/EXTRACTION/EDDY_TRAJECTORIES/CE_Filtered_Trajectories.mat')

%% Load variables & Inputs
%Defaults
load('DefaultData_OPT.mat') %See Paul_FSS_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

%Load colormap and bounds
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);

%% Separate by monsoon
[sumAE, winAE] = monsoonSeparator(AE_trajnew);
[sumCE, winCE] = monsoonSeparator(CE_traj);

%% Getting things we need to plot, using calcEddParams as defined below
gridOpt = 2;
% ROW 1
[genNumSumCE, edNumSumCE, meanRadSumCE, meanAmpSumCE] = calcEddParams(sumCE, X, Y, gridOpt, 0);

% ROW 2
[genNumSumAE, edNumSumAE, meanRadSumAE, meanAmpSumAE] = calcEddParams(sumAE, X, Y, gridOpt, 0);

% ROW 3
[genNumWinCE, edNumWinCE, meanRadWinCE, meanAmpWinCE] = calcEddParams(winCE, X, Y, gridOpt, 1);

% ROW 4
[genNumWinAE, edNumWinAE, meanRadWinAE, meanAmpWinAE] = calcEddParams(winAE, X, Y, gridOpt, 1);

% Place the variables into a master matrix
masterMatrixToPlot = {genNumSumCE, edNumSumCE, meanRadSumCE, meanAmpSumCE;...
    genNumSumAE, edNumSumAE, meanRadSumAE, meanAmpSumAE;...
    genNumWinCE, edNumWinCE, meanRadWinCE, meanAmpWinCE;...
    genNumWinAE, edNumWinAE, meanRadWinAE, meanAmpWinAE};

%% Construct X/Y Grids For Mapping
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%% Plot each field over map geometry
%Initialize Labels
anno = ["(a)", "(b)", "(c)", "(d)";...
    "(e)", "(f)", "(g)", "(h)";...
    "(i)", "(j)", "(k)", "(l)";...
    "(m)", "(n)", "(o)", "(p)";];
ylabs = ["Winter CE", "Winter AE", "Summer CE", "Summer AE"];
xlabs = ["Gen. Num.", "Eddy Num.", "Radius (km)", "Amp. (s^-^1*10^-^4)"];

%Initialize Figure
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.529296875 0.917361111111111], 'Color', 'w'); %fullscreen, invisible figure
t = tiledlayout(4,4,'TileSpacing','none');
load("whiteorangecmap.mat")
%Loop through each tile
maxRows = 4; maxCols = 4;
for nRows = 1:maxRows
    for nCols = 1:maxCols
        %NEXT TILE
        ax = nexttile;

        %MAPPING
        m_proj('mercator','longitude',[WL EL],'latitude',[SL NL]) %initialize map
        hold on
        if nCols == maxCols
            masterMatrixToPlot{nRows,nCols} = masterMatrixToPlot{nRows,nCols} * 1e4;
        end
        m_pcolor(Xgrid,Ygrid,masterMatrixToPlot{nRows,nCols}) %ADT color plot in cm
        shading flat
        m_coast('patch',[.5 .5 .5]);

        %COLORMAP SWITCH
        switch nCols
            % GEN NUM
            case 1
                caxis([0 12]); %0 12
                % NEDDY
            case 2
                caxis([0 160]); %0 80
                % RAD
            case 3
                caxis ([25 100]);
                % AMP
            case 4
                caxis ([0 5]); % 0 0.025
        end

        gridSize = 22;
        xTickInterval = [45 55 65 75];
        yTickInterval = [-8 0 8 16 24];
        % GRID WITH NO X BUT Y
        if ((nRows ~= maxRows) && (nCols == 1))
            m_grid('box', 'fancy','fontsize',gridSize, 'xticklabels',[], 'xtick', xTickInterval, 'ytick', yTickInterval);
            % GRID WITH NO Y BUT X
        elseif ((nRows == maxRows) && (nCols > 1))
            m_grid('box', 'fancy','fontsize',gridSize, 'yticklabels',[], 'xtick', xTickInterval, 'ytick', yTickInterval);
            % GRID WITH BOTH
        elseif ((nRows == maxRows) && (nCols == 1))
            m_grid('box', 'fancy','fontsize',gridSize, 'xtick', xTickInterval, 'ytick', yTickInterval);
            % GRID WITH NEITHER
        else
            m_grid('box', 'fancy','fontsize',gridSize, 'yticklabels',[],'xticklabels',[], 'xtick', xTickInterval, 'ytick', yTickInterval);
        end

        %COLORBARS
        if nRows == 4
            colorbar('southoutside')
            colormap(CustomColormap);
        end

        % FONTSIZE
        set(gca,'fontsize',28);

        % ANNOTATION
        m_text(44,22,anno(nRows,nCols),'color','w','fontsize',28);

        % LABELS - Y
        if nCols == 1
            ylabel(ylabs(nRows));
        end

        % LABELS - X
        if nRows == 1
            title(xlabs(nCols));
        end
    end
end

%% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG

%% Separate trajectories into winter and summer eddies
%14 = datenum
function [winTraj, sumTraj] = monsoonSeparator(trajInput)
winTraj = cell(length(trajInput),21);
sumTraj = cell(length(trajInput),21);
winMonths = [1, 2, 11, 12];
for i = 1:length(trajInput)
    date_num = trajInput{i,14}(1);
    %date num = first eddy thing
    date2use = datestr(date_num,'mm/dd/yyyy');
    month2use = str2double(date2use(1:2));

    %Calibrate date to use - winter or summer
    if ismember(month2use,winMonths)
        winTraj(i,:) = trajInput(i,:);
    else
        sumTraj(i,:) = trajInput(i,:);
    end
end

winTraj = winTraj(~cellfun('isempty',winTraj(:,1)),:);
sumTraj = sumTraj(~cellfun('isempty',sumTraj(:,1)),:);
end

%% Function for fetching each line of the plot (CE/AE/Win/Sum)
%Effectively, we loop thru each eddy state/date and grab the params in each pixel
%Slot params into correct places and spit them back out averaged over timewise area
%gridOpt options: 0 = scaled by radius; 1 = naked points; 2 = boxed into 1 degree boxes
function [genNum, edNum, meanRad, meanAmp] = calcEddParams(trajList, X, Y, gridOpt, winOrSum)

%If this is winter or summer monsoon
if winOrSum
    daysBetween = 120;
else
    daysBetween = 153;
end

%Initialize our end matrices
genNum = zeros(length(X),length(Y)); %Simple sums across geography
edNum = zeros(length(X),length(Y)); %''
meanRad = zeros(length(X),length(Y)); %Averaged temporally
meanAmp = zeros(length(X),length(Y)); %''
eddCheck = zeros(length(X),length(Y)); %For when an eddy is in the area

%Subset our master matrix
idlist = trajList(:,1); x_e = trajList(:,4); y_e = trajList(:,5);
rad_e = trajList(:,6); amp_e = trajList(:,8); date_e = trajList(:,14);
% ampcount = 0; allamps = 0; allamps2 = [];
%loop through all ID'd eddies
for i = 1:length(trajList)
    %extract components for one trajectory
    x = x_e{i,1};
    y = y_e{i,1};
    rad = rad_e{i,1};
    amp = amp_e{i,1};

    %loop thru this particular eddy trajectory
    for j = 1:length(x)
        %subset for this instance of this trajectory
        currX = x(j);
        currY = y(j);
        currRad = rad(j);
        currAmp = amp(j);

        %find points covered by this eddy
        if gridOpt == 0 %scaled by radius
            latoffset = currRad/110.573;
            lonoffset = currRad/(cosd(currY)*111.32);
            NLed = nearestpoint(double(currY+(latoffset)), Y);
            SLed = nearestpoint(double(currY-(latoffset)), Y);
            ELed = nearestpoint(double(currX+(lonoffset)), X);
            WLed = nearestpoint(double(currX-(lonoffset)), X);
            xToSum = WLed:ELed;
            yToSum = SLed:NLed;
        else %naked points or boxed
            xToSum = nearestpoint(double(currX), X);
            yToSum = nearestpoint(double(currY), Y);
        end

        %Gen Number Slotting
        if j == 1 %First iteration we ID this eddy
            genNum(xToSum,yToSum) = genNum(xToSum,yToSum) + 1;
        end
        eddCheck(xToSum,yToSum) = eddCheck(xToSum,yToSum) + 1;

        %Ed Number Slotting
        edNum(xToSum,yToSum) = edNum(xToSum,yToSum) + 1;

        %Rad Slotting
        meanRad(xToSum,yToSum) = meanRad(xToSum,yToSum) + currRad;

        %Amp Slotting
        meanAmp(xToSum,yToSum) = meanAmp(xToSum,yToSum) + (currAmp-gsw_f(currY)); %anomaly

    end %End this eddy
end %End all eddies

% Perform means across time
if gridOpt ~= 2
    meanRad = meanRad./daysBetween;
    meanAmp = meanAmp./daysBetween;
end

if gridOpt == 2 %boxed into 1-degree sums
    genNumTemp = zeros(size(genNum));
    edNumTemp = zeros(size(edNum));
    meanRadTemp = zeros(size(meanRad));
    meanAmpTemp = zeros(size(meanAmp));
    gridSize = 24; %1/12 degree resolution
    gridSizeMin = gridSize-1;
    for i = 0:gridSize:length(X)-gridSize
        for j = 0:gridSize:length(Y)-gridSize
            if i == 0
                i = 1;
            end
            if j == 0
                j = 1;
            end
            %Gridded sums
            genNumTemp(i:i+gridSizeMin, j:j+gridSizeMin) =...
                ones(gridSize,gridSize)*nansum(nansum(genNum(i:i+gridSizeMin, j:j+gridSizeMin)));
            edNumTemp(i:i+gridSizeMin, j:j+gridSizeMin) =...
                ones(gridSize,gridSize)*nansum(nansum(edNum(i:i+gridSizeMin, j:j+gridSizeMin)));
            meanRadTemp(i:i+gridSizeMin, j:j+gridSizeMin) =...
                ones(gridSize,gridSize)*(nansum(nansum(meanRad(i:i+gridSizeMin, j:j+gridSizeMin)))...
                /nansum(nansum(eddCheck(i:i+gridSizeMin, j:j+gridSizeMin)))); %sum across area then divide by num of eddies for total mean
            meanAmpTemp(i:i+gridSizeMin, j:j+gridSizeMin) =...
                ones(gridSize,gridSize)*(nansum(nansum(meanAmp(i:i+gridSizeMin, j:j+gridSizeMin)))...
                /nansum(nansum(eddCheck(i:i+gridSizeMin, j:j+gridSizeMin))));
        end
    end
    %Slot back in
    genNum = genNumTemp;
    edNum = edNumTemp;
    meanRad = meanRadTemp;
    meanAmp = meanAmpTemp;
end %End boxing
genNum(isnan(genNum)) = 0;
edNum(isnan(edNum)) = 0;
meanRad(isnan(meanRad)) = 0;
meanAmp(isnan(meanAmp)) = 0;
end %End eddy parameter calculating function

