%% Paul_OPT_DisplayOptimizedEdgeSmoothAndRatio.m (version 1.1)
%Author: Paul Ernst
%Date Created: 5/27/2022
%Date of Last Update: 9/23/2022
%Update History:
%PE 9/23/2022 - Cleaned for final release
%PE 5/27/2022 - Created
%--------------------------------------
%Purpose: Displays the optimized smoothing and ratios of all edge methods (Figure 4)
%Inputs: Plot sources from Comparison_Main
%Outputs: 5 panels with the following lines:
% -Edge LNAM: SMTH, RATIO    -Edge VORT_T: SMTH, RATIO     -Edge VORT_WA: SMTH, RATIO      -Edge
% OW_T: SMTH, RATIO
% -Edge OW_WA: SMTH, RATIO    -Edge PV_ISO: SMTH, RATIO
%% Inputs
titlestring = "OptimizedEdgeSMTHandRatio";
PlotSources = ["OPT_E_VORT_T_SMTHRAT" "OPT_E_VORT_WA_SMTHRAT" "OPT_E_OW_T_SMTHRAT" "OPT_E_OW_WA_SMTHRAT" "OPT_E_PV_ISO_SMTHRAT"];
xLabels = ["Smoothing Factor", "Smoothing Factor", "Smoothing Factor", "Smoothing Factor", "Smoothing Factor"];
smoothingFactor = 0:15; %what we use in our testing
minRatio = 0:0.02:0.3; %''
load("whiteorangecmap.mat")
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
%% Plots

numPanels = length(PlotSources);
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.64765625 0.917361111111111], 'Color', 'w'); %fullscreen, invisible figure
t = tiledlayout(3,2,'TileSpacing','tight');

panelLabels = ["(a) VORT_T","(b) VORT_W_A","(c) OW_T", "(d) OW_W_A", "(e) PV_I_S_O"];
messageToPrint = strings(5,1);
%Loop thru plots
for nTile = 1:numPanels
    ax = nexttile();

    %load the source of this panel
    thisSource = ["IO_" + PlotSources(nTile) + "_Comparisons.mat"];
    load(thisSource)

    %Grab maximum for highlight later
    [MaxVal,MaxInd] = max(mean(totalSims,1));
    pcolToPlot = reshape(mean(totalSims), [length(smoothingFactor),length(minRatio)]); % x = smooth, y = ratio
    pcolToPlot(isnan(pcolToPlot)) = 0;
    [a, b] = ind2sub([length(smoothingFactor),length(minRatio)],MaxInd);
    maxIndX = smoothingFactor(b);
    maxIndY = minRatio(a);

    %Plot individual panels
    ss = pcolor(smoothingFactor,minRatio,pcolToPlot);
    ss.FaceColor = 'texturemap';
    shading flat
    hold on
    grid on
    xTint = 15/16;
    yTint = 0.3/0.32;
    scatter(maxIndX*xTint+xTint/2,maxIndY*yTint+0.3/16*.5,2000,'kx','linewidth',3)
    xticks(xTint/2:xTint:16)
    xticklabels(0:15)

    colormap(CustomColormap)
    colorbar
    caxis([40 70])

    %labels and limits
    title(panelLabels(nTile))
    ylabel("Centroid Ratio")
    xlabel(xLabels(nTile))
    set(ax,'fontsize',24);
    grid on
    set(gca,'linew',4)

    messageToPrint(nTile,1) = [panelLabels(nTile) + ": Optimum Smoothing Factor: " + num2str(maxIndX) + "; Optimum Ratio: " + num2str(maxIndY)];
    if (nTile == numPanels)
        annotation(gcf,'textbox',[0.51746299553432 0.10024154589372 0.463208685162847 0.17914653784219],...
    'String',[messageToPrint],...
    'FitBoxToText','on', 'FontSize', 26, 'LineWidth', 4);
    end
    
end

%% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG
