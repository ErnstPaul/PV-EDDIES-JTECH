%% Paul_OPT_DisplayOptimizedCenterSmoothing.m (version 1.2)
%Author: Paul Ernst
%Date Created: 5/25/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 5/25/2022 - Created
%--------------------------------------
%Purpose: Displays the optimized smoothing of all center methods (Figure 2)
%Inputs: Plot sources from Comparison_Main
%Outputs: 5 panels with the following lines:
% -Center LNAM: SMTH    -Center VORT_T: SMTH     -Center VORT_WA: SMTH      -Center OW_T: SMTH
% -Center OW_WA: SMTH    -Center PV_ISO: SMTH
%% Inputs
titlestring = "OptimizedCenterSmoothing";
PlotSources = ["OPT_C_VORT_T_SMTH", "OPT_C_VORT_WA_SMTH", "OPT_C_OW_T_SMTH", "OPT_C_OW_WA_SMTH", "OPT_C_PV_ISO_SMTH"];
xLabels = ["Smoothing Factor", "Smoothing Factor", "Smoothing Factor", "Smoothing Factor", "Smoothing Factor"];

load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
%% Plots
numPanels = length(PlotSources);
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.53359375 0.917361111111111], 'Color', 'w'); %fullscreen figure
t = tiledlayout(3,2,'TileSpacing','tight');

panelLabels = ["(a) VORT_T","(b) VORT_W_A","(c) OW_T", "(d) OW_W_A", "(e) PV_I_S_O"];

%Loop thru plots
for nTile = 1:numPanels
    ax = nexttile();

    %load the source of this panel
    thisSource = ["IO_" + PlotSources(nTile) + "_Comparisons.mat"];
    load(thisSource)

    %Grab maximum for line later
    [MaxVal,MaxInd] = nanmax(nanmean(totalSims,1));

    %Plot individual lines
    plot(xToPlot,nanmean(totalNum,1),'--', 'linewidth', 4)
    hold on
    plot(xToPlot,nanmean(totalRatio,1),'--', 'linewidth', 4)
    plot(xToPlot,nanmean(totalSpatialC,1),'--', 'linewidth', 4)
    plot(xToPlot,nanmean(totalSpatialR,1), '--', 'linewidth', 4)
    plot(xToPlot,nanmean(totalSims,1)/100, 'linewidth', 6)
    plot([xToPlot(MaxInd) xToPlot(MaxInd)],[0 1], '--k', 'linewidth', 6)

    %labels and limits
    title(["Optimized at : " + num2str(xToPlot(MaxInd))])
    xlim([nanmin(xToPlot) nanmax(xToPlot)])
    ylim([0 1])
    ylabel("Scores")
    xlabel(xLabels(nTile))
    set(ax,'fontsize',28);
    grid on
    set(gca,'linew',4)

    xTxT = min(xToPlot)+0.02*(max(xToPlot)-min(xToPlot));
    yTxT = 0.88;
    text(xTxT, yTxT, panelLabels(nTile), 'FontSize', 28)

    if (nTile == numPanels)
        lgd = legend("Number Error","Ratio Error", "Spatial Correct Error", "Spatial Incorrect Error", "Aggregate Similarity Score", "Maximum Score");
        lgd.Layout.Tile = 6;
    end
    
end

%% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG
