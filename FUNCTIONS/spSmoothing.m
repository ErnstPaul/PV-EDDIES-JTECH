%% SPATIAL SMOOTHING FUNCTION
function [smoothMat] = spSmoothing(roughMat, windowSize)

[x, y] = size(roughMat); %get dimensions
smoothMat = zeros(x,y); %create return matrix
mask = isnan(roughMat); %create NaN Mask

%Loop across all X's, all Y's
for i = 1:x
    for j = 1:y
        if isnan(roughMat(i,j))
            smoothMat(i,j) = NaN;
        end
        minBridgex = max([i-windowSize 1]);
        maxBridgex = min([i+windowSize x]);
        minBridgey = max([j-windowSize 1]);
        maxBridgey = min([j+windowSize y]);
        intermed = roughMat(minBridgex:maxBridgex, minBridgey:maxBridgey);
        %smooth
        smoothMat(i,j) = nanmean(nanmean(intermed));
    end
end

%mask back in good values so we don't screw it up
smoothMat(mask) = NaN;
end